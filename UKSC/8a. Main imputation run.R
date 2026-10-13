###############################################################################
# 8a. Main imputation run
#
# Purpose
#   Run the multiple imputation model built and tested in File 7 for the
#   full target population. Each imputation is one independent task with
#   its own seed and its own output file, so the run can be split across
#   several computers (or several R sessions on one computer) and resumed
#   after an interruption. File 8b checks the task files and combines them.
#
# Survey weights and the imputation model (a deliberate departure)
#   Silverwood et al. (2024) weighted their imputation models by the initial
#   design weights. That is not done here. The imputation uses R's mice
#   package, which has no weighted imputation methods for categorical
#   variables, and most of the variables imputed here are categorical.
#   Instead:
#     - the log of the initial design weight and the sampling stratum are
#       included as predictors in every imputation equation (File 7). This
#       is a standard way of carrying survey design information into an
#       imputation model (Reiter, Raghunathan and Kinney 2006). As a
#       predictor, strata too small for the imputation models are pooled
#       into one group; the analysis uses the original strata;
#     - the imputation models themselves are unweighted. (The CLS user
#       guide's worked example is from a cohort without sampling weights,
#       so it doesn't settle this either way.) With the design variables
#       as predictors, weighting the imputation equations as well would
#       mainly matter if relationships between variables differ by
#       selection probability, i.e. weight interactions (below);
#     - survey weights, PSUs and strata are applied when each completed
#       dataset is analysed (File 12).
#   Interactions between the weight and other variables (Seaman et al. 2012)
#   are not included. Silverwood et al. could not include them either,
#   because their models became unstable. They note that this should not
#   affect point estimates but may overstate standard errors, which makes
#   conclusions slightly conservative. The same applies here.
#
# What is run
#   Everything comes from File 7's specification: the data, the imputation
#   method for each variable, the predictor matrix, which values are imputed
#   (routed questions are imputed only where they apply), the order in which
#   variables are visited, the number of iterations and the number of
#   imputations. The number of imputations follows the CLS user guide: at
#   least 50, and at least the largest fraction of missing information (as
#   a percentage) estimated in File 7's test run. File 12 checks that rule
#   again on the main imputations.
#   Nothing in the model is changed here.
#
# Running on several computers or R sessions
#   1. Run Files 0-7 once. Copy data_dir/07_imputation_spec.rds to every
#      other computer that will run tasks, into the same place under its
#      own output folder (after running File 0 there). Only copy data
#      between computers covered by the UKDS End User Licence and your
#      institution's data security rules.
#   2. On each computer or R session, set tasks_to_run (Block 1) to a
#      different set of task numbers, e.g. 1:30, 31:60 and 61:90.
#      mice uses one processor core, so a computer with several cores can
#      run several R sessions at once, each with its own tasks_to_run.
#   3. Copy every task file into tasks_dir/08a on one computer and run
#      File 8b there.
#   A task whose file already exists is skipped, so an interrupted run
#   resumes where it stopped. Each task's seed is fixed (base seed + task
#   number), so the results don't depend on which computer ran which task,
#   or in what order. File 0 records package versions on each computer, and
#   File 8b combines tasks only if they were all run with the same versions
#   of R and mice.
#
# Reads
#   data_dir/07_imputation_spec.rds
#
# Writes
#   tasks_dir/08a/imputation_NNN.rds  one per task: the imputation (a mids
#       object with one imputation), plus its task number, seed, run time,
#       warnings, and the R and mice versions (participant level; never
#       shared)
###############################################################################

suppressPackageStartupMessages(library(mice))


###############################################################################
# Block 1. Paths and settings
#
# Folder locations are set once, in File 0, and read here from
# project_paths.rds. Run this script with the working directory set to the
# folder that holds the scripts.
###############################################################################

paths     <- readRDS("project_paths.rds")
data_dir  <- paths$data_dir     # datasets passed between scripts
tasks_dir <- paths$tasks_dir    # one file per task in split runs

# Which tasks this computer (or R session) runs: "all", or task numbers,
# e.g. 1:30. Every task in the table is needed before File 8b can combine.
tasks_to_run <- "all"

base_seed <- 8000   # task t uses seed base_seed + t

task_folder <- file.path(tasks_dir, "08a")
dir.create(task_folder, recursive = TRUE, showWarnings = FALSE)


###############################################################################
# Block 2. Read File 7's specification
###############################################################################

imputation_spec <- readRDS(file.path(data_dir, "07_imputation_spec.rds"))


###############################################################################
# Block 3. Task table
#
# One row per imputation: its task number, seed and output file.
###############################################################################

task_table <- data.frame(task = seq_len(imputation_spec$production_m))
task_table$seed <- base_seed + task_table$task
task_table$file <- file.path(task_folder,
                             sprintf("imputation_%03d.rds", task_table$task))

if (identical(tasks_to_run, "all")) tasks_to_run <- task_table$task
if (!all(tasks_to_run %in% task_table$task)) {
  stop("tasks_to_run includes task numbers outside 1-", nrow(task_table), ".")
}


###############################################################################
# Block 4. Workspace check
#
# mice finds two of its internal objects, "state" and "loggedEvents", by
# name. If the R workspace holds an object with either name (e.g. from other
# code run in the same session), mice picks it up instead and fails with
# "$ operator is invalid for atomic vectors" the first time it logs an
# event. Stop now with a clear message instead.
###############################################################################

for (name in c("state", "loggedEvents")) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) {
    stop("An object called '", name, "' is in the R workspace, and mice would ",
         "mistake it for its own. Remove it with rm(", name, ") or restart R, ",
         "then run this script again.")
  }
}


###############################################################################
# Block 5. Run the tasks
#
# Each task is independent: it imputes once, from its own seed, and saves
# its own file. The file is first written under a temporary name and then
# renamed, so an interruption during saving never leaves a half-written file
# that a later run would skip. Warnings are kept with the task, not printed.
# Progress lines show task numbers and times only.
###############################################################################

cat("\nMain imputation run:", length(tasks_to_run), "of", nrow(task_table),
    "tasks on this computer;", imputation_spec$iterations, "iterations each.\n")

for (t in tasks_to_run) {
  task_file <- task_table$file[t]

  # Skip tasks that are already done.
  if (file.exists(task_file)) {
    cat("  Task", t, "already done; skipped.\n")
    next
  }

  cat("  Task", t, "started", format(Sys.time(), "%H:%M"), "\n")
  start_time <- Sys.time()

  task_warnings <- character()
  imputation <- withCallingHandlers(
    mice(imputation_spec$data,
         m               = 1,
         maxit           = imputation_spec$iterations,
         method          = imputation_spec$method,
         predictorMatrix = imputation_spec$predictor_matrix,
         where           = imputation_spec$where,
         visitSequence   = imputation_spec$visit_sequence,
         seed            = task_table$seed[t],
         printFlag       = FALSE),
    warning = function(w) {
      # Each time a warning occurs, keep its message and let the imputation
      # carry on ("<<-" adds it to task_warnings, outside this small function).
      task_warnings <<- c(task_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  task_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  task_result <- list(
    task         = t,
    seed         = task_table$seed[t],
    imputation   = imputation,
    minutes      = task_minutes,
    warnings     = task_warnings,
    r_version    = as.character(getRversion()),
    mice_version = as.character(packageVersion("mice")),
    computer     = gsub("[^A-Za-z0-9_-]", "_", Sys.info()[["nodename"]])
  )

  temporary_file <- paste0(task_file, ".partial")
  saveRDS(task_result, temporary_file)
  file.rename(temporary_file, task_file)

  cat("  Task", t, "finished in", round(task_minutes, 1), "minutes",
      if (length(task_warnings) > 0) paste0("(", length(task_warnings), " warnings kept)") else "",
      "\n")
}


###############################################################################
# Block 6. Summary
#
# Task counts only. No data values.
###############################################################################

tasks_done <- sum(file.exists(task_table$file))
cat("\nFile 8a finished on this computer.\n")
cat("  Tasks in the table:        ", nrow(task_table), "\n")
cat("  Task files now present:    ", tasks_done, "\n")
if (tasks_done < nrow(task_table)) {
  cat("  Still to run (on this or other computers):",
      nrow(task_table) - tasks_done, "tasks\n")
} else {
  cat("  All tasks done. Next: File 8b (combine), on one computer holding",
      "every task file.\n")
}
