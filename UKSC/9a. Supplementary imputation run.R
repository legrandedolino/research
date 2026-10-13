###############################################################################
# 9a. Supplementary imputation run
#
# Purpose
#   Run the imputation model of Files 7 and 8a again without the CLS
#   predictors of non-response (the user guide's Appendix C variables chosen
#   in File 6). Everything else is identical. Comparing the results with the
#   main imputations (File 12) shows what those predictors add, as in the
#   supplementary analysis of Silverwood et al. (2024), who imputed with and
#   without the predictors of non-response they had identified.
#   Each imputation is one independent task with its own seed and output
#   file, so the run can be split across several computers (or several R
#   sessions on one computer) and resumed after an interruption. File 9b
#   checks the task files and combines them.
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
#   File 7's specification, minus the predictors of non-response:
#     - removed: every auxiliary File 6 kept from the CLS predictors of
#       non-response. Their columns leave the data, and their rows and
#       columns leave the predictor matrix;
#     - kept unchanged: the design predictors, the two prespecified
#       auxiliaries (SRET job and main activity at 25, included for their
#       relation to the outcome rather than to non-response), every analysis
#       variable, the routing (where matrix), the passive calculations, the
#       methods, the visiting order of the remaining variables, the number of
#       iterations and the number of imputations.
#   Removing predictors can't make any predictor category sparser, so File
#   7's check that every category has enough people for the bootstrap
#   logistic models still holds.
#
# Running on several computers or R sessions
#   As for File 8a: run Files 0-7 once, copy data_dir/07_imputation_spec.rds
#   and data_dir/06_auxiliary_list.rds to any other computer that will run
#   tasks (only between computers covered by the UKDS End User Licence and
#   your institution's data security rules), and give each computer or R
#   session its own tasks_to_run (Block 1). mice uses one processor core,
#   so a computer with several cores can run several R sessions at once.
#   Copy every task file into tasks_dir/09a on one computer and run File 9b
#   there. A task whose file already exists is skipped, and each task's
#   seed is fixed (base seed + task number).
#
# Reads
#   data_dir/07_imputation_spec.rds
#   data_dir/06_auxiliary_list.rds
#
# Writes
#   data_dir/09_supplementary_spec.rds  the supplementary specification
#       (the same on every computer; File 9b checks the tasks against it)
#   tasks_dir/09a/imputation_NNN.rds    one per task: the imputation (a mids
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
# e.g. 1:30. Every task in the table is needed before File 9b can combine.
tasks_to_run <- "all"

base_seed <- 9000   # task t uses seed base_seed + t (File 8a used 8000)

task_folder <- file.path(tasks_dir, "09a")
dir.create(task_folder, recursive = TRUE, showWarnings = FALSE)


###############################################################################
# Block 2. Read File 7's specification and File 6's auxiliary screen
###############################################################################

imputation_spec <- readRDS(file.path(data_dir, "07_imputation_spec.rds"))
file6           <- readRDS(file.path(data_dir, "06_auxiliary_list.rds"))


###############################################################################
# Block 3. The supplementary specification
#
# The predictors of non-response to remove are the auxiliaries File 6 kept
# from the group "CLS predictor of non-response".
###############################################################################

nonresponse_predictors <- file6$screen$variable[
  file6$screen$group == "CLS predictor of non-response"
]
to_remove <- intersect(file6$auxiliary_list, nonresponse_predictors)
to_keep   <- setdiff(names(imputation_spec$data), to_remove)

supplementary_spec <- list(
  data             = imputation_spec$data[, to_keep],
  design_data      = imputation_spec$design_data,
  variables        = imputation_spec$variables[
                       imputation_spec$variables$variable %in% to_keep, ],
  method           = imputation_spec$method[to_keep],
  predictor_matrix = imputation_spec$predictor_matrix[to_keep, to_keep],
  where            = imputation_spec$where[, to_keep],
  visit_sequence   = imputation_spec$visit_sequence[
                       imputation_spec$visit_sequence %in% to_keep],
  iterations       = imputation_spec$iterations,
  production_m     = imputation_spec$production_m,
  removed          = to_remove
)

# Checks: the removed variables are all in File 7's model and are all
# predictors of non-response; no analysis variable or prespecified
# auxiliary was removed; no passive calculation uses a removed variable.
if (length(to_remove) == 0) {
  stop("File 6 kept no predictors of non-response, so the supplementary ",
       "imputation would repeat the main one.")
}
if (!all(to_remove %in% names(imputation_spec$data))) {
  stop("The predictors of non-response to remove aren't all in File 7's model.")
}
analysis_and_prespecified <- c(
  imputation_spec$variables$variable[
    !imputation_spec$variables$role %in% c("Auxiliary (File 6)")],
  file6$screen$variable[file6$screen$group == "Prespecified Category b"]
)
if (any(to_remove %in% analysis_and_prespecified)) {
  stop("An analysis variable or prespecified auxiliary would be removed.")
}
passive <- startsWith(supplementary_spec$method, "~")
for (formula_text in supplementary_spec$method[passive]) {
  if (any(all.vars(as.formula(formula_text)) %in% to_remove)) {
    stop("A passive calculation uses a removed predictor of non-response.")
  }
}
if (!identical(names(supplementary_spec$method), colnames(supplementary_spec$predictor_matrix)) ||
    !identical(names(supplementary_spec$data), colnames(supplementary_spec$where))) {
  stop("The supplementary specification's pieces don't line up.")
}

# Several R sessions may start at once and all save this (identical) file,
# so each writes under its own temporary name and then renames it.
spec_file      <- file.path(data_dir, "09_supplementary_spec.rds")
temporary_spec <- paste0(spec_file, ".", Sys.getpid(), ".partial")
saveRDS(supplementary_spec, temporary_spec)
invisible(file.rename(temporary_spec, spec_file))


###############################################################################
# Block 4. Task table
#
# One row per imputation: its task number, seed and output file.
###############################################################################

task_table <- data.frame(task = seq_len(supplementary_spec$production_m))
task_table$seed <- base_seed + task_table$task
task_table$file <- file.path(task_folder,
                             sprintf("imputation_%03d.rds", task_table$task))

if (identical(tasks_to_run, "all")) tasks_to_run <- task_table$task
if (!all(tasks_to_run %in% task_table$task)) {
  stop("tasks_to_run includes task numbers outside 1-", nrow(task_table), ".")
}


###############################################################################
# Block 5. Workspace check
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
# Block 6. Run the tasks
#
# Each task is independent: it imputes once, from its own seed, and saves
# its own file. The file is first written under a temporary name and then
# renamed, so an interruption during saving never leaves a half-written file
# that a later run would skip. Warnings are kept with the task, not printed.
# Progress lines show task numbers and times only.
###############################################################################

cat("\nSupplementary imputation run (", length(to_remove),
    " predictors of non-response removed): ", length(tasks_to_run), " of ",
    nrow(task_table), " tasks on this computer; ", supplementary_spec$iterations,
    " iterations each.\n", sep = "")

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
    mice(supplementary_spec$data,
         m               = 1,
         maxit           = supplementary_spec$iterations,
         method          = supplementary_spec$method,
         predictorMatrix = supplementary_spec$predictor_matrix,
         where           = supplementary_spec$where,
         visitSequence   = supplementary_spec$visit_sequence,
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
# Block 7. Summary
#
# Task counts only. No data values.
###############################################################################

tasks_done <- sum(file.exists(task_table$file))
cat("\nFile 9a finished on this computer.\n")
cat("  Predictors of non-response removed:", length(to_remove), "\n")
cat("  Variables in the supplementary model:", ncol(supplementary_spec$data), "\n")
cat("  Tasks in the table:        ", nrow(task_table), "\n")
cat("  Task files now present:    ", tasks_done, "\n")
if (tasks_done < nrow(task_table)) {
  cat("  Still to run (on this or other computers):",
      nrow(task_table) - tasks_done, "tasks\n")
} else {
  cat("  All tasks done. Next: File 9b (combine), on one computer holding",
      "every task file.\n")
}
