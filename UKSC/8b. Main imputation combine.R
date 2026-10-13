###############################################################################
# 8b. Main imputation combine
#
# Purpose
#   Check every task file from File 8a and combine the imputations into one
#   multiply imputed dataset for Files 10 and 12. Run it on one computer that
#   holds every task file in tasks_dir/08a.
#
# Survey weights and the imputation model (a deliberate departure)
#   Silverwood et al. (2024) weighted their imputation models by the initial
#   design weights. That is not done here: R's mice package has no weighted
#   imputation methods for categorical variables. The log of the initial
#   design weight and the sampling stratum are predictors in every
#   imputation equation instead (File 7), the imputation models are
#   unweighted, and survey weights, PSUs and strata are applied when each
#   completed dataset is analysed (File 12). Weight interactions are not
#   included, which may overstate standard errors (Silverwood et al. 2024).
#   This script only combines the imputations and changes nothing in them.
#
# Checks (the script stops if any fails)
#   1. Every task in the table has a file.
#   2. Each file is the task it claims to be: the right task number and
#      seed, one imputation, and the number of iterations in File 7's
#      specification.
#   3. Each imputation used File 7's specification unchanged: the same data,
#      methods, predictor matrix, where matrix and visit sequence.
#   4. Only structural (not applicable) values are still missing in each
#      completed dataset. Every other missing value was imputed.
#   5. Every task was run with the same versions of R and mice, so results
#      from different computers can be combined.
#   The imputations are then combined with mice::ibind, which checks again
#   that they share the same data and model.
#   6. The combined object gives back every task's completed dataset
#      exactly (after repairing a quirk of ibind; Block 6).
#
# Reads
#   data_dir/07_imputation_spec.rds
#   tasks_dir/08a/imputation_NNN.rds
#
# Writes
#   data_dir/08_imputations_main.rds     the combined imputations (a mids
#                                        object; participant level)
#   review_dir/08_task_log.csv           per task: run time, warnings and
#                                        logged events counted, versions
#                                        (no data values)
#   review_dir/08_logged_events.csv      mice's logged events, all tasks
#                                        (variable and method names only)
#   review_dir/08_warnings.csv           distinct warnings, all tasks
###############################################################################

suppressPackageStartupMessages(library(mice))


###############################################################################
# Block 1. Paths and settings
#
# Folder locations are set once, in File 0, and read here from
# project_paths.rds. Run this script with the working directory set to the
# folder that holds the scripts.
###############################################################################

paths      <- readRDS("project_paths.rds")
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs
tasks_dir  <- paths$tasks_dir    # one file per task in split runs

base_seed   <- 8000              # as in File 8a: task t used seed base_seed + t
task_folder <- file.path(tasks_dir, "08a")


###############################################################################
# Block 2. Read File 7's specification and build the task table
###############################################################################

imputation_spec <- readRDS(file.path(data_dir, "07_imputation_spec.rds"))

task_table <- data.frame(task = seq_len(imputation_spec$production_m))
task_table$seed <- base_seed + task_table$task
task_table$file <- file.path(task_folder,
                             sprintf("imputation_%03d.rds", task_table$task))

# Which values should still be missing after imputation: those that are
# missing and not imputed, i.e. the structural (not applicable) values.
expected_missing <- is.na(imputation_spec$data) & !imputation_spec$where


###############################################################################
# Block 3. Check 1: every task file is present
###############################################################################

missing_tasks <- task_table$task[!file.exists(task_table$file)]
if (length(missing_tasks) > 0) {
  stop(length(missing_tasks), " task files are missing (tasks ",
       paste(missing_tasks, collapse = ", "), "). Run them with File 8a, ",
       "or copy their files into ", task_folder, ".")
}


###############################################################################
# Block 4. Checks 2-4, and combining, one task file at a time
#
# Each file is read and checked, then added to the combined imputations with
# mice::ibind. ibind joins two sets of imputations of the same data and
# model into one, and stops if they differ in their data, missing-data
# pattern or model. Adding the tasks in order gives the imputations in task
# order. A failed check stops the script and names the task, so that task's
# file can be deleted and rerun in File 8a. Nothing is saved until every
# check has passed.
###############################################################################

combined      <- NULL
task_log      <- data.frame()
logged_events <- data.frame()
all_warnings  <- character()

for (t in task_table$task) {
  result <- readRDS(task_table$file[t])
  imputation <- result$imputation

  # 2. The right task, seed, number of imputations and iterations.
  if (!identical(as.integer(result$task), as.integer(t)) ||
      !identical(as.numeric(result$seed), as.numeric(task_table$seed[t])) ||
      !is.mids(imputation) || imputation$m != 1 ||
      imputation$iteration != imputation_spec$iterations) {
    stop("Task ", t, "'s file is not the expected task (task number, seed, ",
         "one imputation or iterations differ). Delete it and rerun the task.")
  }

  # 3. File 7's specification, unchanged. mice stores the predictor matrix
  # with a different number type and drops the row names of the where
  # matrix, so those two are compared by value.
  if (!identical(imputation$data, imputation_spec$data) ||
      !identical(imputation$method, imputation_spec$method) ||
      !isTRUE(all.equal(imputation$predictorMatrix, imputation_spec$predictor_matrix)) ||
      !identical(colnames(imputation$where), colnames(imputation_spec$where)) ||
      any(imputation$where != imputation_spec$where) ||
      !identical(imputation$visitSequence, imputation_spec$visit_sequence)) {
    stop("Task ", t, " was not run with File 7's current specification. ",
         "Delete its file and rerun the task with the current specification.")
  }

  # 4. Only structural values still missing.
  completed <- complete(imputation, 1)
  if (any(is.na(completed) != expected_missing)) {
    stop("Task ", t, "'s completed data has missing values other than the ",
         "structural (not applicable) ones. Delete its file and rerun the task.")
  }

  # Add the imputation to the combined set.
  if (is.null(combined)) {
    combined <- imputation
  } else {
    combined <- ibind(combined, imputation)
  }

  # Record the task (no data values).
  events <- imputation$loggedEvents
  n_events <- if (is.null(events)) 0L else nrow(events)
  if (n_events > 0) {
    events$task <- t
    logged_events <- rbind(logged_events, events)
  }
  all_warnings <- c(all_warnings, result$warnings)

  task_log <- rbind(task_log, data.frame(
    task          = t,
    seed          = result$seed,
    computer      = result$computer,
    minutes       = round(result$minutes, 1),
    warnings      = length(result$warnings),
    logged_events = n_events,
    r_version     = result$r_version,
    mice_version  = result$mice_version
  ))
}


###############################################################################
# Block 5. Check 5: the same R and mice versions for every task
###############################################################################

if (length(unique(task_log$r_version)) > 1 ||
    length(unique(task_log$mice_version)) > 1) {
  stop("Tasks were run with different versions of R (",
       paste(unique(task_log$r_version), collapse = ", "), ") or mice (",
       paste(unique(task_log$mice_version), collapse = ", "), "). Align the ",
       "versions (File 0 records them on each computer) and rerun the tasks ",
       "that used the other versions.")
}


###############################################################################
# Block 6. Repair the combined object, check it, and save it
#
# mice keeps the imputed values in a list with one entry per variable,
# empty for variables with nothing imputed (e.g. the design predictors).
# In mice 3.19, ibind drops those empty entries (in R, assigning an empty
# entry to a list deletes it), and mice::complete() then matches the
# remaining entries to the wrong variables by position. The full list is
# restored here, in the data's column order.
# Check 6 then confirms that the combined object gives back every task's
# completed dataset exactly, whatever mice does internally.
###############################################################################

if (combined$m != nrow(task_table)) {
  stop("The combined imputations number ", combined$m, ", not ",
       nrow(task_table), ".")
}

full_imp <- vector("list", ncol(combined$data))
names(full_imp) <- names(combined$data)
for (variable in names(combined$imp)) {
  full_imp[[variable]] <- combined$imp[[variable]]
}
combined$imp <- full_imp

# 6. Each completed dataset from the combined object equals that task's own.
for (t in task_table$task) {
  task_imputation <- readRDS(task_table$file[t])$imputation
  if (!identical(complete(combined, t), complete(task_imputation, 1))) {
    stop("Completed dataset ", t, " from the combined imputations differs ",
         "from task ", t, "'s own. The combined object is not usable.")
  }
}

saveRDS(combined, file.path(data_dir, "08_imputations_main.rds"))


###############################################################################
# Block 7. Review tables
#
# Disclosure check: these tables hold task numbers, seeds, computer names,
# run times, version numbers, counts of warnings and logged events per task,
# variable and method names, and warning messages. None holds a count of
# people or a data value.
###############################################################################

task_log$disclosure_checked <- TRUE
write.csv(task_log, file.path(review_dir, "08_task_log.csv"), row.names = FALSE)

if (nrow(logged_events) == 0) {
  logged_events <- data.frame(it = integer(), im = integer(), dep = character(),
                              meth = character(), out = character(),
                              task = integer())
}
logged_events$disclosure_checked <- rep(TRUE, nrow(logged_events))
write.csv(logged_events, file.path(review_dir, "08_logged_events.csv"),
          row.names = FALSE)

warning_table <- data.frame(message = unique(all_warnings))
warning_table$times <- as.integer(table(all_warnings)[warning_table$message])
warning_table$disclosure_checked <- rep(TRUE, nrow(warning_table))
write.csv(warning_table, file.path(review_dir, "08_warnings.csv"), row.names = FALSE)


###############################################################################
# Block 8. Summary
#
# Counts of tasks, events and versions only. No data values.
###############################################################################

cat("\nFile 8b finished.\n")
cat("  Imputations combined:                 ", combined$m, "\n")
cat("  Iterations per imputation:            ", combined$iteration, "\n")
cat("  Checks (files, specification, missing values, versions,\n")
cat("    combined equals each task): all passed\n")
cat("  R / mice versions:                    ", task_log$r_version[1], "/",
    task_log$mice_version[1], "\n")
cat("  Computers used:                       ", length(unique(task_log$computer)), "\n")
cat("  Total run time (hours, all tasks):    ", round(sum(task_log$minutes) / 60, 1), "\n")
cat("  Logged events (all tasks):            ", nrow(logged_events), "\n")
cat("  Distinct warnings (all tasks):        ", nrow(warning_table), "\n")
cat("  Saved:", file.path(data_dir, "08_imputations_main.rds"), "\n")
cat("  Next: File 10 (diagnostics on the main imputations) and File 12\n")
cat("    (pooled analysis, which re-checks the number of imputations\n")
cat("    against the fraction of missing information).\n")
