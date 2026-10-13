###############################################################################
# 17b. Repeated nested cross-validation combine
#
# Purpose
#   Check every task file of File 17a and combine them for Files 18 and 19:
#   the out-of-fold predictions, the permuted predictions for importance, and
#   the technical diagnostics. No performance measure is computed here.
#
# Checks (the script stops if any fails)
#   1. Every task in the table has a file, and each file is the task it
#      claims to be (strategy, repeat, fold, landmark).
#   2. Each task's predictions cover exactly its outer test fold, once each,
#      and are complete and between 0 and 1.
#   3. For every strategy x repeat x landmark, the folds together predict
#      every member of the prediction cohort exactly once.
#   4. Every task was run with the same versions of R, mice and glmnet.
#
# Diagnostics summarized across tasks (review tables; no people counted)
#   - run times and versions per task;
#   - imputation convergence: for each strategy x landmark x variable, the
#     share of tasks flagged and the median split R-hat;
#   - the classification-tree fallback: how often each model switched;
#   - the auxiliary screen: how often each auxiliary was kept;
#   - model checks: chosen alpha, lambda at an edge, non-zero coefficients,
#     and the logistic benchmark's convergence and warnings.
#
# Reads
#   data_dir/14_prediction_design.rds
#   tasks_dir/17a/task_NNN.rds
#
# Writes
#   data_dir/17_cv_results.rds         predictions, permuted predictions and
#                                      diagnostics (participant level)
#   review_dir/17_task_log.csv, 17_convergence_summary.csv,
#     17_fallback_summary.csv, 17_auxiliary_summary.csv,
#     17_model_summary.csv, 17_warnings.csv
###############################################################################


###############################################################################
# Block 1. Paths and settings
###############################################################################

paths      <- readRDS("project_paths.rds")
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs
tasks_dir  <- paths$tasks_dir    # one file per task in split runs

task_folder <- file.path(tasks_dir, "17a")


###############################################################################
# Block 2. The design and the task table (as File 17a)
###############################################################################

design   <- readRDS(file.path(data_dir, "14_prediction_design.rds"))
cohort   <- design$cohort
settings <- design$settings

task_table <- rbind(
  expand.grid(landmark = settings$landmark_order, outer_fold = seq_len(settings$n_outer_folds),
              repeat_number = seq_len(settings$n_repeats), strategy = "primary",
              stringsAsFactors = FALSE),
  expand.grid(landmark = settings$landmark_order, outer_fold = seq_len(settings$n_outer_folds),
              repeat_number = 1, strategy = "outcome omitted",
              stringsAsFactors = FALSE)
)
task_table <- task_table[, c("strategy", "repeat_number", "outer_fold", "landmark")]
task_table$task <- seq_len(nrow(task_table))
task_table$file <- file.path(task_folder, sprintf("task_%03d.rds", task_table$task))


###############################################################################
# Block 3. Checks 1, 2 and 4, task by task, and collecting the results
###############################################################################

missing_tasks <- task_table$task[!file.exists(task_table$file)]
if (length(missing_tasks) > 0) {
  stop(length(missing_tasks), " task files are missing (tasks ",
       paste(missing_tasks, collapse = ", "), "). Run them with File 17a.")
}

predictions  <- list()
permuted     <- list()
task_log     <- data.frame()
convergence  <- data.frame()
switches     <- data.frame()
screen       <- data.frame()
model_checks <- data.frame()
warnings_all <- data.frame()

for (t in task_table$task) {
  result <- readRDS(task_table$file[t])

  # 1. The right task.
  if (result$task != t || result$strategy != task_table$strategy[t] ||
      result$repeat_number != task_table$repeat_number[t] ||
      result$outer_fold != task_table$outer_fold[t] ||
      result$landmark != task_table$landmark[t]) {
    stop("Task file ", t, " doesn't match the task table.")
  }

  # 2. Predictions cover exactly the outer test fold, once each.
  test_ids <- cohort$NSID[cohort[[paste0("outer_fold_r", result$repeat_number)]] == result$outer_fold]
  p <- result$predictions
  if (!setequal(p$NSID, test_ids) || anyDuplicated(p$NSID) > 0 || nrow(p) != length(test_ids)) {
    stop("Task ", t, ": predictions don't cover exactly its test fold.")
  }
  probabilities <- as.matrix(p[, c("p_elastic_net", "p_logistic", "p_null")])
  if (anyNA(probabilities) || any(probabilities < 0 | probabilities > 1) ||
      anyNA(result$permuted_predictions)) {
    stop("Task ", t, ": predictions are missing or outside 0-1.")
  }

  predictions[[t]] <- p
  permuted[[t]] <- list(NSID = p$NSID, units = result$permutation_units,
                        predictions = result$permuted_predictions)

  task_label <- data.frame(task = t, strategy = result$strategy,
                           repeat_number = result$repeat_number,
                           outer_fold = result$outer_fold, landmark = result$landmark)
  task_log <- rbind(task_log, data.frame(
    task_label,
    imputation_minutes = round(result$imputation_minutes, 1),
    total_minutes = round(result$total_minutes, 1),
    fallback_models = if (nrow(result$switches) > 0) length(unique(result$switches$model_for)) else 0L,
    auxiliaries_kept = sum(result$screen$kept),
    convergence_flags = sum(result$convergence$flag),
    convergence_rows = nrow(result$convergence),
    logged_events = result$logged_events,
    warnings = sum(result$warnings),
    r_version = result$r_version, mice_version = result$mice_version,
    glmnet_version = result$glmnet_version, computer = result$computer
  ))
  if (nrow(result$convergence) > 0) convergence <- rbind(convergence, data.frame(task_label, result$convergence))
  if (nrow(result$switches) > 0) switches <- rbind(switches, data.frame(task_label, result$switches))
  screen <- rbind(screen, data.frame(task_label, result$screen))
  model_checks <- rbind(model_checks, data.frame(task_label, result$model_checks))
  if (length(result$warnings) > 0) {
    warnings_all <- rbind(warnings_all, data.frame(
      task_label, message = names(result$warnings), times = as.integer(result$warnings)))
  }
}

# 4. The same versions for every task.
if (length(unique(task_log$r_version)) > 1 || length(unique(task_log$mice_version)) > 1 ||
    length(unique(task_log$glmnet_version)) > 1) {
  stop("Tasks were run with different versions of R, mice or glmnet (",
       paste(unique(paste(task_log$r_version, task_log$mice_version, task_log$glmnet_version)),
             collapse = "; "), ").")
}

# Join the tasks' predictions into one table, in task order (do.call passes
# every table in the list to rbind at once).
predictions <- do.call(rbind, predictions)


###############################################################################
# Block 4. Check 3: each strategy x repeat x landmark predicts everyone once
###############################################################################

combinations <- unique(predictions[, c("strategy", "repeat_number", "landmark")])
for (k in seq_len(nrow(combinations))) {
  rows <- predictions$strategy == combinations$strategy[k] &
    predictions$repeat_number == combinations$repeat_number[k] &
    predictions$landmark == combinations$landmark[k]
  if (!setequal(predictions$NSID[rows], cohort$NSID) || anyDuplicated(predictions$NSID[rows]) > 0) {
    stop("Out-of-fold predictions don't cover the cohort exactly once (",
         paste(combinations[k, ], collapse = ", "), ").")
  }
}

saveRDS(list(task_table = task_table, predictions = predictions, permuted = permuted,
             convergence = convergence, switches = switches, screen = screen,
             model_checks = model_checks, task_log = task_log),
        file.path(data_dir, "17_cv_results.rds"))


###############################################################################
# Block 5. Review tables
#
# Disclosure check: run times, versions, computer names, counts of models,
# auxiliaries, flags and warnings (not people), shares of tasks, convergence
# statistics of imputed values (each from at least 10 imputed values), model
# settings and warning messages. No prediction or performance is written.
###############################################################################

task_log$disclosure_checked <- TRUE
write.csv(task_log, file.path(review_dir, "17_task_log.csv"), row.names = FALSE)

# Convergence: share of tasks flagged, and median split R-hat.
convergence_summary <- data.frame()
if (nrow(convergence) > 0) {
  key <- paste(convergence$strategy, convergence$landmark, convergence$variable,
               convergence$statistic, sep = " | ")
  for (k in unique(key)) {
    rows <- convergence[key == k, ]
    convergence_summary <- rbind(convergence_summary, data.frame(
      strategy = rows$strategy[1], landmark = rows$landmark[1], variable = rows$variable[1],
      statistic = rows$statistic[1], tasks = nrow(rows),
      share_flagged = round(mean(rows$flag), 2),
      median_split_rhat = round(median(rows$split_rhat, na.rm = TRUE), 3),
      share_trend_beyond_spread = round(mean(rows$trend_p < 0.01 & rows$drift_vs_spread > 1,
                                             na.rm = TRUE), 2)))
  }
}
convergence_summary$disclosure_checked <- rep(TRUE, nrow(convergence_summary))
write.csv(convergence_summary, file.path(review_dir, "17_convergence_summary.csv"), row.names = FALSE)

# Classification-tree fallback: tasks in which each model switched.
fallback_summary <- data.frame()
if (nrow(switches) > 0) {
  switched <- unique(switches[, c("task", "strategy", "landmark", "model_for")])
  key <- paste(switched$strategy, switched$landmark, switched$model_for, sep = " | ")
  for (k in unique(key)) {
    rows <- switched[key == k, ]
    fallback_summary <- rbind(fallback_summary, data.frame(
      strategy = rows$strategy[1], landmark = rows$landmark[1], model_for = rows$model_for[1],
      tasks_switched = nrow(rows),
      tasks_at_landmark = sum(task_table$strategy == rows$strategy[1] &
                                task_table$landmark == rows$landmark[1]),
      sparse_predictors = paste(unique(switches$predictor[paste(switches$strategy, switches$landmark,
                                                                 switches$model_for, sep = " | ") == k]),
                                collapse = " | ")))
  }
}
fallback_summary$disclosure_checked <- rep(TRUE, nrow(fallback_summary))
write.csv(fallback_summary, file.path(review_dir, "17_fallback_summary.csv"), row.names = FALSE)

# Auxiliary screen: share of tasks in which each auxiliary was kept.
auxiliary_summary <- data.frame()
key <- paste(screen$strategy, screen$landmark, screen$auxiliary, sep = " | ")
for (k in unique(key)) {
  rows <- screen[key == k, ]
  auxiliary_summary <- rbind(auxiliary_summary, data.frame(
    strategy = rows$strategy[1], landmark = rows$landmark[1], auxiliary = rows$auxiliary[1],
    screened_against = paste(unique(rows$screened_against), collapse = " | "),
    tasks = nrow(rows), share_kept = round(mean(rows$kept), 2)))
}
auxiliary_summary$disclosure_checked <- TRUE
write.csv(auxiliary_summary, file.path(review_dir, "17_auxiliary_summary.csv"), row.names = FALSE)

# Model checks by strategy x landmark.
model_summary <- data.frame()
key <- paste(model_checks$strategy, model_checks$landmark, sep = " | ")
for (k in unique(key)) {
  rows <- model_checks[key == k, ]
  model_summary <- rbind(model_summary, data.frame(
    strategy = rows$strategy[1], landmark = rows$landmark[1], fits = nrow(rows),
    alpha_0 = sum(rows$chosen_alpha == 0), alpha_0.25 = sum(rows$chosen_alpha == 0.25),
    alpha_0.5 = sum(rows$chosen_alpha == 0.5), alpha_0.75 = sum(rows$chosen_alpha == 0.75),
    alpha_1 = sum(rows$chosen_alpha == 1),
    lambda_at_largest = sum(rows$lambda_at_largest),
    lambda_at_smallest = sum(rows$lambda_at_smallest),
    median_nonzero_coefficients = median(rows$nonzero_coefficients),
    median_columns = median(rows$columns),
    logistic_not_converged = sum(!rows$logistic_converged),
    logistic_with_warnings = sum(rows$logistic_warnings > 0)))
}
model_summary$disclosure_checked <- TRUE
write.csv(model_summary, file.path(review_dir, "17_model_summary.csv"), row.names = FALSE)

# Warnings, totalled by message.
warning_summary <- data.frame(message = character(), tasks = integer(), times = integer())
if (nrow(warnings_all) > 0) {
  for (m in unique(warnings_all$message)) {
    rows <- warnings_all[warnings_all$message == m, ]
    warning_summary <- rbind(warning_summary, data.frame(
      message = m, tasks = length(unique(rows$task)), times = sum(rows$times)))
  }
}
warning_summary$disclosure_checked <- rep(TRUE, nrow(warning_summary))
write.csv(warning_summary, file.path(review_dir, "17_warnings.csv"), row.names = FALSE)


###############################################################################
# Block 6. Summary
#
# Counts only. No predictions or performance.
###############################################################################

cat("\nFile 17b finished.\n")
cat("  Tasks combined:                         ", nrow(task_table), "\n")
cat("  Checks (task identity, test-fold coverage, cohort coverage, versions): all passed\n")
cat("  R / mice / glmnet versions:             ", task_log$r_version[1], "/",
    task_log$mice_version[1], "/", task_log$glmnet_version[1], "\n")
cat("  Total run time (hours):                 ", round(sum(task_log$total_minutes) / 60, 1), "\n")
cat("  Tasks with a classification-tree fallback:", sum(task_log$fallback_models > 0), "\n")
cat("  Convergence rows flagged:               ", sum(task_log$convergence_flags), "of",
    sum(task_log$convergence_rows), "\n")
cat("  Review: 17_task_log.csv, 17_convergence_summary.csv, 17_fallback_summary.csv,\n")
cat("    17_auxiliary_summary.csv, 17_model_summary.csv, 17_warnings.csv\n")
