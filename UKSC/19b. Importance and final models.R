###############################################################################
# 19b. Importance and final models
#
# Purpose
#   1. Permutation importance of the elastic net (the primary model). In
#      each outer test fold, File 17a shuffled each feature, and each block
#      of features, among the test rows and predicted again. Here the fall
#      in performance is scored with File 18's measures, in each of the 25
#      primary outer splits (5 repeats x 5 folds) of every landmark, and
#      summarized over the splits: the mean fall, its standard deviation
#      between splits, and the number of splits in which performance fell.
#   2. The final models of File 19a: each landmark's elastic-net
#      coefficients, averaged over the imputed datasets, with the tuning
#      and imputation diagnostics.
#
# Interpretation
#   Importance describes what the models relied on among correlated
#   features, and the elastic net shrinks coefficients on purpose. Neither
#   estimates an association or an effect of a characteristic on careers.
#   A feature that is shuffled can be partly replaced by features related
#   to it, so block importance is reported beside feature importance.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   - Importance: scored only in outer test folds with at least 10 events
#     and 10 non-events (checked); only summaries over splits are written.
#   - Final models: shown only with at least 10 events, 10 non-events and
#     10 residual degrees of freedom. A coefficient for a 0/1 term (a
#     category, or a two-category feature) is shown only if, in every
#     imputed dataset, at least 10 people have the term, with at least 10
#     events and 10 non-events among them. The counts stay in a private
#     table.
#
# Reads
#   data_dir/14_prediction_design.rds
#   data_dir/17_cv_results.rds
#   data_dir/18_performance.rds (to check the measures are File 18's)
#   tasks_dir/19a/task_N.rds
#
# Writes
#   data_dir/19_final_models.rds             importance by split, summaries,
#                                            coefficients and diagnostics
#   output_dir/private_19_term_counts.csv    counts behind the coefficient check
#   review_dir/19_importance.csv, 19_importance_plots.pdf,
#     19_final_model_coefficients.csv, 19_final_model_tuning.csv,
#     19_disclosure_log.csv
#
# References
#   Collins GS, et al. (2024). TRIPOD+AI statement. BMJ 385: e078378.
#   Kapoor S, et al. (2024). REFORMS. Science Advances 10: eadk3452.
###############################################################################


###############################################################################
# Block 1. Paths and settings
###############################################################################

paths      <- readRDS("project_paths.rds")
output_dir <- paths$output_dir   # exact (private) tables are written here
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs
tasks_dir  <- paths$tasks_dir    # one file per task in split runs

disclosure_threshold <- 10     # smallest count behind a shown result


###############################################################################
# Block 2. Read the design, the cross-validation results and File 18
###############################################################################

design         <- readRDS(file.path(data_dir, "14_prediction_design.rds"))
settings       <- design$settings
landmark_order <- settings$landmark_order

cv_results  <- readRDS(file.path(data_dir, "17_cv_results.rds"))
predictions <- cv_results$predictions
task_table  <- cv_results$task_table

file18 <- readRDS(file.path(data_dir, "18_performance.rds"))


###############################################################################
# Block 3. Read the final-model tasks (File 19a), with checks
#
# Checks (the script stops if any fails): every landmark's task file is
# present and is the landmark it claims to be; every task used the same R,
# mice and glmnet versions as each other and as the cross-validation.
###############################################################################

final_files <- file.path(tasks_dir, "19a", sprintf("task_%d.rds", seq_along(landmark_order)))
if (any(!file.exists(final_files))) {
  stop("Final-model task files are missing (tasks ",
       paste(which(!file.exists(final_files)), collapse = ", "), "). Run them with File 19a.")
}

final_tasks <- list()
for (t in seq_along(landmark_order)) {
  final_tasks[[t]] <- readRDS(final_files[t])
  if (final_tasks[[t]]$task != t || final_tasks[[t]]$landmark != landmark_order[t]) {
    stop("Final-model task file ", t, " isn't the task it should be.")
  }
  versions_t <- paste(final_tasks[[t]]$r_version, final_tasks[[t]]$mice_version,
                      final_tasks[[t]]$glmnet_version)
  versions_cv <- paste(cv_results$task_log$r_version[1], cv_results$task_log$mice_version[1],
                       cv_results$task_log$glmnet_version[1])
  if (versions_t != versions_cv) {
    stop("Final-model task ", t, " used R / mice / glmnet ", versions_t,
         ", but the cross-validation used ", versions_cv, ".")
  }
}


###############################################################################
# Block 4. The performance measures (the same definitions as File 18)
#
# Copied from File 18, Block 5; the check below confirms that they give
# File 18's fold-by-fold results exactly.
###############################################################################

performance_measures <- function(y, p, w) {
  # Prevalence: the (weighted) share with the outcome; the PR-AUC of a
  # model with no skill.
  prevalence <- sum(w * y) / sum(w)

  # Brier score: the mean squared difference between prediction and outcome.
  brier <- sum(w * (p - y)^2) / sum(w)

  # Log loss: the mean negative log-likelihood. Probabilities are kept a
  # hair away from 0 and 1 so that the logarithm is finite.
  p_bounded <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  log_loss  <- -sum(w * (y * log(p_bounded) + (1 - y) * log(1 - p_bounded))) / sum(w)

  # PR-AUC and ROC-AUC. People are sorted from the highest prediction to the
  # lowest. Each distinct prediction is one cut-off; people with the same
  # prediction are passed together.
  sorted     <- order(p, decreasing = TRUE)
  p_sorted   <- p[sorted]
  events     <- cumsum(w[sorted] * y[sorted])          # events at or above each cut-off
  non_events <- cumsum(w[sorted] * (1 - y[sorted]))    # non-events at or above each cut-off
  last_at_cut_off <- c(p_sorted[-1] != p_sorted[-length(p_sorted)], TRUE)
  events     <- events[last_at_cut_off]
  non_events <- non_events[last_at_cut_off]
  all_events     <- events[length(events)]
  all_non_events <- non_events[length(non_events)]

  # PR-AUC (average precision): the precision at each cut-off, weighted by
  # the share of all events first reached at that cut-off.
  precision <- events / (events + non_events)
  recall    <- events / all_events
  pr_auc    <- sum(diff(c(0, recall)) * precision)

  # ROC-AUC: the chance that an event has a higher prediction than a
  # non-event, with ties counting one half.
  new_events     <- diff(c(0, events))
  new_non_events <- diff(c(0, non_events))
  roc_auc <- sum(new_events * (all_non_events - non_events + new_non_events / 2)) /
    (all_events * all_non_events)

  # Calibration: logistic regressions of the outcome on the log-odds of the
  # prediction. Intercept: with the slope fixed at 1. Slope: the coefficient
  # of the log-odds. A fit that doesn't converge gives NA. The quasibinomial
  # family gives the same estimates as the binomial but accepts weights that
  # aren't whole numbers without a warning.
  log_odds <- qlogis(p_bounded)
  intercept_fit <- suppressWarnings(glm.fit(x = matrix(1, length(y), 1), y = y, weights = w,
                                            offset = log_odds, family = quasibinomial()))
  slope_fit <- suppressWarnings(glm.fit(x = cbind(1, log_odds), y = y, weights = w,
                                        family = quasibinomial()))
  calibration_intercept <- if (intercept_fit$converged) unname(intercept_fit$coefficients[1]) else NA
  calibration_slope     <- if (slope_fit$converged) unname(slope_fit$coefficients[2]) else NA

  c(prevalence = prevalence, pr_auc = pr_auc, roc_auc = roc_auc, brier = brier,
    log_loss = log_loss, calibration_intercept = calibration_intercept,
    calibration_slope = calibration_slope)
}

# Check: the measures reproduce File 18's fold-by-fold results (primary
# analysis, repeat 1, elastic net, weighted).
fold_check <- file18$fold_comparison
for (k in seq_len(nrow(fold_check))) {
  rows <- predictions$strategy == "primary" & predictions$repeat_number == 1 &
    predictions$outer_fold == fold_check$outer_fold[k] & predictions$landmark == fold_check$landmark[k]
  value <- performance_measures(predictions$outcome[rows], predictions$p_elastic_net[rows],
                                predictions$w9_final_weight[rows])[fold_check$measure[k]]
  if (!isTRUE(all.equal(unname(value), fold_check$primary[k], tolerance = 1e-8))) {
    stop("The measures don't reproduce File 18 (", fold_check$landmark[k], ", fold ",
         fold_check$outer_fold[k], ", ", fold_check$measure[k], "). Copy File 18's Block 5 here.")
  }
}


###############################################################################
# Block 5. Importance in each outer split
#
# For each primary task (repeat x fold x landmark): the measures of the
# test fold's predictions, and of the predictions with each feature or
# block shuffled. Importance is the fall in PR-AUC and ROC-AUC and the rise
# in Brier score and log loss (positive = performance got worse without the
# feature's information). Weighted by the Wave-9 final weight.
###############################################################################

primary_tasks <- task_table$task[task_table$strategy == "primary"]
split_importance <- data.frame()

for (t in primary_tasks) {
  rows <- predictions$strategy == "primary" &
    predictions$repeat_number == task_table$repeat_number[t] &
    predictions$outer_fold == task_table$outer_fold[t] &
    predictions$landmark == task_table$landmark[t]
  task_predictions <- predictions[rows, ]
  permuted <- cv_results$permuted[[t]]

  # Checks: the permuted predictions belong to the same people, in the
  # same order; the fold has at least 10 events and 10 non-events.
  if (!identical(permuted$NSID, task_predictions$NSID)) {
    stop("Task ", t, ": the permuted predictions don't match the task's test rows.")
  }
  if (sum(task_predictions$outcome) < disclosure_threshold ||
      sum(1 - task_predictions$outcome) < disclosure_threshold) {
    stop("Task ", t, ": the test fold has fewer than 10 events or non-events.")
  }

  original <- performance_measures(task_predictions$outcome, task_predictions$p_elastic_net,
                                   task_predictions$w9_final_weight)
  for (k in seq_len(nrow(permuted$units))) {
    shuffled <- performance_measures(task_predictions$outcome,
                                     permuted$predictions[, permuted$units$unit[k]],
                                     task_predictions$w9_final_weight)
    split_importance <- rbind(split_importance, data.frame(
      task = t, repeat_number = task_table$repeat_number[t], outer_fold = task_table$outer_fold[t],
      landmark = task_table$landmark[t], unit = permuted$units$unit[k], type = permuted$units$type[k],
      pr_auc_drop       = unname(original["pr_auc"] - shuffled["pr_auc"]),
      roc_auc_drop      = unname(original["roc_auc"] - shuffled["roc_auc"]),
      brier_increase    = unname(shuffled["brier"] - original["brier"]),
      log_loss_increase = unname(shuffled["log_loss"] - original["log_loss"])))
  }
}


###############################################################################
# Block 6. Importance summarized over the splits
#
# For each landmark x feature (or block): the mean over the 25 splits of
# each importance measure, the standard deviation of the PR-AUC fall
# between splits, and the number of splits in which PR-AUC fell. Sorted
# within each landmark and type by the mean PR-AUC fall.
###############################################################################

units <- unique(split_importance[, c("landmark", "unit", "type")])
importance <- data.frame()
for (k in seq_len(nrow(units))) {
  rows <- split_importance[split_importance$landmark == units$landmark[k] &
                             split_importance$unit == units$unit[k], ]
  importance <- rbind(importance, data.frame(
    landmark = units$landmark[k], type = units$type[k],
    unit = sub("^(Feature|Block): ", "", units$unit[k]),
    splits = nrow(rows),
    mean_pr_auc_drop = mean(rows$pr_auc_drop),
    sd_pr_auc_drop = sd(rows$pr_auc_drop),
    splits_pr_auc_fell = sum(rows$pr_auc_drop > 0),
    mean_roc_auc_drop = mean(rows$roc_auc_drop),
    mean_brier_increase = mean(rows$brier_increase),
    mean_log_loss_increase = mean(rows$log_loss_increase)))
}
importance <- importance[order(match(importance$landmark, landmark_order), importance$type,
                               -importance$mean_pr_auc_drop), ]
rownames(importance) <- NULL


###############################################################################
# Block 7. Final-model coefficients
#
# For each landmark and model term: the coefficient averaged over the
# imputed datasets (a term that is absent from a dataset's model counts as
# 0 there; imputations_with_term says how often it was present), the
# standardized coefficient (coefficient x the column's SD, averaged), the
# odds ratio, and in how many datasets the elastic net kept the term (non-
# zero). Coefficients are on the log-odds scale for one unit of the
# feature, or for the category against the reference category.
###############################################################################

all_coefficients <- data.frame()
for (t in seq_along(landmark_order)) {
  all_coefficients <- rbind(all_coefficients,
                            data.frame(landmark = landmark_order[t], final_tasks[[t]]$coefficients))
}

terms <- unique(all_coefficients[, c("landmark", "term", "feature")])
final_coefficients <- data.frame()
for (k in seq_len(nrow(terms))) {
  landmark_rows <- all_coefficients$landmark == terms$landmark[k]
  rows <- all_coefficients[landmark_rows & all_coefficients$term == terms$term[k], ]
  n_imputations <- length(unique(all_coefficients$imputation[landmark_rows]))
  final_coefficients <- rbind(final_coefficients, data.frame(
    landmark = terms$landmark[k], feature = terms$feature[k], term = terms$term[k],
    coefficient = sum(rows$coefficient) / n_imputations,
    standardized_coefficient = sum(rows$coefficient * rows$column_sd) / n_imputations,
    imputations_with_term = nrow(rows),
    imputations_total = n_imputations,
    imputations_kept = sum(rows$coefficient != 0),
    binary = any(rows$binary),
    # For the disclosure check: the smallest counts over the imputed datasets.
    fewest_people_with_term = min(rows$people_with_term),
    fewest_events = min(rows$events),
    fewest_non_events = min(rows$non_events)))
}
final_coefficients$odds_ratio <- exp(final_coefficients$coefficient)
final_coefficients$standardized_coefficient[final_coefficients$term == "(Intercept)"] <- NA

# Disclosure check for models: at least 10 events, 10 non-events and 10
# residual degrees of freedom (people minus model columns).
model_check <- data.frame()
for (t in seq_along(landmark_order)) {
  people <- final_tasks[[t]]$events + final_tasks[[t]]$non_events
  model_check <- rbind(model_check, data.frame(
    landmark = landmark_order[t],
    model_shown = final_tasks[[t]]$events >= disclosure_threshold &
      final_tasks[[t]]$non_events >= disclosure_threshold &
      people - max(final_tasks[[t]]$tuning$columns) - 1 >= disclosure_threshold))
}

# Disclosure check for 0/1 terms (see the header).
final_coefficients$shown <- model_check$model_shown[match(final_coefficients$landmark, model_check$landmark)]
binary_too_small <- final_coefficients$binary &
  (final_coefficients$imputations_with_term < final_coefficients$imputations_total |
     final_coefficients$fewest_people_with_term < disclosure_threshold |
     final_coefficients$fewest_events < disclosure_threshold |
     final_coefficients$fewest_non_events < disclosure_threshold)
final_coefficients$shown[binary_too_small] <- FALSE


###############################################################################
# Block 8. Final-model tuning and diagnostics, by landmark
###############################################################################

final_tuning <- data.frame()
for (t in seq_along(landmark_order)) {
  task_t <- final_tasks[[t]]
  final_tuning <- rbind(final_tuning, data.frame(
    landmark = landmark_order[t],
    imputations = nrow(task_t$tuning),
    chosen_alpha = paste(task_t$tuning$chosen_alpha, collapse = " | "),
    median_lambda = median(task_t$tuning$lambda),
    lambda_at_an_edge = sum(task_t$tuning$lambda_at_largest | task_t$tuning$lambda_at_smallest),
    model_columns = max(task_t$tuning$columns),
    median_nonzero_coefficients = median(task_t$tuning$nonzero_coefficients),
    auxiliaries_kept = paste(task_t$screen$auxiliary[task_t$screen$kept], collapse = " | "),
    tree_fallback_models = length(unique(task_t$switches$model_for)),
    convergence_rows_flagged = sum(task_t$convergence$flag),
    convergence_rows = nrow(task_t$convergence),
    logged_events = task_t$logged_events,
    minutes = round(task_t$total_minutes, 1)))
}


###############################################################################
# Block 9. Importance plots
#
# For each landmark: the mean fall in PR-AUC when each block, and each
# feature, is shuffled, with bars of one standard deviation between the 25
# splits (spread between splits, not a confidence interval).
###############################################################################

type_labels <- c(block = "Blocks of features", feature = "Single features")

pdf(file.path(review_dir, "19_importance_plots.pdf"), width = 8.5, height = 11)
for (landmark in landmark_order) {
  par(mfrow = c(2, 1), mar = c(4.5, 17, 3, 1))
  for (unit_type in c("block", "feature")) {
    rows <- importance[importance$landmark == landmark & importance$type == unit_type, ]
    rows <- rows[order(rows$mean_pr_auc_drop), ]
    lower <- rows$mean_pr_auc_drop - rows$sd_pr_auc_drop
    upper <- rows$mean_pr_auc_drop + rows$sd_pr_auc_drop
    dotchart(rows$mean_pr_auc_drop, labels = rows$unit, pch = 19, xlim = range(c(lower, upper, 0)),
             xlab = "Mean fall in PR-AUC when shuffled (bars: 1 SD between splits)",
             main = paste(landmark, "-", type_labels[unit_type]))
    segments(lower, seq_len(nrow(rows)), upper, seq_len(nrow(rows)))
    abline(v = 0, lty = 2, col = "grey50")
  }
}
invisible(dev.off())


###############################################################################
# Block 10. Save and write the review tables
#
# Disclosure check (before writing):
#   - importance: summaries over outer test folds that each have at least
#     10 events and 10 non-events (checked in Block 5); no counts;
#   - coefficients: models and 0/1 terms checked in Block 7; the counts
#     used for the check stay in the private table;
#   - tuning: model settings and diagnostic counts (no people counted).
###############################################################################

saveRDS(list(split_importance = split_importance, importance = importance,
             final_coefficients = final_coefficients, model_check = model_check,
             final_tuning = final_tuning),
        file.path(data_dir, "19_final_models.rds"))

write.csv(final_coefficients[, c("landmark", "feature", "term", "binary", "fewest_people_with_term",
                                 "fewest_events", "fewest_non_events", "shown")],
          file.path(output_dir, "private_19_term_counts.csv"), row.names = FALSE)

importance_review <- importance
rounded <- c("mean_pr_auc_drop", "sd_pr_auc_drop", "mean_roc_auc_drop", "mean_brier_increase",
             "mean_log_loss_increase")
importance_review[, rounded] <- round(importance_review[, rounded], 5)
importance_review$disclosure_checked <- TRUE
write.csv(importance_review, file.path(review_dir, "19_importance.csv"), row.names = FALSE)

coefficients_review <- final_coefficients[, c("landmark", "feature", "term", "coefficient",
                                              "odds_ratio", "standardized_coefficient",
                                              "imputations_with_term", "imputations_kept", "shown")]
coefficients_review[!coefficients_review$shown,
                    c("coefficient", "odds_ratio", "standardized_coefficient")] <- NA
coefficients_review$note <- ""
coefficients_review$note[!coefficients_review$shown] <-
  "Suppressed: fewer than 10 people, events or non-events with this term"
coefficients_review[, c("coefficient", "odds_ratio", "standardized_coefficient")] <-
  round(coefficients_review[, c("coefficient", "odds_ratio", "standardized_coefficient")], 4)
coefficients_review$disclosure_checked <- TRUE
write.csv(coefficients_review, file.path(review_dir, "19_final_model_coefficients.csv"), row.names = FALSE)

tuning_review <- final_tuning
tuning_review$median_lambda <- signif(tuning_review$median_lambda, 3)
tuning_review$disclosure_checked <- TRUE
write.csv(tuning_review, file.path(review_dir, "19_final_model_tuning.csv"), row.names = FALSE)

disclosure_log <- data.frame(
  table = c("19_importance.csv", "19_importance_plots.pdf", "19_final_model_coefficients.csv",
            "19_final_model_tuning.csv"),
  rule = c("Summaries over outer test folds with at least 10 events and 10 non-events",
           "Drawn from 19_importance.csv",
           "Models with at least 10 events, non-events and residual df; 0/1 terms with at least 10 people, events and non-events in every imputed dataset",
           "Model settings and diagnostic counts only"),
  rows_suppressed = c(0, 0, sum(!final_coefficients$shown), 0)
)
disclosure_log$disclosure_checked <- TRUE
write.csv(disclosure_log, file.path(review_dir, "19_disclosure_log.csv"), row.names = FALSE)


###############################################################################
# Block 11. Summary
#
# Counts and checks only. The results are in the review tables.
###############################################################################

cat("\nFile 19b finished.\n")
cat("  Final-model tasks read and checked:      ", length(final_tasks), "\n")
cat("  Measures reproduce File 18:               yes\n")
cat("  Splits scored for importance:            ", length(primary_tasks), "\n")
cat("  Features and blocks summarized:          ", nrow(importance), "\n")
cat("  Final models shown:                      ", sum(model_check$model_shown), "of",
    nrow(model_check), "\n")
cat("  Coefficients suppressed (0/1 terms):     ", sum(!final_coefficients$shown), "of",
    nrow(final_coefficients), "\n")
cat("  Final-model convergence rows flagged:    ", sum(final_tuning$convergence_rows_flagged), "of",
    sum(final_tuning$convergence_rows), "\n")
cat("  Review: 19_importance.csv, 19_importance_plots.pdf, 19_final_model_coefficients.csv,\n")
cat("    19_final_model_tuning.csv, 19_disclosure_log.csv\n")
