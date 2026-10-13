###############################################################################
# 21. Manuscript tables, figures and final disclosure check
#
# Purpose
#   1. Assemble the manuscript's main tables and figures from the
#      disclosure-checked review tables of Files 4-20.
#   2. Write an index of every manuscript item and its source file
#      (supplementary tables are the review tables themselves).
#   3. Write the reporting-checklist map: where each REFORMS item (Kapoor
#      et al. 2024) and each TRIPOD+AI item (Collins et al. 2024) is
#      reported. Item descriptions are short paraphrases; see the
#      published checklists for their wording.
#   4. Run one final disclosure check over everything in the review folder,
#      which is everything that may be shared (UK Data Service output rules;
#      threshold 10).
#
# This script reads only the review folder: aggregate tables that have
# already passed each script's disclosure check. It reads no participant
# data.
#
# Main tables
#   Table 1  samples: original cohort, target population, Wave-9
#            respondents, prediction cohort; SRET prevalence
#   Table 2  science risk differences and the mathematics odds ratio
#            across the sequential models M0-M4 (imputation, primary)
#   Table 3  the baseline model (M0) across analyses: complete cases,
#            primary imputation, and the sensitivity analyses
#   Table 4  tipping points for data missing not at random
#   Table 5  predictive performance by landmark
#   Table 6  predictive performance by sex, ethnic group and family NS-SEC
#   Table 7  association versus prediction, stage by stage
#
# Figures
#   Figure 1  participant flow
#   Figure 2  missing-data handling in the two analyses (flow chart)
#   Figure 3  discrimination by landmark
#   Figure 4  association versus prediction (File 20's figure)
#   Figure 5  missing-not-at-random sensitivity (File 13's figure)
#
# Final disclosure check (every file in the review folder)
#   - no file named private_*, and only tables (.csv), figures (.pdf) and
#     text (.txt);
#   - every table carries disclosure_checked = TRUE, except the
#     disclosure logs themselves and the predictor matrix (variable
#     structure only);
#   - no person identifier;
#   - no count of people between 1 and 9 in any column that counts people
#     (listed in Block 1). Zeros are listed for a manual look: they are
#     allowed only where they are structural (for example, no refusals at
#     the first sweep).
#   The results are written to 21_final_disclosure_check.csv, and the
#   script stops if any file fails.
#
# Reads
#   review_dir/ (the review tables and figures of Files 0-20)
#
# Writes
#   review_dir/21_table1_samples.csv ... 21_table7_bridge.csv,
#     21_figure1_participant_flow.pdf ... 21_figure5_mnar.pdf,
#     21_output_index.csv, 21_reporting_checklist.csv,
#     21_final_disclosure_check.csv
#
# References
#   Collins GS, et al. (2024). TRIPOD+AI statement. BMJ 385: e078378.
#   Kapoor S, et al. (2024). REFORMS. Science Advances 10: eadk3452.
###############################################################################


###############################################################################
# Block 1. Paths and settings
###############################################################################

paths      <- readRDS("project_paths.rds")
review_dir <- paths$review_dir   # disclosure-checked outputs

disclosure_threshold <- 10

primary_analysis     <- "Primary: target population, design weight"
respondents_analysis <- "Sensitivity: Wave-9 respondents, Wave-9 weight"
observed_analysis    <- "Sensitivity: observed-outcome sample, Wave-9 weight"
supplementary_analysis <- "Supplementary imputation: target population, design weight"

# Columns, in any review table, that count people (for the final check).
people_count_columns <- c("n", "total", "events", "non_events", "n_observed", "n_missing",
                          "n_imputed", "cohort_n", "complete_case_n", "complete_case_events",
                          "complete_case_non_events", "people", "y_events_used",
                          "y_non_events_used", "available_people", "available_events")

# Tables with no data values, which need no disclosure_checked column.
structure_only_files <- c("07_predictor_matrix.csv")


###############################################################################
# Block 2. Read the review tables used here
#
# The script stops if any is missing (run the script that writes it).
###############################################################################

needed <- c("04_sample_flow.csv", "04_landmarks.csv", "12_prevalence.csv",
            "12_risk_differences.csv", "12_coefficients.csv", "12_comparison.csv",
            "12_imputation_number_check.csv", "13_tipping_points.csv", "13_mnar_plots.pdf",
            "14_landmark_design.csv", "18_performance.csv", "18_comparisons.csv",
            "20_stage_changes.csv", "20_bridge_plots.pdf")
missing_files <- needed[!file.exists(file.path(review_dir, needed))]
if (length(missing_files) > 0) {
  stop("These review files are missing: ", paste(missing_files, collapse = ", "))
}

sample_flow       <- read.csv(file.path(review_dir, "04_sample_flow.csv"))
landmarks_04      <- read.csv(file.path(review_dir, "04_landmarks.csv"))
prevalence_12     <- read.csv(file.path(review_dir, "12_prevalence.csv"))
risk_diff_12      <- read.csv(file.path(review_dir, "12_risk_differences.csv"))
coefficients_12   <- read.csv(file.path(review_dir, "12_coefficients.csv"))
comparison_12     <- read.csv(file.path(review_dir, "12_comparison.csv"))
imputation_check  <- read.csv(file.path(review_dir, "12_imputation_number_check.csv"))
tipping_points    <- read.csv(file.path(review_dir, "13_tipping_points.csv"))
landmark_design   <- read.csv(file.path(review_dir, "14_landmark_design.csv"))
performance_18    <- read.csv(file.path(review_dir, "18_performance.csv"))
comparisons_18    <- read.csv(file.path(review_dir, "18_comparisons.csv"))
stage_changes_20  <- read.csv(file.path(review_dir, "20_stage_changes.csv"))

landmark_order <- landmark_design$landmark
models <- c("M0 baseline", "M1 + course route (15)", "M2 + intentions (16)",
            "M3 + HE pathway (20)", "M4 + degree pathway (25)")
contrast_labels <- c("Science liking: like it a lot vs don't like it at all",
                     "Perceived science competence: very good vs no good at all",
                     "Science is a favourite subject: yes vs no")


###############################################################################
# Block 3. Table 1: samples
#
# The counts come from File 4 (original cohort, prediction cohort) and
# File 12 (target population, Wave-9 respondents, observed-outcome sample).
# Check: the observed-outcome sample of File 12 is the prediction cohort.
# Disclosure check: the differences between successive samples (shown in
# Figure 1) must each be 0 or at least 10, so no small group can be found by
# subtraction.
###############################################################################

original_cohort   <- sample_flow$total[sample_flow$sweep == 1][1]
target_population <- prevalence_12$n[prevalence_12$analysis == primary_analysis &
                                       prevalence_12$outcome == "SRET overall"]
w9_respondents    <- prevalence_12$n[prevalence_12$analysis == respondents_analysis &
                                       prevalence_12$outcome == "SRET overall"]
observed_outcome  <- prevalence_12$n[prevalence_12$analysis == observed_analysis &
                                       prevalence_12$outcome == "SRET overall"]
prediction_cohort <- landmarks_04$cohort_n[landmarks_04$landmark == landmark_order[1]]
prediction_events <- landmarks_04$events[landmarks_04$landmark == landmark_order[1]]

if (!isTRUE(all.equal(observed_outcome, prediction_cohort))) {
  stop("File 12's observed-outcome sample and File 4's prediction cohort differ in size.")
}

differences <- c(ineligible_at_32 = original_cohort - target_population,
                 not_responding_at_32 = target_population - w9_respondents,
                 no_observed_outcome = w9_respondents - prediction_cohort)
if (any(differences < 0) || any(differences > 0 & differences < disclosure_threshold)) {
  stop("A difference between successive samples is negative or between 1 and 9: ",
       paste(names(differences), collapse = ", "), ". It can't be shown.")
}

prevalence_primary <- prevalence_12[prevalence_12$analysis == primary_analysis &
                                      prevalence_12$outcome == "SRET overall", ]
prevalence_cohort <- performance_18[performance_18$analysis == "Primary" &
                                      performance_18$landmark == landmark_order[1] &
                                      performance_18$model == "Elastic net" &
                                      performance_18$group == "All" &
                                      performance_18$measure == "prevalence", ]
weighted_cohort   <- prevalence_cohort[prevalence_cohort$weighting == "Weighted", ]
unweighted_cohort <- prevalence_cohort[prevalence_cohort$weighting == "Unweighted", ]

table1 <- data.frame(
  sample = c("Original cohort, interviewed at age 14 (Wave 1)",
             "Target population at 32 (excluding Wave-9 ineligible)",
             "Wave-9 respondents with a Wave-9 weight",
             "Prediction cohort (observed outcome at 32)"),
  used_for = c("", "Imputation analysis (primary)", "Imputation analysis (sensitivity)",
               "Prediction analysis; imputation analysis (sensitivity)"),
  n = c(original_cohort, target_population, w9_respondents, prediction_cohort),
  events = c(NA, NA, NA, prediction_events),
  sret_percent = c("",
                   sprintf("%.1f (%.1f, %.1f), imputed, design weight",
                           prevalence_primary$weighted_percent,
                           prevalence_primary$weighted_percent_low,
                           prevalence_primary$weighted_percent_high),
                   "",
                   sprintf("%.1f unweighted; %.1f (%.1f, %.1f) with the Wave-9 weight",
                           100 * unweighted_cohort$estimate, 100 * weighted_cohort$estimate,
                           100 * weighted_cohort$ci_lower, 100 * weighted_cohort$ci_upper))
)
table1$disclosure_checked <- TRUE
write.csv(table1, file.path(review_dir, "21_table1_samples.csv"), row.names = FALSE, na = "")


###############################################################################
# Block 4. Table 2: associations across the sequential models
#
# Primary imputation analysis. Risk differences in percentage points; the
# mathematics odds ratio per point of liking. Values File 12 suppressed are
# shown as "suppressed".
###############################################################################

table2 <- data.frame(measure = c(paste("Risk difference (points):", contrast_labels),
                                 "Odds ratio: liking mathematics (per point)"))
for (model in models) {
  column <- c()
  for (contrast in contrast_labels) {
    row <- risk_diff_12[risk_diff_12$analysis == primary_analysis & risk_diff_12$model == model &
                          risk_diff_12$contrast == contrast, ]
    value <- "not available"
    if (nrow(row) == 1) value <- sprintf("%.1f (%.1f, %.1f)", 100 * row$risk_difference,
                                         100 * row$low, 100 * row$high)
    column <- c(column, value)
  }
  row <- coefficients_12[coefficients_12$analysis == primary_analysis &
                           coefficients_12$model == model & coefficients_12$term == "math_liking14", ]
  value <- "not available"
  if (nrow(row) == 1) value <- sprintf("%.2f (%.2f, %.2f)", row$odds_ratio, row$odds_ratio_low,
                                       row$odds_ratio_high)
  column <- c(column, value)
  column[grepl("NA", column)] <- "suppressed"
  table2[[model]] <- column
}
table2$disclosure_checked <- TRUE
write.csv(table2, file.path(review_dir, "21_table2_associations.csv"), row.names = FALSE)


###############################################################################
# Block 5. Table 3: the baseline model across analyses
#
# M0 in the complete cases and in each imputation analysis (File 12's
# comparison table).
###############################################################################

table3_analyses <- c("Complete cases: M0's own sample", primary_analysis, respondents_analysis,
                     observed_analysis, supplementary_analysis)
table3_quantities <- c(paste("Risk difference:", contrast_labels), "Odds ratio: math_liking14")

table3 <- data.frame(measure = c(paste("Risk difference (points):", contrast_labels),
                                 "Odds ratio: liking mathematics (per point)"))
for (analysis in table3_analyses) {
  column <- c()
  for (quantity in table3_quantities) {
    row <- comparison_12[comparison_12$model == "M0 baseline" & comparison_12$analysis == analysis &
                           comparison_12$quantity == quantity, ]
    value <- "not available"
    if (nrow(row) == 1 && startsWith(quantity, "Risk difference")) {
      value <- sprintf("%.1f (%.1f, %.1f)", 100 * row$estimate, 100 * row$low, 100 * row$high)
    }
    if (nrow(row) == 1 && startsWith(quantity, "Odds ratio")) {
      value <- sprintf("%.2f (%.2f, %.2f)", row$estimate, row$low, row$high)
    }
    column <- c(column, value)
  }
  column[grepl("NA", column)] <- "suppressed"
  table3[[analysis]] <- column
}
table3$disclosure_checked <- TRUE
write.csv(table3, file.path(review_dir, "21_table3_analyses.csv"), row.names = FALSE)


###############################################################################
# Block 6. Table 4: tipping points for data missing not at random (File 13)
###############################################################################

table4 <- data.frame(
  shift_applied_to = tipping_points$family,
  estimand = tipping_points$estimand,
  conclusion_under_missing_at_random = tipping_points$conclusion_under_mar,
  tipping_delta_below_0 = tipping_points$tipping_delta_negative,
  tipping_odds_ratio_below_1 = tipping_points$tipping_odds_ratio_negative,
  tipping_delta_above_0 = tipping_points$tipping_delta_positive,
  tipping_odds_ratio_above_1 = tipping_points$tipping_odds_ratio_positive
)
table4$disclosure_checked <- TRUE
write.csv(table4, file.path(review_dir, "21_table4_tipping_points.csv"), row.names = FALSE, na = "")


###############################################################################
# Block 7. Table 5: predictive performance by landmark
#
# Primary analysis, weighted by the Wave-9 final weight. Elastic net and
# logistic regression; gains are each landmark minus the one before
# (elastic net, paired 95% intervals).
###############################################################################

table5 <- data.frame()
for (landmark in landmark_order) {
  rows <- performance_18[performance_18$analysis == "Primary" & performance_18$weighting == "Weighted" &
                           performance_18$group == "All" & performance_18$landmark == landmark, ]
  en <- rows[rows$model == "Elastic net", ]
  lr <- rows[rows$model == "Logistic regression (unpenalized)", ]
  gains <- comparisons_18[comparisons_18$comparison == "Landmark minus the landmark before" &
                            comparisons_18$analysis == "Primary" & comparisons_18$weighting == "Weighted" &
                            comparisons_18$first == paste(landmark, "Elastic net"), ]
  pr_gain  <- gains[gains$measure == "pr_auc", ]
  roc_gain <- gains[gains$measure == "roc_auc", ]
  # The first landmark has no landmark before it, so no gain.
  pr_gain_text  <- ""
  roc_gain_text <- ""
  if (nrow(pr_gain) == 1) {
    pr_gain_text <- sprintf("%+.3f (%+.3f, %+.3f)", pr_gain$difference, pr_gain$ci_lower,
                            pr_gain$ci_upper)
  }
  if (nrow(roc_gain) == 1) {
    roc_gain_text <- sprintf("%+.3f (%+.3f, %+.3f)", roc_gain$difference, roc_gain$ci_lower,
                             roc_gain$ci_upper)
  }
  table5 <- rbind(table5, data.frame(
    landmark = landmark,
    prevalence_percent = sprintf("%.1f", 100 * en$estimate[en$measure == "prevalence"]),
    elastic_net_pr_auc = sprintf("%.3f (%.3f, %.3f)", en$estimate[en$measure == "pr_auc"],
                                 en$ci_lower[en$measure == "pr_auc"], en$ci_upper[en$measure == "pr_auc"]),
    elastic_net_roc_auc = sprintf("%.3f (%.3f, %.3f)", en$estimate[en$measure == "roc_auc"],
                                  en$ci_lower[en$measure == "roc_auc"], en$ci_upper[en$measure == "roc_auc"]),
    elastic_net_brier = sprintf("%.4f", en$estimate[en$measure == "brier"]),
    elastic_net_calibration_intercept = sprintf("%.2f (%.2f, %.2f)",
                                                en$estimate[en$measure == "calibration_intercept"],
                                                en$ci_lower[en$measure == "calibration_intercept"],
                                                en$ci_upper[en$measure == "calibration_intercept"]),
    elastic_net_calibration_slope = sprintf("%.2f (%.2f, %.2f)",
                                            en$estimate[en$measure == "calibration_slope"],
                                            en$ci_lower[en$measure == "calibration_slope"],
                                            en$ci_upper[en$measure == "calibration_slope"]),
    logistic_pr_auc = sprintf("%.3f (%.3f, %.3f)", lr$estimate[lr$measure == "pr_auc"],
                              lr$ci_lower[lr$measure == "pr_auc"], lr$ci_upper[lr$measure == "pr_auc"]),
    logistic_calibration_slope = sprintf("%.2f (%.2f, %.2f)",
                                         lr$estimate[lr$measure == "calibration_slope"],
                                         lr$ci_lower[lr$measure == "calibration_slope"],
                                         lr$ci_upper[lr$measure == "calibration_slope"]),
    pr_auc_gain = pr_gain_text,
    roc_auc_gain = roc_gain_text))
}
table5$disclosure_checked <- TRUE
write.csv(table5, file.path(review_dir, "21_table5_performance.csv"), row.names = FALSE)


###############################################################################
# Block 8. Table 6: performance by subgroup
#
# Elastic net, primary analysis, weighted, at L14 and L25a. Groups that
# File 18 suppressed are shown as "suppressed".
###############################################################################

subgroup_rows <- performance_18[performance_18$analysis == "Primary" & performance_18$model == "Elastic net" &
                                  performance_18$weighting == "Weighted" &
                                  performance_18$group_variable %in% c("All", "Sex", "Ethnic group",
                                                                       "Family NS-SEC at 14"), ]
groups <- unique(subgroup_rows[, c("group_variable", "group")])
table6 <- data.frame()
for (k in seq_len(nrow(groups))) {
  new_row <- data.frame(group_variable = groups$group_variable[k], group = groups$group[k])
  for (landmark in c(landmark_order[1], "L25a")) {
    rows <- subgroup_rows[subgroup_rows$group_variable == groups$group_variable[k] &
                            subgroup_rows$group == groups$group[k] & subgroup_rows$landmark == landmark, ]
    prevalence <- rows[rows$measure == "prevalence", ]
    pr_auc     <- rows[rows$measure == "pr_auc", ]
    roc_auc    <- rows[rows$measure == "roc_auc", ]
    intercept  <- rows[rows$measure == "calibration_intercept", ]
    new_row[[paste(landmark, "prevalence_percent")]] <- sprintf("%.1f", 100 * prevalence$estimate)
    new_row[[paste(landmark, "pr_auc")]] <- sprintf("%.3f (%.3f, %.3f)", pr_auc$estimate,
                                                    pr_auc$ci_lower, pr_auc$ci_upper)
    new_row[[paste(landmark, "roc_auc")]] <- sprintf("%.3f (%.3f, %.3f)", roc_auc$estimate,
                                                     roc_auc$ci_lower, roc_auc$ci_upper)
    new_row[[paste(landmark, "calibration_intercept")]] <- sprintf("%.2f (%.2f, %.2f)", intercept$estimate,
                                                                   intercept$ci_lower, intercept$ci_upper)
  }
  table6 <- rbind(table6, new_row)
}
for (column in names(table6)[-(1:2)]) {
  table6[[column]][grepl("NA", table6[[column]])] <- "suppressed"
}
table6$disclosure_checked <- TRUE
write.csv(table6, file.path(review_dir, "21_table6_subgroups.csv"), row.names = FALSE)


###############################################################################
# Block 9. Table 7: association versus prediction (File 20, Result 2)
#
# The stages with a prediction landmark, and the change from M0 to M4.
# Risk-difference changes in percentage points.
###############################################################################

bridge_rows <- stage_changes_20[!(stage_changes_20$change %in% c("Course route at 15 (M0 to M1)",
                                                               "Intentions at 16 (M1 to M2)")), ]
table7 <- data.frame(
  stage = bridge_rows$change,
  contrast = bridge_rows$contrast,
  change_in_risk_difference_points = sprintf("%+.2f (%+.2f, %+.2f)",
                                             100 * bridge_rows$change_in_risk_difference,
                                             100 * bridge_rows$change_low, 100 * bridge_rows$change_high),
  association = bridge_rows$association_result,
  landmark = bridge_rows$landmark,
  pr_auc_gain = sprintf("%+.3f (%+.3f, %+.3f)", bridge_rows$pr_auc_gain, bridge_rows$pr_auc_gain_low,
                        bridge_rows$pr_auc_gain_high),
  added_block_importance = sprintf("%.3f (%.3f, %.3f)", bridge_rows$added_block_importance,
                                   bridge_rows$added_block_importance_low,
                                   bridge_rows$added_block_importance_high),
  prediction = bridge_rows$prediction_result,
  agreement = bridge_rows$agreement
)
for (column in c("change_in_risk_difference_points", "pr_auc_gain", "added_block_importance")) {
  table7[[column]][grepl("NA", table7[[column]])] <- ""
}
table7$disclosure_checked <- TRUE
write.csv(table7, file.path(review_dir, "21_table7_bridge.csv"), row.names = FALSE, na = "")


###############################################################################
# Block 10. Figure 1: participant flow
#
# Boxes for the four samples of Table 1, with the people who leave between
# them (checked in Block 3).
###############################################################################

pdf(file.path(review_dir, "21_figure1_participant_flow.pdf"), width = 8.5, height = 9)
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 10))

main_x <- 2.8   # centre of the main boxes (half-width 2.6)
side_x <- 7.9   # centre of the side boxes (half-width 1.95)
box_y  <- c(9, 6.6, 4.2, 1.8)
main_text <- c(
  paste0("Original cohort interviewed at age 14 (Wave 1)\nn = ", format(original_cohort, big.mark = ",")),
  paste0("Target population at age 32\nn = ", format(target_population, big.mark = ","),
         "\nImputation analysis (primary)"),
  paste0("Wave-9 respondents with a Wave-9 weight\nn = ", format(w9_respondents, big.mark = ",")),
  paste0("Prediction cohort: observed outcome at 32\nn = ", format(prediction_cohort, big.mark = ","),
         ", SRET = ", format(prediction_events, big.mark = ","),
         "\nPrediction analysis")
)
side_text <- c(
  paste0("Ineligible at Wave 9\n(died, in prison or on probation)\nn = ",
         format(differences["ineligible_at_32"], big.mark = ",")),
  paste0("Not interviewed at Wave 9, or no Wave-9 weight\n(outcome imputed in the imputation analysis)\nn = ",
         format(differences["not_responding_at_32"], big.mark = ",")),
  paste0("No observed outcome at 32\nn = ", format(differences["no_observed_outcome"], big.mark = ","))
)
for (k in seq_along(box_y)) {
  rect(main_x - 2.6, box_y[k] - 0.8, main_x + 2.6, box_y[k] + 0.8)
  text(main_x, box_y[k], main_text[k], cex = 0.85)
}
for (k in seq_along(side_text)) {
  arrows(main_x, box_y[k] - 0.8, main_x, box_y[k + 1] + 0.8, length = 0.1)
  side_y <- (box_y[k] + box_y[k + 1]) / 2
  segments(main_x, side_y, side_x - 1.95, side_y)
  rect(side_x - 1.95, side_y - 0.65, side_x + 1.95, side_y + 0.65)
  text(side_x, side_y, side_text[k], cex = 0.75)
}
invisible(dev.off())


###############################################################################
# Block 11. Figure 2: missing-data handling in the two analyses
#
# A flow chart of the steps (REFORMS 4c; CLS user guide). The only numbers
# are settings: the number of imputations (File 12's check table) and the
# fold design.
###############################################################################

m_imputations <- max(imputation_check$imputations[imputation_check$analysis == primary_analysis])

imputation_steps <- c(
  "Target population: analysis variables, outcome at 32,\nCLS predictors of non-response, design weight and stratum",
  "Complete-case analysis first (reference)",
  paste0("Chained-equation imputation of everyone: ", m_imputations,
         " imputations x 30 iterations;\nrouted questions in parts; outcome in parts"),
  "Diagnostics: trace plots, convergence table, observed vs imputed;\nrepresentativeness check (masking at age 14)",
  "Survey-weighted models in each completed dataset;\nRubin's rules (primary and sensitivity samples)",
  "Missing not at random: delta adjustment to a tipping point"
)
prediction_steps <- c(
  "Prediction cohort: observed outcome at 32\n(Wave-9 respondents, Wave-9 weight)",
  "5 repeats x 5 outer folds, whole schools;\n5 inner folds for tuning",
  "In each training set: auxiliary screen; 5 imputations x 30 iterations,\nmodels fitted to training rows; test outcomes masked",
  "Elastic net tuned by inner cross-validation;\nlogistic regression and null model as benchmarks",
  "Predictions averaged over imputations;\nperformance in the test folds, weighted",
  "Bootstrap of schools: intervals, landmark gains,\nsubgroups, calibration, importance"
)

pdf(file.path(review_dir, "21_figure2_missing_data_flow.pdf"), width = 11, height = 8.5)
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 10))
text(2.5, 9.7, "Imputation analysis (population association)", font = 2)
text(7.5, 9.7, "Prediction analysis (predictability)", font = 2)
step_y <- seq(8.7, 1.2, length.out = length(imputation_steps))
for (k in seq_along(step_y)) {
  rect(0.2, step_y[k] - 0.55, 4.8, step_y[k] + 0.55)
  text(2.5, step_y[k], imputation_steps[k], cex = 0.7)
  rect(5.2, step_y[k] - 0.55, 9.8, step_y[k] + 0.55)
  text(7.5, step_y[k], prediction_steps[k], cex = 0.7)
}
for (k in seq_len(length(step_y) - 1)) {
  arrows(2.5, step_y[k] - 0.55, 2.5, step_y[k + 1] + 0.55, length = 0.08)
  arrows(7.5, step_y[k] - 0.55, 7.5, step_y[k + 1] + 0.55, length = 0.08)
}
text(5, 0.4, "Completed datasets from the imputation analysis are never used in the prediction analysis.",
     cex = 0.75, font = 3)
invisible(dev.off())


###############################################################################
# Block 12. Figure 3: discrimination by landmark
#
# PR-AUC (with the prevalence, a model with no skill, as a dashed line) and
# ROC-AUC, for the elastic net and logistic regression, with 95% intervals.
###############################################################################

figure_rows <- performance_18[performance_18$analysis == "Primary" & performance_18$weighting == "Weighted" &
                                performance_18$group == "All", ]
x <- seq_along(landmark_order)

pdf(file.path(review_dir, "21_figure3_performance.pdf"), width = 10, height = 4.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
for (measure in c("pr_auc", "roc_auc")) {
  en <- figure_rows[figure_rows$model == "Elastic net" & figure_rows$measure == measure, ]
  lr <- figure_rows[figure_rows$model == "Logistic regression (unpenalized)" & figure_rows$measure == measure, ]
  en <- en[match(landmark_order, en$landmark), ]
  lr <- lr[match(landmark_order, lr$landmark), ]
  reference <- 0.5
  if (measure == "pr_auc") {
    reference <- figure_rows$estimate[figure_rows$model == "Elastic net" & figure_rows$measure == "prevalence" &
                                        figure_rows$landmark == landmark_order[1]]
  }
  y_range <- range(c(reference, en$ci_lower, en$ci_upper, lr$ci_lower, lr$ci_upper))
  plot(x - 0.08, en$estimate, pch = 19, ylim = y_range, xlim = c(0.5, length(x) + 0.5), xaxt = "n",
       xlab = "Landmark", ylab = toupper(sub("_", "-", measure)),
       main = paste(toupper(sub("_", "-", measure)), "by landmark (weighted)"), cex.main = 0.95)
  axis(1, at = x, labels = landmark_order)
  segments(x - 0.08, en$ci_lower, x - 0.08, en$ci_upper)
  points(x + 0.08, lr$estimate, pch = 1)
  segments(x + 0.08, lr$ci_lower, x + 0.08, lr$ci_upper, lty = 3)
  abline(h = reference, lty = 2, col = "grey50")
  legend("topleft", legend = c("Elastic net", "Logistic regression", "No skill"),
         pch = c(19, 1, NA), lty = c(1, 3, 2), col = c("black", "black", "grey50"), bty = "n", cex = 0.8)
}
invisible(dev.off())


###############################################################################
# Block 13. Figures 4 and 5: copies of File 20's and File 13's figures
###############################################################################

invisible(file.copy(file.path(review_dir, "20_bridge_plots.pdf"),
                    file.path(review_dir, "21_figure4_bridge.pdf"), overwrite = TRUE))
invisible(file.copy(file.path(review_dir, "13_mnar_plots.pdf"),
                    file.path(review_dir, "21_figure5_mnar.pdf"), overwrite = TRUE))


###############################################################################
# Block 14. Index of manuscript items
###############################################################################

output_index <- data.frame(
  item = c(paste("Table", 1:7), paste("Figure", 1:5),
           paste0("Supplementary table S", 1:12), paste0("Supplementary figure S", 1:4)),
  title = c("Samples and SRET prevalence",
            "Science associations across the sequential models",
            "The baseline model across analyses",
            "Tipping points for data missing not at random",
            "Predictive performance by landmark",
            "Predictive performance by subgroup",
            "Association versus prediction, stage by stage",
            "Participant flow",
            "Missing-data handling in the two analyses",
            "Discrimination by landmark",
            "Association versus prediction",
            "Sensitivity to data missing not at random",
            "Missingness by variable",
            "Complete-case and incomplete-case comparison",
            "Auxiliary variables",
            "Imputation diagnostics (convergence; observed vs imputed)",
            "Representativeness check",
            "Number of imputations (CLS rule) and all coefficients",
            "Sample size by landmark",
            "Landmark variables and fold balance",
            "All performance measures and comparisons (incl. outcome omitted, features observed vs imputed)",
            "Permutation importance",
            "Final-model coefficients",
            "The science signal across stages",
            "Trace plots (main and supplementary imputations)",
            "Calibration plots",
            "Importance plots",
            "Cross-validation diagnostics"),
  file = c("21_table1_samples.csv", "21_table2_associations.csv", "21_table3_analyses.csv",
           "21_table4_tipping_points.csv", "21_table5_performance.csv", "21_table6_subgroups.csv",
           "21_table7_bridge.csv",
           "21_figure1_participant_flow.pdf", "21_figure2_missing_data_flow.pdf",
           "21_figure3_performance.pdf", "21_figure4_bridge.pdf", "21_figure5_mnar.pdf",
           "04_missingness.csv", "05_complete_vs_incomplete.csv", "06_auxiliary_screen.csv",
           "10_convergence.csv; 10_observed_vs_imputed.csv", "11_representativeness.csv",
           "12_imputation_number_check.csv; 12_coefficients.csv", "15_sample_size.csv",
           "14_landmark_variables.csv; 14_fold_balance.csv",
           "18_performance.csv; 18_comparisons.csv; 18_strategy_comparison_by_fold.csv; 18_school_differences.csv",
           "19_importance.csv", "19_final_model_coefficients.csv", "20_orientation_across_stages.csv",
           "10_trace_plots_main.pdf; 10_trace_plots_supplementary.pdf", "18_calibration_plots.pdf",
           "19_importance_plots.pdf", "17_convergence_summary.csv; 17_fallback_summary.csv; 17_model_summary.csv")
)
output_index$disclosure_checked <- TRUE
write.csv(output_index, file.path(review_dir, "21_output_index.csv"), row.names = FALSE)


###############################################################################
# Block 15. Reporting-checklist map (REFORMS; TRIPOD+AI)
#
# Short paraphrases of each item and where it is reported: the Methods
# section numbers of the manuscript draft, the tables and figures above, or
# the script. "Not applicable" is stated with the reason.
###############################################################################

reforms <- data.frame(
  guideline = "REFORMS",
  item = c("1a", "1b", "1c", "2a", "2b", "2c", "2d", "2e", "3a", "3b", "3c", "3d", "3e", "3f", "3g",
           "4a", "4b", "4c", "5a", "5b", "5c", "5d", "5e", "5f", "6a", "6b", "6c", "7a", "7b", "7c",
           "8a", "8b"),
  topic = c("Population of the claim", "Why this population", "Why prediction methods",
            "Dataset and identifier", "Code and version", "Computing environment", "README",
            "Reproduction script", "Data sources and timing", "Sampling frame",
            "Why the data suit the task", "Outcome definition and summary",
            "Sample size and outcome frequency", "Missing data by outcome class",
            "Representativeness", "Exclusions", "Impossible or corrupt values",
            "Transformations and flow chart", "Models trained", "Why these models",
            "Evaluation design", "Model selection", "Tuning", "Baselines",
            "Preprocessing on training data only", "Dependencies between splits",
            "Legitimate features", "Metrics", "Uncertainty", "Statistical tests",
            "External validity", "Where findings may not hold"),
  where_reported = c("Methods 6.1", "Methods 6.1", "Methods 6.1", "Methods 1 (UKDS SN 5545)",
                     "Methods 9", "Methods 9; File 0 package versions", "README.md",
                     "Reproduction script (runs Files 0-21 in order)", "Methods 1",
                     "Methods 1", "Methods 1, 6.1", "Methods 2; Table 1", "Table 1; Figure 1",
                     "04_landmarks.csv (complete cases by outcome); 04_missingness.csv",
                     "Methods 5.6, 6.1; Table 3; Supplementary table S5",
                     "Methods 4; Figure 1", "Methods 5.4 (deposited-data correction); File 2",
                     "Methods 3, 5, 6.5; Figure 2", "Methods 6.6", "Methods 6.1, 6.6",
                     "Methods 6.4", "Methods 6.6", "Methods 6.6",
                     "Methods 6.6; Table 5", "Methods 6.5; File 16 (leakage test)",
                     "Methods 6.4 (whole schools in folds)", "Methods 6.2; File 14 checks",
                     "Methods 6.6, 6.7", "Methods 6.7", "Methods 6.7 (intervals; no tests)",
                     "Methods 13 (internal validation only)", "Methods 13")
)

tripod <- data.frame(
  guideline = "TRIPOD+AI",
  item = c("1", "2", "3a", "3b", "3c", "4", "5a", "5b", "6a", "6b", "6c", "7", "8a", "8b", "8c",
           "9a", "9b", "9c", "10", "11", "12a", "12b", "12c", "12d", "12e", "12f", "12g", "13",
           "14", "15", "16", "17", "18a", "18b", "18c", "18d", "18e", "18f", "19", "20a", "20b",
           "20c", "21", "22", "23a", "23b", "24", "25", "26", "27a", "27b", "27c"),
  topic = c("Title", "Abstract", "Context and existing models", "Target population and purpose",
            "Known inequalities", "Objectives", "Data sources and representativeness",
            "Data collection dates", "Setting and centres", "Eligibility", "Treatments",
            "Preprocessing and quality checks", "Outcome and time horizon", "Outcome assessors",
            "Blinding of outcome", "Choice of predictors", "Predictor definitions and timing",
            "Predictor assessors", "Sample size", "Missing data", "Use of the data",
            "Handling of predictors", "Model type, building, tuning, validation",
            "Heterogeneity across clusters", "Performance measures", "Model updating",
            "How predictions were calculated", "Class imbalance", "Fairness", "Model output",
            "Development versus evaluation data", "Ethics", "Funding", "Conflicts of interest",
            "Protocol", "Registration", "Data sharing", "Code sharing", "Public involvement",
            "Participant flow", "Characteristics and missing data",
            "Evaluation versus development predictors", "Participants and events per analysis",
            "Full model specification", "Performance with intervals, incl. subgroups",
            "Heterogeneity across clusters (results)", "Updated model", "Interpretation",
            "Limitations", "Poor or missing inputs in use", "Users and expertise",
            "Next steps"),
  where_reported = c("Title", "Abstract", "Introduction",
                     "Methods 6.1 (not intended for individual use)", "Introduction; Methods 11",
                     "Introduction; Methods 6.1", "Methods 1, 5.6", "Methods 1",
                     "Methods 1 (schools in England)", "Methods 4", "Not applicable (no treatments)",
                     "Methods 3, 4; Files 2-3", "Methods 2; Methods 11",
                     "Methods 2 (occupation coded to SOC2010 by the survey)", "Methods 11 (not blinded)",
                     "Methods 3, 6.2", "Methods 3, 6.2; 14_landmark_variables.csv",
                     "Not applicable (survey responses)", "Methods 6.3; 15_sample_size.csv",
                     "Methods 5, 6.5; Figure 2", "Methods 6.4", "Methods 6.6",
                     "Methods 6.4-6.6", "Methods 6.8", "Methods 6.7",
                     "Methods 11 (none)", "Methods 6.5 (averaged over imputations)",
                     "Methods 11 (none used)", "Methods 6.8", "Methods 11 (probabilities; no thresholds)",
                     "Methods 6.4 (internal validation, same cohort)", "Methods 10", "Methods 10",
                     "Methods 10", "Methods 10 (protocol)", "Methods 10 (not registered)",
                     "Methods 9 (UKDS End User Licence)", "Methods 9", "Methods 10 (none)",
                     "Figure 1; Table 1", "Table 1; Supplementary table S1",
                     "Not applicable (internal validation)", "Table 1; 14_fold_balance.csv; 15_sample_size.csv",
                     "Supplementary table S11 (averaged over imputations)", "Tables 5-6; Figure 3",
                     "18_school_differences.csv (aggregate only)", "Not applicable (no updating)",
                     "Discussion", "Discussion; Methods 13",
                     "Methods 11 (not intended for use)", "Not applicable (not intended for use)",
                     "Discussion")
)
reporting_checklist <- rbind(reforms, tripod)
reporting_checklist$disclosure_checked <- TRUE
write.csv(reporting_checklist, file.path(review_dir, "21_reporting_checklist.csv"), row.names = FALSE)


###############################################################################
# Block 16. Final disclosure check of the review folder
#
# Every file, including this script's own. The rules are in the header.
###############################################################################

review_files <- sort(list.files(review_dir))
final_check <- data.frame()

for (f in review_files) {
  extension <- tolower(tools::file_ext(f))
  problem <- ""
  zeros <- 0
  if (startsWith(f, "private_")) problem <- "Private file in the review folder"
  if (!(extension %in% c("csv", "pdf", "txt"))) problem <- "Not a table, figure or text file"

  if (extension == "csv" && problem == "") {
    table_f <- read.csv(file.path(review_dir, f), check.names = FALSE)
    needs_flag <- !endsWith(f, "_disclosure_log.csv") && !(f %in% structure_only_files)
    has_flag <- "disclosure_checked" %in% names(table_f) && all(table_f$disclosure_checked %in% TRUE)
    if (needs_flag && !has_flag) problem <- "No disclosure_checked = TRUE on every row"
    if ("NSID" %in% names(table_f)) problem <- "Contains a person identifier"
    count_values <- suppressWarnings(as.numeric(unlist(table_f[intersect(people_count_columns,
                                                                         names(table_f))])))
    if (any(count_values >= 1 & count_values < disclosure_threshold, na.rm = TRUE)) {
      problem <- "A count of people between 1 and 9"
    }
    zeros <- sum(count_values == 0, na.rm = TRUE)
  }

  result <- "Passed"
  if (problem != "") result <- "Failed"
  final_check <- rbind(final_check, data.frame(
    file = f, type = extension, result = result, problem = problem,
    zero_counts_to_look_at = zeros))
}
final_check$disclosure_checked <- TRUE
write.csv(final_check, file.path(review_dir, "21_final_disclosure_check.csv"), row.names = FALSE)

if (any(final_check$result == "Failed")) {
  stop("The final disclosure check failed for: ",
       paste(final_check$file[final_check$result == "Failed"], collapse = ", "),
       ". See 21_final_disclosure_check.csv; remove or fix these files before sharing.")
}


###############################################################################
# Block 17. Summary
###############################################################################

cat("\nFile 21 finished.\n")
cat("  Main tables written:                 7 (21_table1 ... 21_table7)\n")
cat("  Figures written:                     5 (21_figure1 ... 21_figure5)\n")
cat("  Reporting-checklist items mapped:   ", nrow(reporting_checklist), "(REFORMS",
    nrow(reforms), "; TRIPOD+AI", nrow(tripod), ")\n")
cat("  Final disclosure check:              all", nrow(final_check), "files passed\n")
cat("  Files with zero counts to look at:  ", sum(final_check$zero_counts_to_look_at > 0), "\n")
if (any(final_check$zero_counts_to_look_at > 0)) {
  cat("   ", paste(final_check$file[final_check$zero_counts_to_look_at > 0], collapse = ", "), "\n")
}
cat("  Review: 21_output_index.csv lists every manuscript item and its source file.\n")
