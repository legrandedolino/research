###############################################################################
# 18. Performance and landmark comparisons
#
# Purpose
#   Measure how well the prediction models of File 17a distinguish people in
#   SRET employment at 32 from those who are not, at each landmark, and how
#   much each landmark adds to the one before. This is a study of
#   predictability: the measures describe the prediction task in the Wave-9
#   target population (through the Wave-9 final weight). They are not a tool
#   for individuals, and no individual prediction is written.
#
# How the measures are computed (REFORMS 7; TRIPOD+AI 12e, 23)
#   In each repeat of the cross-validation, the five outer test folds
#   together give one out-of-fold prediction for every person in the
#   cohort. Each measure is computed on a repeat's predictions and then
#   averaged over the repeats.
#     - PR-AUC (average precision): the primary measure, because the outcome
#       is rare (Saito and Rehmsmeier 2015). A model with no skill scores
#       the outcome's prevalence, which is reported beside it.
#     - ROC-AUC (the C-statistic).
#     - Brier score and log loss: the overall accuracy of the probabilities.
#     - Calibration intercept (calibration in the large; ideal 0) and
#       calibration slope (ideal 1; below 1 means the predictions are too
#       extreme) (Van Calster et al. 2019).
#   Measures are weighted by the Wave-9 final weight (primary) and also
#   given unweighted. Models: the elastic net (primary), unpenalized
#   logistic regression (benchmark) and the null model (the training
#   prevalence). Calibration is reported descriptively; landmarks are
#   compared on discrimination.
#
# Uncertainty
#   95% confidence intervals come from a bootstrap that keeps schools whole:
#   the primary sampling units (schools at age 14) are resampled with
#   replacement within each stratum, and every measure is recomputed on the
#   out-of-fold predictions of the resampled people (the models are not
#   refitted). Every comparison uses the same resampled people on both sides,
#   so comparisons are paired. Differences between folds are not used as
#   confidence intervals.
#
# Results
#   1. Performance by analysis, landmark, model and weighting.
#   2. Paired comparisons: each landmark against the one before (what the
#      later landmark adds), and the models against each other.
#   3. The primary analysis against the outcome-omitted sensitivity analysis
#      (File 17a), on repeat 1, which both share: overall (bootstrap) and
#      fold by fold. The shift in calibration slope between the two
#      estimates the effect of using the outcome in the training-set
#      imputation (Sisk et al. 2023); the outcome-omitted slope's departure
#      from 1 estimates overfitting. The primary analysis is not switched.
#   4. Calibration split between test rows with every feature observed and
#      rows with at least one imputed feature.
#   5. Binned calibration plots (repeat 1). Each bin has at least 10 events
#      and 10 non-events; no individual points are drawn.
#   6. Performance by sex, ethnic group and family NS-SEC at 14 (elastic
#      net, primary analysis, weighted).
#   7. Differences between schools, in aggregate only: the standard
#      deviation of school-level calibration intercepts, from a random-
#      intercept logistic model. No result is shown for any single school.
#   8. For the comparison with the imputation analysis (File 20): the
#      elastic net's fall in performance when the science and mathematics
#      orientation block, each science item, and each landmark's added
#      block are shuffled (permutation importance; shuffled predictions from
#      File 17a), with paired bootstrap intervals.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   A result for a group of people is shown only if the group has at least
#   10 events and 10 non-events (unweighted). No counts of people, no
#   individual predictions and no minima or maxima are written to the
#   review folder.
#
# Run time
#   About 20-35 minutes (1,000 bootstrap replicates). On a Mac, keep the
#   computer awake while it runs (for example, run "caffeinate -dimsu" in
#   Terminal).
#
# Reads
#   data_dir/14_prediction_design.rds
#   data_dir/17_cv_results.rds
#
# Writes
#   data_dir/18_performance.rds              estimates, bootstrap replicates
#                                            and the tables below
#   output_dir/private_18_group_counts.csv   events and non-events per group
#   review_dir/18_performance.csv, 18_comparisons.csv,
#     18_strategy_comparison_by_fold.csv, 18_calibration_bins.csv,
#     18_calibration_plots.pdf, 18_school_differences.csv,
#     18_disclosure_log.csv
#
# References
#   Collins GS, et al. (2024). TRIPOD+AI statement. BMJ 385: e078378.
#   Kapoor S, et al. (2024). REFORMS. Science Advances 10: eadk3452.
#   Saito T, Rehmsmeier M (2015). PLoS ONE 10(3): e0118432.
#   Sisk R, et al. (2023). Statistical Methods in Medical Research 32(8):
#   1461-1477.
#   Van Calster B, et al. (2019). Calibration: the Achilles heel of
#   predictive analytics. BMC Medicine 17: 230.
###############################################################################

suppressPackageStartupMessages(library(lme4))


###############################################################################
# Block 1. Paths and settings
###############################################################################

paths      <- readRDS("project_paths.rds")
output_dir <- paths$output_dir   # exact (private) tables are written here
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs

disclosure_threshold <- 10     # fewest events (and non-events) behind a shown result
bootstrap_replicates <- 1000   # bootstrap samples of schools
bootstrap_seed       <- 1800   # replicate b uses seed bootstrap_seed + b
calibration_bins     <- 10     # starting number of bins in the calibration plots


###############################################################################
# Block 2. Read the prediction design (File 14) and the predictions (File 17b)
###############################################################################

design         <- readRDS(file.path(data_dir, "14_prediction_design.rds"))
cohort         <- design$cohort
settings       <- design$settings
landmark_order <- settings$landmark_order

cv_results  <- readRDS(file.path(data_dir, "17_cv_results.rds"))
predictions <- cv_results$predictions

n <- nrow(cohort)
y <- as.numeric(cohort[[settings$outcome]])   # 1 = SRET at 32, 0 = not
w <- cohort$w9_final_weight                   # Wave-9 final weight

# Schools: the primary sampling unit within its stratum.
school <- interaction(cohort$stratum, cohort$psu, drop = TRUE)

models <- c("elastic_net", "logistic", "null")
model_labels <- c(elastic_net = "Elastic net",
                  logistic    = "Logistic regression (unpenalized)",
                  null        = "Null model (training prevalence)")


###############################################################################
# Block 3. Predictions in wide form
#
# One row per person (in the cohort's order) and one column per analysis x
# repeat x landmark x model, named e.g. "primary | 1 | L14 | elastic_net".
# Also, for each landmark, whether every feature was observed for the
# person (no feature or routed part imputed).
###############################################################################

splits <- unique(predictions[, c("strategy", "repeat_number", "landmark")])

prediction_matrix <- matrix(nrow = n, ncol = 0)
all_observed <- matrix(NA, nrow = n, ncol = length(landmark_order),
                       dimnames = list(NULL, landmark_order))

for (k in seq_len(nrow(splits))) {
  rows <- predictions$strategy == splits$strategy[k] &
    predictions$repeat_number == splits$repeat_number[k] &
    predictions$landmark == splits$landmark[k]
  split_predictions <- predictions[rows, ]
  position <- match(cohort$NSID, split_predictions$NSID)   # put into the cohort's order

  # Check: everyone is predicted once, with the cohort's outcome.
  if (anyNA(position) || !identical(as.numeric(split_predictions$outcome[position]), y)) {
    stop("The predictions of ", paste(splits[k, ], collapse = ", "),
         " don't match the cohort.")
  }

  for (model in models) {
    prediction_matrix <- cbind(prediction_matrix, split_predictions[[paste0("p_", model)]][position])
    colnames(prediction_matrix)[ncol(prediction_matrix)] <-
      paste(splits$strategy[k], splits$repeat_number[k], splits$landmark[k], model, sep = " | ")
  }

  # Whether every feature was observed depends only on the data, so it
  # must be the same in every repeat and analysis.
  observed_here <- split_predictions$all_features_observed[position]
  landmark <- splits$landmark[k]
  if (all(is.na(all_observed[, landmark]))) {
    all_observed[, landmark] <- observed_here
  } else if (!identical(all_observed[, landmark], observed_here)) {
    stop("Which features were observed differs between repeats at ", landmark, ".")
  }
}


###############################################################################
# Block 3b. Shuffled predictions for the comparison with File 20
#
# File 17a saved, for every outer test fold, the elastic net's predictions
# with each feature and each block of features shuffled among the test rows.
# File 20 compares the imputation analysis with these units, so they get the
# same bootstrap as the models:
#   - at every landmark: the science and mathematics orientation block and
#     each of the three science items;
#   - at each later landmark: the block it adds, and at L20 and L25a the
#     single feature that matches the imputation analysis's stage (HE
#     pathway, degree pathway).
# Primary analysis, every repeat; each column is assembled from the five
# test folds of a repeat, like the unshuffled predictions.
###############################################################################

orientation_units <- c("Block: Science and mathematics orientation", "Feature: science_liking14",
                       "Feature: science_competence14", "Feature: science_favourite14")
bridge_units <- rbind(
  expand.grid(landmark = landmark_order, unit = orientation_units, stringsAsFactors = FALSE),
  data.frame(landmark = c("L16", "L20", "L25a", "L25b", "L20", "L25a"),
             unit = c("Block: Course route and intentions", "Block: HE and early activity",
                      "Block: Attained education and activity",
                      "Block: Occupation at 25 (sensitivity)",
                      "Feature: he_pathway20", "Feature: degree_pathway25"))
)
bridge_units$model <- paste("shuffled", bridge_units$unit)
task_table <- cv_results$task_table

# One column per unit x repeat, each filled from the repeat's five test folds.
bridge_columns <- expand.grid(repeat_number = seq_len(settings$n_repeats),
                              unit_row = seq_len(nrow(bridge_units)))
for (k in seq_len(nrow(bridge_columns))) {
  r        <- bridge_columns$repeat_number[k]
  unit     <- bridge_units$unit[bridge_columns$unit_row[k]]
  landmark <- bridge_units$landmark[bridge_columns$unit_row[k]]
  column   <- rep(NA, n)
  tasks_r  <- which(task_table$strategy == "primary" & task_table$repeat_number == r &
                      task_table$landmark == landmark)
  for (t in tasks_r) {
    permuted <- cv_results$permuted[[t]]
    if (!(unit %in% colnames(permuted$predictions))) {
      stop("Task ", t, " has no shuffled predictions for ", unit, ".")
    }
    column[match(permuted$NSID, cohort$NSID)] <- permuted$predictions[, unit]
  }
  if (anyNA(column)) {
    stop("The shuffled predictions for ", unit, " (repeat ", r, ") don't cover the cohort.")
  }
  prediction_matrix <- cbind(prediction_matrix, column)
  colnames(prediction_matrix)[ncol(prediction_matrix)] <-
    paste("primary", r, landmark, bridge_units$model[bridge_columns$unit_row[k]], sep = " | ")
}
model_labels[bridge_units$model] <- paste0("Elastic net, ", bridge_units$unit, " shuffled")


###############################################################################
# Block 4. Groups for the subgroup results
#
# Sex, ethnic group and family NS-SEC at 14 use the observed values. People
# whose value is missing form a "Not known" group (their imputed values
# differ between folds and imputations). Groups are set before any result
# is seen; a group with fewer than 10 events or 10 non-events is suppressed.
###############################################################################

sex_group <- as.character(cohort$sex14)
sex_group[is.na(sex_group)] <- "Not known"

# Ethnic group: the deposited eight categories in five broad groups.
ethnic_groups <- c(
  "White"           = "White",
  "Mixed"           = "Mixed",
  "Indian"          = "Asian or Asian British",
  "Pakistani"       = "Asian or Asian British",
  "Bangladeshi"     = "Asian or Asian British",
  "Black Caribbean" = "Black or Black British",
  "Black African"   = "Black or Black British",
  "Other"           = "Other ethnic group"
)
not_grouped <- setdiff(levels(cohort$ethnicity14), names(ethnic_groups))
if (length(not_grouped) > 0) {
  stop("Ethnic-group categories without a broad group: ", paste(not_grouped, collapse = ", "),
       ". Add them to ethnic_groups.")
}
ethnic_group <- unname(ethnic_groups[as.character(cohort$ethnicity14)])
ethnic_group[is.na(ethnic_group)] <- "Not known"

# Family NS-SEC: the ONS three-class version, by class number. The
# categories are in class order (1-8), then "Not applicable" (no
# classifiable family occupation), which joins class 8.
nssec_levels <- levels(cohort$nssec_family14)
if (length(nssec_levels) != 9 || nssec_levels[9] != "Not applicable") {
  stop("Family NS-SEC doesn't have eight classes and 'Not applicable': ",
       paste(nssec_levels, collapse = ", "))
}
nssec_class <- as.integer(cohort$nssec_family14)   # 1-8; 9 = not applicable
nssec_group <- rep("Not known", n)
nssec_group[nssec_class %in% 1:2] <- "Managerial and professional (classes 1-2)"
nssec_group[nssec_class %in% 3:4] <- "Intermediate (classes 3-4)"
nssec_group[nssec_class %in% 5:7] <- "Routine and manual (classes 5-7)"
nssec_group[nssec_class %in% 8:9] <- "Never worked, long-term unemployed or not applicable"

group_values <- data.frame(check.names = FALSE,
                           "Sex"                 = sex_group,
                           "Ethnic group"        = ethnic_group,
                           "Family NS-SEC at 14" = nssec_group)


###############################################################################
# Block 5. The performance measures
#
# Written once and used for every result below. y: outcomes (0/1); p:
# predicted probabilities; w: weights (the Wave-9 final weight, or 1 for
# everyone when unweighted). Returns the seven measures described in the
# header.
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

measure_names <- c("prevalence", "pr_auc", "roc_auc", "brier", "log_loss",
                   "calibration_intercept", "calibration_slope")

# Checks of the measures on small made-up examples (not study data).
example_y <- c(0, 0, 1, 0, 1, 0, 0, 1, 0, 0)
example_w <- rep(1, 10)

perfect <- performance_measures(example_y, ifelse(example_y == 1, 0.9, 0.1), example_w)
if (perfect["pr_auc"] != 1 || perfect["roc_auc"] != 1) {
  stop("Measure check failed: a perfect ranking should give PR-AUC and ROC-AUC of 1.")
}
reversed <- performance_measures(example_y, ifelse(example_y == 1, 0.1, 0.9), example_w)
if (reversed["roc_auc"] != 0) stop("Measure check failed: a reversed ranking should give ROC-AUC 0.")
constant <- performance_measures(example_y, rep(0.3, 10), example_w)
if (abs(constant["pr_auc"] - 0.3) > 1e-12 || constant["roc_auc"] != 0.5 ||
    abs(constant["calibration_intercept"]) > 1e-6) {
  stop("Measure check failed: a constant prediction equal to the prevalence should give ",
       "PR-AUC = prevalence, ROC-AUC 0.5 and calibration intercept 0.")
}
example_p <- c(0.05, 0.40, 0.70, 0.10, 0.40, 0.20, 0.15, 0.30, 0.40, 0.02)
pairs_above <- outer(example_p[example_y == 1], example_p[example_y == 0], ">")
pairs_tied  <- outer(example_p[example_y == 1], example_p[example_y == 0], "==")
if (abs(performance_measures(example_y, example_p, example_w)["roc_auc"] -
        mean(pairs_above + pairs_tied / 2)) > 1e-12) {
  stop("Measure check failed: ROC-AUC differs from the count over all event/non-event pairs.")
}
# A weight of 2 must act like counting that person twice.
doubled <- performance_measures(example_y, example_p, c(2, rep(1, 9)))
repeated <- performance_measures(c(example_y, example_y[1]), c(example_p, example_p[1]), rep(1, 11))
if (!isTRUE(all.equal(doubled, repeated, tolerance = 1e-6))) {
  stop("Measure check failed: a weight of 2 should equal counting a person twice.")
}


###############################################################################
# Block 6. The table of results to estimate
#
# One row per result: analysis x landmark x model x weighting x group.
#   Analyses: "Primary" (repeats 1-5 averaged), and "Primary, repeat 1" and
#   "Outcome omitted, repeat 1" for the prespecified comparison.
# For each row, the people in the group and the unweighted events and
# non-events (distinct people), which decide whether the row is shown.
###############################################################################

analyses <- data.frame(
  analysis    = c("Primary", "Primary, repeat 1", "Outcome omitted, repeat 1"),
  strategy    = c("primary", "primary", "outcome omitted"),
  repeat_from = c(1, 1, 1),
  repeat_to   = c(settings$n_repeats, 1, 1)
)

estimates <- rbind(
  # 1. Every model, the whole cohort: primary weighted and unweighted; the
  #    repeat-1 analyses weighted.
  expand.grid(analysis = "Primary", landmark = landmark_order, model = models,
              weighting = c("Weighted", "Unweighted"), group_variable = "All", group = "All",
              stringsAsFactors = FALSE),
  expand.grid(analysis = analyses$analysis[2:3], landmark = landmark_order, model = models,
              weighting = "Weighted", group_variable = "All", group = "All",
              stringsAsFactors = FALSE),
  # 2. Elastic net, by whether every feature was observed at test time.
  expand.grid(analysis = analyses$analysis, landmark = landmark_order, model = "elastic_net",
              weighting = "Weighted", group_variable = "Features at test time",
              group = c("All features observed", "Any feature imputed"),
              stringsAsFactors = FALSE),
  # 3. Elastic net, primary analysis, by sex, ethnic group and NS-SEC.
  expand.grid(analysis = "Primary", landmark = landmark_order, model = "elastic_net",
              weighting = "Weighted", group_variable = "Sex", group = sort(unique(sex_group)),
              stringsAsFactors = FALSE),
  expand.grid(analysis = "Primary", landmark = landmark_order, model = "elastic_net",
              weighting = "Weighted", group_variable = "Ethnic group",
              group = sort(unique(ethnic_group)), stringsAsFactors = FALSE),
  expand.grid(analysis = "Primary", landmark = landmark_order, model = "elastic_net",
              weighting = "Weighted", group_variable = "Family NS-SEC at 14",
              group = sort(unique(nssec_group)), stringsAsFactors = FALSE),
  # 4. For File 20: the elastic net with each unit of Block 3b shuffled
  #    (primary analysis, weighted, whole cohort).
  data.frame(analysis = "Primary", landmark = bridge_units$landmark, model = bridge_units$model,
             weighting = "Weighted", group_variable = "All", group = "All")
)
estimates$key <- paste(estimates$analysis, estimates$landmark, estimates$model,
                       estimates$weighting, estimates$group_variable, estimates$group, sep = " | ")
estimates$strategy    <- analyses$strategy[match(estimates$analysis, analyses$analysis)]
estimates$repeat_from <- analyses$repeat_from[match(estimates$analysis, analyses$analysis)]
estimates$repeat_to   <- analyses$repeat_to[match(estimates$analysis, analyses$analysis)]

# Who is in each row's group (one column per row of the table).
members <- matrix(FALSE, nrow = n, ncol = nrow(estimates))
for (k in seq_len(nrow(estimates))) {
  if (estimates$group_variable[k] == "All") {
    members[, k] <- TRUE
  }
  if (estimates$group_variable[k] == "Features at test time") {
    members[, k] <- all_observed[, estimates$landmark[k]] ==
      (estimates$group[k] == "All features observed")
  }
  if (estimates$group_variable[k] %in% names(group_values)) {
    members[, k] <- group_values[[estimates$group_variable[k]]] == estimates$group[k]
  }
}

# Disclosure check: a row is shown only with at least 10 events and 10
# non-events (unweighted people). The counts stay in the private table.
estimates$events     <- colSums(members * y)
estimates$non_events <- colSums(members * (1 - y))
estimates$shown      <- estimates$events >= disclosure_threshold &
  estimates$non_events >= disclosure_threshold


###############################################################################
# Block 7. Estimates and bootstrap replicates
#
# Replicate 0 is the cohort itself and gives the estimates. Replicates 1 to
# bootstrap_replicates resample schools with replacement within each
# stratum (as many schools as the stratum has). A stratum with a single
# school always keeps it, so it adds no variation. Every row of the table
# is computed on the same resampled people, which pairs the comparisons.
###############################################################################

# The schools of each stratum, and the people of each school.
stratum_schools <- list()
for (s in levels(cohort$stratum)) {
  stratum_schools[[s]] <- unique(as.character(school[cohort$stratum == s]))
}
school_rows <- list()
for (school_name in levels(school)) {
  school_rows[[school_name]] <- which(school == school_name)
}

results <- array(NA, dim = c(bootstrap_replicates + 1, nrow(estimates), length(measure_names)),
                 dimnames = list(paste0("replicate_", 0:bootstrap_replicates),
                                 estimates$key, measure_names))
rows_to_compute <- which(estimates$shown)

cat("\nBootstrap:", bootstrap_replicates, "replicates of", length(stratum_schools),
    "strata and", length(school_rows), "schools.\n")
start_time <- Sys.time()

for (b in 0:bootstrap_replicates) {
  # The people in this replicate (a person appears once per time their
  # school is drawn).
  if (b == 0) {
    rows_b <- seq_len(n)
  } else {
    set.seed(bootstrap_seed + b)
    drawn_schools <- c()
    for (s in names(stratum_schools)) {
      schools_s <- stratum_schools[[s]]
      drawn_schools <- c(drawn_schools, schools_s[sample.int(length(schools_s), replace = TRUE)])
    }
    rows_b <- unlist(school_rows[drawn_schools], use.names = FALSE)
  }

  # Every shown row: the measures in each of its repeats, averaged.
  for (k in rows_to_compute) {
    rows_k <- rows_b[members[rows_b, k]]
    weights_k <- if (estimates$weighting[k] == "Weighted") w[rows_k] else rep(1, length(rows_k))
    repeats <- estimates$repeat_from[k]:estimates$repeat_to[k]
    measure_sum <- 0
    for (r in repeats) {
      column <- paste(estimates$strategy[k], r, estimates$landmark[k], estimates$model[k], sep = " | ")
      measure_sum <- measure_sum + performance_measures(y[rows_k], prediction_matrix[rows_k, column], weights_k)
    }
    results[b + 1, k, ] <- measure_sum / length(repeats)
  }

  if (b %% 50 == 0 && b > 0) {
    minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat("  Replicate", b, "of", bootstrap_replicates, "done;", round(minutes, 1),
        "minutes so far, about", round(minutes / b * (bootstrap_replicates - b)), "to go.\n")
  }
}
bootstrap_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

# The null model predicts each training set's prevalence, so its
# predictions barely vary and its calibration slope has no meaning.
results[, estimates$model == "null", "calibration_slope"] <- NA


###############################################################################
# Block 8. Performance table
#
# Estimate (replicate 0) and 95% percentile interval (replicates 1 onward).
# replicates_used counts the replicates in which the measure could be
# computed (e.g. a calibration fit that didn't converge is left out).
###############################################################################

performance <- data.frame()
for (k in seq_len(nrow(estimates))) {
  for (m in measure_names) {
    replicate_values <- results[-1, k, m]
    performance <- rbind(performance, data.frame(
      analysis = estimates$analysis[k], landmark = estimates$landmark[k],
      model = model_labels[estimates$model[k]], weighting = estimates$weighting[k],
      group_variable = estimates$group_variable[k], group = estimates$group[k],
      measure = m,
      estimate = results[1, k, m],
      ci_lower = unname(quantile(replicate_values, 0.025, na.rm = TRUE)),
      ci_upper = unname(quantile(replicate_values, 0.975, na.rm = TRUE)),
      replicates_used = sum(!is.na(replicate_values)),
      note = if (estimates$shown[k]) "" else "Suppressed: fewer than 10 events or 10 non-events"
    ))
  }
}
rownames(performance) <- NULL


###############################################################################
# Block 9. Paired comparisons
#
# Each comparison is first minus second, computed in every replicate on the
# same people, so its interval is paired.
#   a. Each landmark minus the one before (same analysis, model and
#      weighting): what the later landmark adds.
#   b. At each landmark: elastic net minus logistic regression, and each of
#      them minus the null model.
#   c. Primary minus outcome-omitted analysis, on repeat 1 (elastic net and
#      logistic regression, weighted).
#   d. For File 20: the elastic net minus the elastic net with a unit
#      shuffled (permutation importance: a positive PR-AUC difference means
#      the model relied on the shuffled information).
###############################################################################

overall <- estimates[estimates$group_variable == "All", ]
comparison_pairs <- data.frame()

# a. Later landmark minus earlier landmark.
for (k in which(overall$model %in% c("elastic_net", "logistic") &
                overall$landmark != landmark_order[1])) {
  earlier <- landmark_order[match(overall$landmark[k], landmark_order) - 1]
  comparison_pairs <- rbind(comparison_pairs, data.frame(
    comparison = "Landmark minus the landmark before", analysis = overall$analysis[k],
    weighting = overall$weighting[k], landmark = overall$landmark[k],
    first = paste(overall$landmark[k], model_labels[overall$model[k]]),
    second = paste(earlier, model_labels[overall$model[k]]),
    first_key = overall$key[k],
    second_key = paste(overall$analysis[k], earlier, overall$model[k], overall$weighting[k],
                       "All", "All", sep = " | ")))
}

# b. Models against each other at the same landmark.
model_pairs <- data.frame(first  = c("elastic_net", "elastic_net", "logistic"),
                          second = c("logistic", "null", "null"))
for (k in which(overall$model == "elastic_net")) {
  for (j in seq_len(nrow(model_pairs))) {
    comparison_pairs <- rbind(comparison_pairs, data.frame(
      comparison = "Model minus model", analysis = overall$analysis[k],
      weighting = overall$weighting[k], landmark = overall$landmark[k],
      first = model_labels[model_pairs$first[j]], second = model_labels[model_pairs$second[j]],
      first_key = paste(overall$analysis[k], overall$landmark[k], model_pairs$first[j],
                        overall$weighting[k], "All", "All", sep = " | "),
      second_key = paste(overall$analysis[k], overall$landmark[k], model_pairs$second[j],
                         overall$weighting[k], "All", "All", sep = " | ")))
  }
}

# c. Primary minus outcome omitted (repeat 1).
for (k in which(overall$analysis == "Primary, repeat 1" &
                overall$model %in% c("elastic_net", "logistic"))) {
  comparison_pairs <- rbind(comparison_pairs, data.frame(
    comparison = "Primary minus outcome omitted (repeat 1)", analysis = "Repeat 1",
    weighting = overall$weighting[k], landmark = overall$landmark[k],
    first = paste("Primary,", model_labels[overall$model[k]]),
    second = paste("Outcome omitted,", model_labels[overall$model[k]]),
    first_key = overall$key[k],
    second_key = paste("Outcome omitted, repeat 1", overall$landmark[k], overall$model[k],
                       overall$weighting[k], "All", "All", sep = " | ")))
}

# d. Elastic net minus elastic net with a unit shuffled (for File 20).
for (k in seq_len(nrow(bridge_units))) {
  comparison_pairs <- rbind(comparison_pairs, data.frame(
    comparison = "Elastic net minus elastic net with a unit shuffled", analysis = "Primary",
    weighting = "Weighted", landmark = bridge_units$landmark[k],
    first = model_labels["elastic_net"], second = model_labels[bridge_units$model[k]],
    first_key = paste("Primary", bridge_units$landmark[k], "elastic_net", "Weighted",
                      "All", "All", sep = " | "),
    second_key = paste("Primary", bridge_units$landmark[k], bridge_units$model[k], "Weighted",
                       "All", "All", sep = " | ")))
}

# The differences, with paired 95% intervals.
comparison_measures <- setdiff(measure_names, "prevalence")   # same people, so prevalence doesn't differ
comparisons <- data.frame()
for (k in seq_len(nrow(comparison_pairs))) {
  for (m in comparison_measures) {
    difference <- results[, comparison_pairs$first_key[k], m] - results[, comparison_pairs$second_key[k], m]
    comparisons <- rbind(comparisons, data.frame(
      comparison_pairs[k, c("comparison", "analysis", "weighting", "landmark", "first", "second")],
      measure = m,
      difference = difference[1],
      ci_lower = unname(quantile(difference[-1], 0.025, na.rm = TRUE)),
      ci_upper = unname(quantile(difference[-1], 0.975, na.rm = TRUE)),
      replicates_used = sum(!is.na(difference[-1]))
    ))
  }
}
rownames(comparisons) <- NULL


###############################################################################
# Block 10. Primary against outcome omitted, fold by fold (repeat 1)
#
# The prespecified pairing by fold x landmark: both analyses share repeat
# 1's folds. Elastic net, weighted. Every outer test fold has at least 10
# events and 10 non-events (File 14); checked again here.
###############################################################################

fold_table <- expand.grid(outer_fold = seq_len(settings$n_outer_folds), landmark = landmark_order,
                          stringsAsFactors = FALSE)
fold_measures <- c("pr_auc", "roc_auc", "brier", "calibration_intercept", "calibration_slope")
fold_comparison <- data.frame()

for (k in seq_len(nrow(fold_table))) {
  in_fold <- cohort$outer_fold_r1 == fold_table$outer_fold[k]
  if (sum(y[in_fold]) < disclosure_threshold || sum(1 - y[in_fold]) < disclosure_threshold) {
    stop("Outer fold ", fold_table$outer_fold[k], " of repeat 1 has fewer than 10 events or non-events.")
  }
  primary <- performance_measures(
    y[in_fold], prediction_matrix[in_fold, paste("primary", 1, fold_table$landmark[k], "elastic_net", sep = " | ")],
    w[in_fold])
  omitted <- performance_measures(
    y[in_fold], prediction_matrix[in_fold, paste("outcome omitted", 1, fold_table$landmark[k], "elastic_net", sep = " | ")],
    w[in_fold])
  fold_comparison <- rbind(fold_comparison, data.frame(
    landmark = fold_table$landmark[k], outer_fold = fold_table$outer_fold[k],
    measure = fold_measures, primary = primary[fold_measures],
    outcome_omitted = omitted[fold_measures],
    difference = primary[fold_measures] - omitted[fold_measures]))
}
rownames(fold_comparison) <- NULL


###############################################################################
# Block 11. Binned calibration (repeat 1)
#
# For each analysis x landmark x model (elastic net and logistic
# regression), people are split into equal-sized bins by predicted
# probability. Going up from the lowest bin, a bin with fewer than 10 events
# or 10 non-events is merged into the next one; the top bin, if still
# short, is merged into the one below. Each bin gives the (weighted) mean
# prediction and observed proportion. Repeat 1 is used because both
# analyses share it; each person appears once.
###############################################################################

calibration_sets <- expand.grid(landmark = landmark_order, model = c("elastic_net", "logistic"),
                                strategy = c("primary", "outcome omitted"),
                                stringsAsFactors = FALSE)
calibration_table <- data.frame()

for (k in seq_len(nrow(calibration_sets))) {
  p <- prediction_matrix[, paste(calibration_sets$strategy[k], 1, calibration_sets$landmark[k],
                                 calibration_sets$model[k], sep = " | ")]

  # Equal-sized starting bins.
  breaks <- unique(quantile(p, probs = seq(0, 1, length.out = calibration_bins + 1)))
  bin <- cut(p, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  # Merge short bins upwards, then the top bin downwards if needed.
  for (b in sort(unique(bin))) {
    in_bin <- bin == b
    too_small <- sum(y[in_bin]) < disclosure_threshold || sum(1 - y[in_bin]) < disclosure_threshold
    if (too_small && b < max(bin)) bin[in_bin] <- b + 1
  }
  top <- max(bin)
  too_small <- sum(y[bin == top]) < disclosure_threshold || sum(1 - y[bin == top]) < disclosure_threshold
  if (too_small && length(unique(bin)) > 1) bin[bin == top] <- max(bin[bin != top])

  # Disclosure check and the bin summaries.
  bins <- sort(unique(bin))
  for (b in bins) {
    in_bin <- bin == b
    if (sum(y[in_bin]) < disclosure_threshold || sum(1 - y[in_bin]) < disclosure_threshold) {
      stop("A calibration bin has fewer than 10 events or non-events.")
    }
    calibration_table <- rbind(calibration_table, data.frame(
      analysis = if (calibration_sets$strategy[k] == "primary") "Primary, repeat 1" else "Outcome omitted, repeat 1",
      landmark = calibration_sets$landmark[k], model = model_labels[calibration_sets$model[k]],
      bin = match(b, bins),
      mean_predicted = weighted.mean(p[in_bin], w[in_bin]),
      observed_proportion = weighted.mean(y[in_bin], w[in_bin])))
  }
}
rownames(calibration_table) <- NULL

# Calibration plots: bins only (no individual points or rug marks).
pdf(file.path(review_dir, "18_calibration_plots.pdf"), width = 11, height = 7.5)
for (model in c("elastic_net", "logistic")) {
  par(mfrow = c(2, 3))
  for (landmark in landmark_order) {
    primary_bins <- calibration_table[calibration_table$model == model_labels[model] &
                                        calibration_table$landmark == landmark &
                                        calibration_table$analysis == "Primary, repeat 1", ]
    omitted_bins <- calibration_table[calibration_table$model == model_labels[model] &
                                        calibration_table$landmark == landmark &
                                        calibration_table$analysis == "Outcome omitted, repeat 1", ]
    axis_end <- 1.05 * max(c(primary_bins$mean_predicted, primary_bins$observed_proportion,
                             omitted_bins$mean_predicted, omitted_bins$observed_proportion))
    plot(primary_bins$mean_predicted, primary_bins$observed_proportion, type = "b", pch = 19,
         xlim = c(0, axis_end), ylim = c(0, axis_end),
         xlab = "Mean predicted probability in bin", ylab = "Observed proportion in bin (weighted)",
         main = paste(landmark, "-", model_labels[model]))
    lines(omitted_bins$mean_predicted, omitted_bins$observed_proportion, type = "b", pch = 2, lty = 3)
    abline(0, 1, col = "grey50", lty = 2)
    legend("topleft", legend = c("Primary", "Outcome omitted", "Perfect calibration"),
           pch = c(19, 2, NA), lty = c(1, 3, 2), col = c("black", "black", "grey50"), bty = "n")
  }
  plot.new()
  text(0.5, 0.5, paste("Out-of-fold predictions of repeat 1, in bins of predicted probability.",
                       "Each bin has at least 10 events and 10 non-events.",
                       "Weighted by the Wave-9 final weight.", sep = "\n"), cex = 0.9)
}
invisible(dev.off())


###############################################################################
# Block 12. Differences between schools, in aggregate
#
# A random-intercept logistic model of the outcome, with the log-odds of the
# prediction as an offset and a random intercept for each school, gives the
# standard deviation of school-level calibration intercepts: how much
# schools' observed outcomes vary around what the model predicts. For the
# null model it is the variation in the outcome between schools; the
# elastic net's value shows how much of it remains. Primary analysis,
# unweighted, each repeat, averaged. Only the aggregate is written.
###############################################################################

school_fits <- expand.grid(repeat_number = seq_len(settings$n_repeats), model = c("null", "elastic_net"),
                           landmark = landmark_order, stringsAsFactors = FALSE)
school_fits$school_sd <- NA
school_fits$singular <- NA
school_fits$convergence_messages <- NA

for (k in seq_len(nrow(school_fits))) {
  p <- prediction_matrix[, paste("primary", school_fits$repeat_number[k], school_fits$landmark[k],
                                 school_fits$model[k], sep = " | ")]
  school_data <- data.frame(y = y, school = school, log_odds = qlogis(p))
  # A singular fit (school SD estimated as 0) is recorded in the table
  # instead of printed.
  fit <- suppressMessages(glmer(y ~ 1 + offset(log_odds) + (1 | school), data = school_data,
                                family = binomial))
  school_fits$school_sd[k] <- attr(VarCorr(fit)$school, "stddev")
  school_fits$singular[k] <- isSingular(fit)
  # lme4 stores its convergence messages here; the singular-fit notice is
  # already counted above, so it isn't counted again.
  messages <- fit@optinfo$conv$lme4$messages
  school_fits$convergence_messages[k] <- sum(!grepl("singular", messages))
}

school_differences <- data.frame()
for (landmark in landmark_order) {
  for (model in c("null", "elastic_net")) {
    rows <- school_fits$landmark == landmark & school_fits$model == model
    school_differences <- rbind(school_differences, data.frame(
      landmark = landmark, model = model_labels[model],
      school_sd_log_odds = mean(school_fits$school_sd[rows]),
      repeats = sum(rows),
      singular_fits = sum(school_fits$singular[rows]),
      fits_with_convergence_messages = sum(school_fits$convergence_messages[rows] > 0),
      schools = nlevels(school)))
  }
}
rownames(school_differences) <- NULL


###############################################################################
# Block 13. Save and write the review tables
#
# Disclosure check (before writing):
#   - performance: rows of groups with fewer than 10 events or 10
#     non-events have no values (Block 6); no counts of people are written;
#   - comparisons: whole-cohort rows only (every group shown);
#   - fold by fold: every fold has at least 10 events and 10 non-events
#     (checked in Block 10);
#   - calibration bins: every bin has at least 10 events and 10 non-events
#     (checked in Block 11); no individual points;
#   - school differences: models of the whole cohort with more than 10
#     schools; aggregate only.
###############################################################################

if (nlevels(school) < disclosure_threshold) stop("Fewer than 10 schools; the school model can't be shown.")

saveRDS(list(estimates = estimates, results = results, performance = performance,
             comparisons = comparisons, bridge_units = bridge_units,
             fold_comparison = fold_comparison,
             calibration_table = calibration_table, school_fits = school_fits,
             settings = list(bootstrap_replicates = bootstrap_replicates,
                             bootstrap_seed = bootstrap_seed, measure_names = measure_names)),
        file.path(data_dir, "18_performance.rds"))

# Exact counts behind each row (private; never copied for review).
write.csv(estimates[, c("analysis", "landmark", "model", "weighting", "group_variable", "group",
                        "events", "non_events", "shown")],
          file.path(output_dir, "private_18_group_counts.csv"), row.names = FALSE)

performance_review <- performance
performance_review[performance_review$note != "", c("estimate", "ci_lower", "ci_upper")] <- NA
performance_review[, c("estimate", "ci_lower", "ci_upper")] <- round(performance_review[, c("estimate", "ci_lower", "ci_upper")], 5)
performance_review$disclosure_checked <- TRUE
write.csv(performance_review, file.path(review_dir, "18_performance.csv"), row.names = FALSE)

comparisons_review <- comparisons
comparisons_review[, c("difference", "ci_lower", "ci_upper")] <- round(comparisons_review[, c("difference", "ci_lower", "ci_upper")], 5)
comparisons_review$disclosure_checked <- TRUE
write.csv(comparisons_review, file.path(review_dir, "18_comparisons.csv"), row.names = FALSE)

fold_review <- fold_comparison
fold_review[, c("primary", "outcome_omitted", "difference")] <- round(fold_review[, c("primary", "outcome_omitted", "difference")], 5)
fold_review$disclosure_checked <- TRUE
write.csv(fold_review, file.path(review_dir, "18_strategy_comparison_by_fold.csv"), row.names = FALSE)

calibration_review <- calibration_table
calibration_review[, c("mean_predicted", "observed_proportion")] <- round(calibration_review[, c("mean_predicted", "observed_proportion")], 4)
calibration_review$disclosure_checked <- TRUE
write.csv(calibration_review, file.path(review_dir, "18_calibration_bins.csv"), row.names = FALSE)

school_review <- school_differences
school_review$school_sd_log_odds <- round(school_review$school_sd_log_odds, 3)
school_review$disclosure_checked <- TRUE
write.csv(school_review, file.path(review_dir, "18_school_differences.csv"), row.names = FALSE)

disclosure_log <- data.frame(
  table = c("18_performance.csv", "18_comparisons.csv", "18_strategy_comparison_by_fold.csv",
            "18_calibration_bins.csv", "18_calibration_plots.pdf", "18_school_differences.csv"),
  rule = c("Groups with fewer than 10 events or 10 non-events suppressed; no counts of people",
           "Whole-cohort comparisons only",
           "Every outer fold has at least 10 events and 10 non-events",
           "Bins merged until each has at least 10 events and 10 non-events; no individual points",
           "Drawn from the checked bins only",
           "Whole-cohort models; aggregate over schools only"),
  rows_suppressed = c(sum(!estimates$shown), 0, 0, 0, 0, 0)
)
disclosure_log$disclosure_checked <- TRUE
write.csv(disclosure_log, file.path(review_dir, "18_disclosure_log.csv"), row.names = FALSE)


###############################################################################
# Block 14. Summary
#
# Counts and run time only. The results are in the review tables.
###############################################################################

suppressed_groups <- estimates[!estimates$shown, ]
calibration_values <- results[-1, estimates$shown & estimates$model != "null", "calibration_slope"]

cat("\nFile 18 finished.\n")
cat("  Measure checks on made-up examples:     passed\n")
cat("  Results estimated (rows):              ", sum(estimates$shown), "of", nrow(estimates), "\n")
cat("  Rows suppressed (fewer than 10 events or non-events):", nrow(suppressed_groups), "\n")
if (nrow(suppressed_groups) > 0) {
  cat("    Groups suppressed:", paste(unique(paste(suppressed_groups$group_variable, suppressed_groups$group,
                                                  sep = ": ")), collapse = "; "), "\n")
}
cat("  Bootstrap replicates:                  ", bootstrap_replicates, "(",
    round(bootstrap_minutes, 1), "minutes )\n")
cat("  Calibration slopes not estimable in a replicate:", sum(is.na(calibration_values)), "of",
    length(calibration_values), "\n")
cat("  Comparisons:                           ", nrow(comparison_pairs), "\n")
cat("  Calibration bins per plot (fewest):    ", min(table(paste(calibration_table$analysis,
                                                               calibration_table$landmark,
                                                               calibration_table$model))), "\n")
cat("  School models: singular fits", sum(school_fits$singular), "; fits with convergence messages",
    sum(school_fits$convergence_messages > 0), "of", nrow(school_fits), "\n")
cat("  Review: 18_performance.csv, 18_comparisons.csv, 18_strategy_comparison_by_fold.csv,\n")
cat("    18_calibration_bins.csv, 18_calibration_plots.pdf, 18_school_differences.csv,\n")
cat("    18_disclosure_log.csv\n")
