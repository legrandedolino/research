###############################################################################
# 20. Association versus prediction
#
# Purpose
#   The study's integrating question: how much does adolescent science
#   orientation tell us about professional SRET employment at 32? First as a
#   population association that holds up under attrition (the imputation
#   analysis, Files 7-13), and second as information that tells individuals
#   apart (the prediction analysis, Files 14-19). This script sets the two
#   side by side, stage by stage. Neither is causal. Their magnitudes are on
#   different scales (percentage points and PR-AUC), so only their patterns
#   are compared.
#
# The stages
#   Stage                            Imputation model   Prediction landmark
#   Background + orientation at 14   M0                 L14
#   Course route at 15               M1                 (L16 adds 15 and 16 together)
#   Intentions at 16                 M2                 L16
#   HE by 20                         M3                 L20 (also activity 16-19)
#   Degree by 25                     M4                 L25a (also activity at 25)
#   SRET job at 25                   (none)             L25b (prediction only)
#
# Result 1: does the science signal persist across the stages?
#   For each science contrast (liking, perceived competence, favourite
#   subject): its standardized risk difference under M0-M4 (File 12, primary
#   analysis), beside the elastic net's permutation importance of the same
#   item, and of the whole orientation block, at each landmark (File 18).
#
# Result 2: where the association changes, does the stage add predictive
# value?
#   For each stage: the change in each science risk difference from the
#   model before (M0 to M2 for L16, M2 to M3 for L20, M3 to M4 for L25a),
#   beside the landmark's gain in PR-AUC over the landmark before and the
#   importance of the block it adds (File 18). M0 to M1, M1 to M2 and M0 to
#   M4 are also shown. Risk differences, not odds ratios, are compared across
#   models, because odds ratios change across nested logistic models even
#   without confounding (Mood 2010).
#
# Uncertainty
#   Risk differences: File 12's estimates and intervals.
#   Changes between models: File 12 fits each model separately, so a change
#   needs both models fitted to the same resampled data. In each completed
#   dataset, schools are resampled within strata (Rao and Wu 1988: in a
#   stratum with n schools, n - 1 are drawn with replacement, and weights are
#   multiplied by n / (n - 1) and by the number of times the school was
#   drawn; a stratum with one school keeps it). The same draws are used in
#   every completed dataset. M0-M4 are refitted and each change recomputed.
#   The within-dataset variance is the variance over the replicates, and
#   Rubin's rules combine the datasets as in File 12.
#   Prediction: File 18's paired bootstrap of schools.
#
# Rules, fixed before any result was seen
#   Association change: "weakened" (opposite in sign to the earlier
#     estimate) or "strengthened" if its 95% CI excludes 0; otherwise "no
#     clear change".
#   Prediction: "adds predictive value" if the landmark's PR-AUC gain has a
#     95% CI above 0; "lowers performance" if below 0; otherwise "no clear
#     added value". Block importance is shown beside it.
#   Agreement (L16, L20, L25a): "agree" if both lanes show something or
#     neither does; otherwise "differ".
#   Result 1: an association is "present" if the risk difference's CI
#     excludes 0; predictive information is "present" if the item's
#     importance CI is above 0.
#
# Populations
#   Association: the primary imputation analysis (target population, design
#   weight), with the observed-outcome analysis (the prediction cohort's
#   people, Wave-9 weight) beside it as point estimates. Prediction: the
#   elastic net, weighted by the Wave-9 final weight. All association
#   results assume missing at random (File 13 examines departures for M0).
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   Every quantity is a pooled or whole-cohort estimate. Risk differences
#   and changes are shown only for models with at least 10 events and 10
#   non-events in every completed dataset (File 12's rule).
#
# Run time
#   About 40 minutes (90 completed datasets x 5 models x 101 fits).
#
# Reads
#   data_dir/07_imputation_spec.rds, 08_imputations_main.rds,
#   12_pooled_estimates.rds, 14_prediction_design.rds, 18_performance.rds
#
# Writes
#   data_dir/20_bridge.rds
#   review_dir/20_orientation_across_stages.csv, 20_stage_changes.csv,
#     20_bridge_plots.pdf, 20_disclosure_log.csv
#
# References
#   Barnard J, Rubin DB (1999). Biometrika 86(4): 948-955.
#   Mood C (2010). Logistic regression: why we cannot do what we think we
#   can do, and what we can do about it. European Sociological Review 26(1):
#   67-82.
#   Rao JNK, Wu CFJ (1988). Resampling inference with complex survey data.
#   Journal of the American Statistical Association 83(401): 231-241.
#   Rubin DB (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
###############################################################################

suppressPackageStartupMessages(library(mice))


###############################################################################
# Block 1. Paths and settings
###############################################################################

paths      <- readRDS("project_paths.rds")
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs

disclosure_threshold <- 10    # fewest events (and non-events) behind a shown model
bootstrap_replicates <- 100   # school resamples per completed dataset
bootstrap_seed       <- 2000  # replicate r uses seed bootstrap_seed + r

primary_analysis  <- "Primary: target population, design weight"
observed_analysis <- "Sensitivity: observed-outcome sample, Wave-9 weight"


###############################################################################
# Block 2. Read the imputations, File 12's results and File 18's results
###############################################################################

imputation_spec <- readRDS(file.path(data_dir, "07_imputation_spec.rds"))
design_data     <- imputation_spec$design_data      # same row order as the imputations
imputations     <- readRDS(file.path(data_dir, "08_imputations_main.rds"))

file12        <- readRDS(file.path(data_dir, "12_pooled_estimates.rds"))
contrasts     <- file12$contrasts
models        <- file12$models[file12$models$model != "M0 baseline, employed only", ]   # M0-M4
by_imputation <- file12$results_by_imputation
risk_differences_12 <- file12$risk_differences

landmark_order <- readRDS(file.path(data_dir, "14_prediction_design.rds"))$settings$landmark_order
file18 <- readRDS(file.path(data_dir, "18_performance.rds"))

if (nrow(models) != 5 || any(models$outcome != "sret_professional_32")) {
  stop("File 12's models aren't the expected M0-M4 with the overall SRET outcome.")
}
n_imputations <- imputations$m


###############################################################################
# Block 3. The stages and the changes between models
###############################################################################

stages <- data.frame(
  stage    = c("Background and orientation at 14", "Course route at 15", "Intentions at 16",
               "HE by 20", "Degree by 25", "SRET job at 25"),
  model    = c(models$model, NA),
  landmark = c("L14", NA, "L16", "L20", "L25a", "L25b")
)

changes <- data.frame(
  change = c("Course route at 15 (M0 to M1)", "Intentions at 16 (M1 to M2)",
             "Course route and intentions (M0 to M2)", "HE by 20 (M2 to M3)",
             "Degree by 25 (M3 to M4)", "All later stages (M0 to M4)"),
  earlier  = models$model[c(1, 2, 1, 3, 4, 1)],
  later    = models$model[c(2, 3, 3, 4, 5, 5)],
  landmark = c(NA, NA, "L16", "L20", "L25a", NA)
)

# The permutation unit of each science contrast in the prediction analysis.
contrasts$unit <- paste("Feature:", contrasts$variable)


###############################################################################
# Block 4. School resamples (the same in every completed dataset)
#
# A school is a primary sampling unit within its stratum, as in File 12's
# survey design. multipliers[school, r] is what replicate r multiplies the
# design weight by (0 for a school not drawn).
###############################################################################

school  <- as.character(interaction(design_data$stratum, design_data$psu, drop = TRUE))
schools <- unique(data.frame(school = school, stratum = as.character(design_data$stratum)))

multipliers <- matrix(1, nrow = nrow(schools), ncol = bootstrap_replicates,
                      dimnames = list(schools$school, NULL))
for (r in seq_len(bootstrap_replicates)) {
  set.seed(bootstrap_seed + r)
  for (s in unique(schools$stratum)) {
    schools_s <- schools$school[schools$stratum == s]
    n_s <- length(schools_s)
    if (n_s < 2) next    # a stratum with one school keeps it (no variation)
    drawn <- schools_s[sample.int(n_s, n_s - 1, replace = TRUE)]
    times_drawn <- as.numeric(table(factor(drawn, levels = schools_s)))
    multipliers[schools_s, r] <- times_drawn * n_s / (n_s - 1)
  }
}
person_school <- match(school, rownames(multipliers))   # each person's row in multipliers


###############################################################################
# Block 5. Risk differences in every completed dataset and replicate
#
# One row of fit_table per completed dataset x model. For each:
#   a. the completed dataset, prepared as in File 12 (once per dataset);
#   b. the model's design matrix, and the matrices with each science item
#      set to its "high" and its "low" level (File 5's contrasts);
#   c. for replicate 0 (the design weight itself) and each school
#      replicate: a weighted logistic regression, and each contrast's
#      standardized risk difference, weighted by the same weights:
#        RD = weighted mean of p(high) - weighted mean of p(low).
#      Replicates start from replicate 0's coefficients, which only speeds
#      up the fit.
# rd_values[dataset, model, replicate, contrast] keeps every value.
###############################################################################

fit_table <- expand.grid(model = models$model, imputation = seq_len(n_imputations),
                         stringsAsFactors = FALSE)
rd_values <- array(NA, dim = c(n_imputations, nrow(models), bootstrap_replicates + 1, nrow(contrasts)),
                   dimnames = list(NULL, models$model, paste0("replicate_", 0:bootstrap_replicates),
                                   contrasts$contrast))
fits_not_converged <- 0

cat("\nRefitting M0-M4 in", n_imputations, "completed datasets x", bootstrap_replicates + 1,
    "weightings.\n")
start_time <- Sys.time()

for (k in seq_len(nrow(fit_table))) {
  i <- fit_table$imputation[k]
  m <- match(fit_table$model[k], models$model)

  # a. The completed dataset, as in File 12.
  if (m == 1) {
    analysis_data <- cbind(design_data, complete(imputations, i))
    analysis_data$he_likelihood16 <- factor(analysis_data$he_likelihood16, ordered = FALSE)
  }

  # b. Design matrices.
  model_formula   <- as.formula(paste(models$outcome[m], "~", models$predictors[m]))
  predictor_terms <- delete.response(terms(model_formula))
  X <- model.matrix(predictor_terms, data = analysis_data)
  y <- analysis_data[[models$outcome[m]]]
  if (nrow(X) != nrow(analysis_data)) stop("Rows were dropped from ", models$model[m], ".")
  X_high <- list()
  X_low  <- list()
  for (j in seq_len(nrow(contrasts))) {
    data_high <- analysis_data
    data_low  <- analysis_data
    variable  <- contrasts$variable[j]
    if (is.numeric(data_high[[variable]])) {
      data_high[[variable]] <- as.numeric(contrasts$high[j])
      data_low[[variable]]  <- as.numeric(contrasts$low[j])
    } else {
      data_high[[variable]] <- factor(contrasts$high[j], levels = levels(analysis_data[[variable]]))
      data_low[[variable]]  <- factor(contrasts$low[j],  levels = levels(analysis_data[[variable]]))
    }
    X_high[[j]] <- model.matrix(predictor_terms, data = data_high)
    X_low[[j]]  <- model.matrix(predictor_terms, data = data_low)
  }

  # c. Replicate 0 (design weight), then the school replicates.
  start_values <- NULL
  for (r in 0:bootstrap_replicates) {
    w <- design_data$design_weight
    if (r > 0) w <- w * multipliers[person_school, r]
    fit <- glm.fit(X, y, weights = w, family = quasibinomial(), start = start_values)
    if (!fit$converged) fits_not_converged <- fits_not_converged + 1
    beta <- fit$coefficients
    beta[is.na(beta)] <- 0            # a column for an empty category has no coefficient
    if (r == 0) start_values <- beta
    for (j in seq_len(nrow(contrasts))) {
      p_high <- plogis(as.vector(X_high[[j]] %*% beta))
      p_low  <- plogis(as.vector(X_low[[j]] %*% beta))
      rd_values[i, m, r + 1, j] <- sum(w * (p_high - p_low)) / sum(w)
    }
  }

  if (m == nrow(models) && i %% 10 == 0) {
    minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat("  Completed dataset", i, "of", n_imputations, "done;", round(minutes, 1),
        "minutes so far, about", round(minutes / i * (n_imputations - i)), "to go.\n")
  }
}
bootstrap_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))


###############################################################################
# Block 6. Checks
#
#   1. Replicate 0 reproduces File 12's risk difference in every completed
#      dataset, model and contrast (the script stops otherwise).
#   2. Sanity check (reported, not a stop): the bootstrap standard errors
#      of the risk differences against File 12's delta-method standard
#      errors, as a median ratio over datasets, models and contrasts.
###############################################################################

check_table <- expand.grid(imputation = seq_len(n_imputations), model = models$model,
                           contrast = contrasts$contrast, stringsAsFactors = FALSE)
se_ratio <- rep(NA, nrow(check_table))
for (k in seq_len(nrow(check_table))) {
  row_12 <- by_imputation$analysis == primary_analysis & by_imputation$model == check_table$model[k] &
    by_imputation$quantity == "risk difference" & by_imputation$term == check_table$contrast[k] &
    by_imputation$imputation == check_table$imputation[k]
  if (sum(row_12) != 1) stop("File 12 has no single risk difference for ", paste(check_table[k, ], collapse = ", "))
  values <- rd_values[check_table$imputation[k], check_table$model[k], , check_table$contrast[k]]
  if (abs(values[1] - by_imputation$estimate[row_12]) > 1e-6) {
    stop("The refit doesn't reproduce File 12's risk difference (",
         paste(check_table[k, ], collapse = ", "), ").")
  }
  se_ratio[k] <- sd(values[-1]) / sqrt(by_imputation$variance[row_12])
}
median_se_ratio <- median(se_ratio)


###############################################################################
# Block 7. The changes between models, with Rubin's rules
#
# In each completed dataset: the change (later model minus earlier model)
# at replicate 0, and its variance over the school replicates. Then, as in
# File 12: pooled estimate, within and between variance, Barnard-Rubin
# degrees of freedom with the survey design's degrees of freedom as the
# complete-data value, and a 95% CI.
###############################################################################

df_complete <- by_imputation$df_complete[by_imputation$analysis == primary_analysis &
                                           by_imputation$quantity == "risk difference"][1]

change_table <- expand.grid(change = changes$change, contrast = contrasts$contrast,
                            stringsAsFactors = FALSE)
change_results <- data.frame()
for (k in seq_len(nrow(change_table))) {
  row_change <- match(change_table$change[k], changes$change)
  earlier <- changes$earlier[row_change]
  later   <- changes$later[row_change]
  per_replicate <- rd_values[, later, , change_table$contrast[k]] -
    rd_values[, earlier, , change_table$contrast[k]]        # datasets x replicates

  estimates_i <- per_replicate[, 1]                         # the change in each dataset
  variances_i <- rep(NA, n_imputations)                     # its variance over replicates
  for (i in seq_len(n_imputations)) variances_i[i] <- var(per_replicate[i, -1])

  q_bar   <- mean(estimates_i)
  within  <- mean(variances_i)
  between <- var(estimates_i)
  total   <- within + (1 + 1 / n_imputations) * between
  rvi     <- (1 + 1 / n_imputations) * between / within
  lambda  <- (1 + 1 / n_imputations) * between / total
  df_observed <- (df_complete + 1) / (df_complete + 3) * df_complete * (1 - lambda)
  df <- df_observed
  if (lambda > 0) {
    df_old <- (n_imputations - 1) / lambda^2
    df     <- df_old * df_observed / (df_old + df_observed)
  }
  fmi <- (rvi + 2 / (df + 3)) / (rvi + 1)
  se  <- sqrt(total)

  change_results <- rbind(change_results, data.frame(
    change = change_table$change[k], contrast = change_table$contrast[k],
    earlier_model = earlier, later_model = later,
    change_estimate = q_bar,
    change_low  = q_bar - qt(0.975, df) * se,
    change_high = q_bar + qt(0.975, df) * se,
    df = df, fmi = fmi))
}


###############################################################################
# Block 8. The prediction side (File 18)
#
# For each landmark after the first: the gain in PR-AUC over the landmark
# before (elastic net, primary, weighted). For each unit File 18 shuffled
# for this script: the fall in PR-AUC when it is shuffled (importance). Each
# with its paired 95% interval. Each must be found exactly once.
###############################################################################

comparisons_18 <- file18$comparisons
bridge_units   <- file18$bridge_units

prediction <- rbind(
  data.frame(landmark = landmark_order[-1], quantity = "PR-AUC gain", unit = ""),
  data.frame(landmark = bridge_units$landmark, quantity = "Importance", unit = bridge_units$unit)
)
prediction$estimate <- NA
prediction$ci_lower <- NA
prediction$ci_upper <- NA

for (k in seq_len(nrow(prediction))) {
  if (prediction$quantity[k] == "PR-AUC gain") {
    rows <- comparisons_18$comparison == "Landmark minus the landmark before" &
      comparisons_18$analysis == "Primary" & comparisons_18$weighting == "Weighted" &
      comparisons_18$first == paste(prediction$landmark[k], "Elastic net") &
      comparisons_18$measure == "pr_auc"
  } else {
    rows <- comparisons_18$comparison == "Elastic net minus elastic net with a unit shuffled" &
      comparisons_18$landmark == prediction$landmark[k] &
      comparisons_18$second == paste0("Elastic net, ", prediction$unit[k], " shuffled") &
      comparisons_18$measure == "pr_auc"
  }
  if (sum(rows) != 1) {
    stop("File 18 has no single result for ", prediction$quantity[k], " at ",
         prediction$landmark[k], " ", prediction$unit[k], ".")
  }
  prediction$estimate[k] <- comparisons_18$difference[rows]
  prediction$ci_lower[k] <- comparisons_18$ci_lower[rows]
  prediction$ci_upper[k] <- comparisons_18$ci_upper[rows]
}


###############################################################################
# Block 9. Result 1: the science signal across the stages
#
# One row per stage x science contrast. Disclosure: a risk difference is
# shown only if its model has at least 10 events and 10 non-events in every
# completed dataset (File 12's rule).
###############################################################################

# Rows of NAs, used where a stage has no model or no landmark.
no_prediction <- prediction[1, ]
no_prediction[1, ] <- NA
no_risk_difference <- risk_differences_12[1, ]
no_risk_difference[1, ] <- NA

stage_contrast <- expand.grid(stage = stages$stage, contrast = contrasts$contrast,
                              stringsAsFactors = FALSE)
across_stages <- data.frame()

for (k in seq_len(nrow(stage_contrast))) {
  row_stage    <- match(stage_contrast$stage[k], stages$stage)
  row_contrast <- match(stage_contrast$contrast[k], contrasts$contrast)
  model    <- stages$model[row_stage]
  landmark <- stages$landmark[row_stage]

  # Association (File 12): primary analysis, and the observed-outcome sample.
  primary_rd <- risk_differences_12[risk_differences_12$analysis == primary_analysis &
                                      risk_differences_12$model %in% model &
                                      risk_differences_12$term == stage_contrast$contrast[k], ]
  observed_rd <- risk_differences_12[risk_differences_12$analysis == observed_analysis &
                                       risk_differences_12$model %in% model &
                                       risk_differences_12$term == stage_contrast$contrast[k], ]
  if (nrow(primary_rd) == 0)  primary_rd  <- no_risk_difference
  if (nrow(observed_rd) == 0) observed_rd <- no_risk_difference

  # Prediction (File 18): the same item, and the whole orientation block.
  item  <- prediction[prediction$quantity == "Importance" & prediction$landmark %in% landmark &
                        prediction$unit == contrasts$unit[row_contrast], ]
  block <- prediction[prediction$quantity == "Importance" & prediction$landmark %in% landmark &
                        prediction$unit == "Block: Science and mathematics orientation", ]
  if (nrow(item) == 0)  item  <- no_prediction
  if (nrow(block) == 0) block <- no_prediction

  # The rules for Result 1.
  association <- NA
  if (!is.na(primary_rd$estimate)) association <- "Not clear"
  if (!is.na(primary_rd$estimate) && (primary_rd$low > 0 || primary_rd$high < 0)) association <- "Present"
  predictive_information <- NA
  if (!is.na(item$estimate)) predictive_information <- "Not clear"
  if (!is.na(item$estimate) && item$ci_lower > 0) predictive_information <- "Present"

  # Disclosure: File 12's model rule, in every completed dataset.
  shown <- is.na(primary_rd$estimate) ||
    (primary_rd$smallest_events >= disclosure_threshold &&
       primary_rd$smallest_non_events >= disclosure_threshold)
  observed_estimate <- observed_rd$estimate
  if (!is.na(observed_estimate) && (observed_rd$smallest_events < disclosure_threshold ||
                                    observed_rd$smallest_non_events < disclosure_threshold)) {
    observed_estimate <- NA
  }

  across_stages <- rbind(across_stages, data.frame(
    stage = stage_contrast$stage[k], model = model, landmark = landmark,
    contrast = stage_contrast$contrast[k],
    risk_difference = primary_rd$estimate, rd_low = primary_rd$low, rd_high = primary_rd$high,
    rd_observed_outcome_sample = observed_estimate,
    association = association,
    item_importance = item$estimate, item_importance_low = item$ci_lower,
    item_importance_high = item$ci_upper,
    block_importance = block$estimate, block_importance_low = block$ci_lower,
    block_importance_high = block$ci_upper,
    predictive_information = predictive_information,
    shown = shown))
}


###############################################################################
# Block 10. Result 2: the change at each stage and its predictive value
#
# One row per change x science contrast, plus the prediction-only stage
# (L25b). The rules are in the header.
###############################################################################

stage_changes <- data.frame()
for (k in seq_len(nrow(change_results))) {
  row_change <- match(change_results$change[k], changes$change)
  landmark   <- changes$landmark[row_change]

  # Association result.
  earlier_rd <- risk_differences_12$estimate[risk_differences_12$analysis == primary_analysis &
                                               risk_differences_12$model == change_results$earlier_model[k] &
                                               risk_differences_12$term == change_results$contrast[k]]
  association_result <- "No clear change"
  if (change_results$change_low[k] > 0 | change_results$change_high[k] < 0) {
    association_result <- "Strengthened"
    if (sign(change_results$change_estimate[k]) == -sign(earlier_rd)) association_result <- "Weakened"
  }

  # Prediction: the landmark gain, the importance of the block the landmark
  # adds, and (L20, L25a) of the matching single feature.
  gain  <- prediction[prediction$quantity == "PR-AUC gain" & prediction$landmark %in% landmark, ]
  added <- prediction[prediction$quantity == "Importance" & prediction$landmark %in% landmark &
                        startsWith(prediction$unit, "Block:") &
                        prediction$unit != "Block: Science and mathematics orientation", ]
  matched <- prediction[prediction$quantity == "Importance" & prediction$landmark %in% landmark &
                          prediction$unit %in% c("Feature: he_pathway20", "Feature: degree_pathway25"), ]
  if (nrow(gain) == 0)    gain    <- no_prediction
  if (nrow(added) == 0)   added   <- no_prediction
  if (nrow(matched) == 0) matched <- no_prediction

  # The rules for Result 2 (only where the stage has a landmark).
  prediction_result <- NA
  agreement <- NA
  if (!is.na(gain$estimate)) {
    prediction_result <- "No clear added value"
    if (gain$ci_lower > 0) prediction_result <- "Adds predictive value"
    if (gain$ci_upper < 0) prediction_result <- "Lowers performance"
    agreement <- "Differ"
    association_shows <- association_result != "No clear change"
    prediction_shows  <- prediction_result != "No clear added value"
    if (association_shows == prediction_shows) agreement <- "Agree"
  }

  # Observed-outcome sample: the change as a difference of pooled estimates.
  observed_earlier <- risk_differences_12$estimate[risk_differences_12$analysis == observed_analysis &
                                                     risk_differences_12$model == change_results$earlier_model[k] &
                                                     risk_differences_12$term == change_results$contrast[k]]
  observed_later   <- risk_differences_12$estimate[risk_differences_12$analysis == observed_analysis &
                                                     risk_differences_12$model == change_results$later_model[k] &
                                                     risk_differences_12$term == change_results$contrast[k]]

  stage_changes <- rbind(stage_changes, data.frame(
    change = change_results$change[k], contrast = change_results$contrast[k],
    change_in_risk_difference = change_results$change_estimate[k],
    change_low = change_results$change_low[k], change_high = change_results$change_high[k],
    change_fmi = change_results$fmi[k],
    change_observed_outcome_sample = observed_later - observed_earlier,
    association_result = association_result,
    landmark = landmark,
    pr_auc_gain = gain$estimate, pr_auc_gain_low = gain$ci_lower, pr_auc_gain_high = gain$ci_upper,
    added_block_importance = added$estimate, added_block_importance_low = added$ci_lower,
    added_block_importance_high = added$ci_upper,
    matched_feature_importance = matched$estimate,
    matched_feature_importance_low = matched$ci_lower,
    matched_feature_importance_high = matched$ci_upper,
    prediction_result = prediction_result,
    agreement = agreement))
}

# The prediction-only stage (SRET job at 25, L25b).
gain  <- prediction[prediction$quantity == "PR-AUC gain" & prediction$landmark == "L25b", ]
added <- prediction[prediction$quantity == "Importance" & prediction$landmark == "L25b" &
                      prediction$unit == "Block: Occupation at 25 (sensitivity)", ]
prediction_result <- "No clear added value"
if (gain$ci_lower > 0) prediction_result <- "Adds predictive value"
if (gain$ci_upper < 0) prediction_result <- "Lowers performance"
stage_changes <- rbind(stage_changes, data.frame(
  change = "SRET job at 25 (prediction only)", contrast = NA,
  change_in_risk_difference = NA, change_low = NA, change_high = NA, change_fmi = NA,
  change_observed_outcome_sample = NA, association_result = NA, landmark = "L25b",
  pr_auc_gain = gain$estimate, pr_auc_gain_low = gain$ci_lower, pr_auc_gain_high = gain$ci_upper,
  added_block_importance = added$estimate, added_block_importance_low = added$ci_lower,
  added_block_importance_high = added$ci_upper,
  matched_feature_importance = NA, matched_feature_importance_low = NA,
  matched_feature_importance_high = NA,
  prediction_result = prediction_result, agreement = NA))


###############################################################################
# Block 11. Disclosure check for the changes
#
# A change is shown only if both of its models pass File 12's model rule.
###############################################################################

model_check <- expand.grid(analysis = c(primary_analysis, observed_analysis), model = models$model,
                           stringsAsFactors = FALSE)
model_check$passes <- NA
for (k in seq_len(nrow(model_check))) {
  rows <- risk_differences_12[risk_differences_12$analysis == model_check$analysis[k] &
                                risk_differences_12$model == model_check$model[k], ]
  model_check$passes[k] <- all(rows$smallest_events >= disclosure_threshold) &&
    all(rows$smallest_non_events >= disclosure_threshold)
}
primary_passes  <- model_check$passes[model_check$analysis == primary_analysis]
observed_passes <- model_check$passes[model_check$analysis == observed_analysis]
names(primary_passes)  <- models$model
names(observed_passes) <- models$model

change_shown   <- rep(TRUE, nrow(stage_changes))
observed_shown <- rep(TRUE, nrow(stage_changes))
for (k in seq_len(nrow(change_results))) {
  change_shown[k] <- primary_passes[change_results$earlier_model[k]] &&
    primary_passes[change_results$later_model[k]]
  observed_shown[k] <- observed_passes[change_results$earlier_model[k]] &&
    observed_passes[change_results$later_model[k]]
}
association_columns <- c("change_in_risk_difference", "change_low", "change_high", "change_fmi",
                         "change_observed_outcome_sample", "association_result", "agreement")
stage_changes[!change_shown, association_columns] <- NA
stage_changes$change_observed_outcome_sample[!observed_shown] <- NA
across_stages[!across_stages$shown,
              c("risk_difference", "rd_low", "rd_high", "rd_observed_outcome_sample", "association")] <- NA


###############################################################################
# Block 12. Plots
#
# One row per science contrast. Left: the risk difference under M0-M4
# (95% CI). Right: the importance of the same item (filled) and of the
# orientation block (open) at each landmark (95% CI). Both panels use the
# same six stages, left to right; a dash marks a stage with no model (L25b)
# or no landmark of its own (M1).
###############################################################################

x <- seq_len(nrow(stages))
model_labels    <- c("M0", "M1", "M2", "M3", "M4", "-")
landmark_labels <- c("L14", "-", "L16", "L20", "L25a", "L25b")

pdf(file.path(review_dir, "20_bridge_plots.pdf"), width = 11, height = 11)
par(mfrow = c(nrow(contrasts), 2), mar = c(4.5, 4.5, 3, 1))
for (j in seq_len(nrow(contrasts))) {
  rows <- across_stages[across_stages$contrast == contrasts$contrast[j], ]
  rows <- rows[match(stages$stage, rows$stage), ]      # in stage order

  y_range <- range(c(0, rows$rd_low, rows$rd_high), na.rm = TRUE)
  plot(x, rows$risk_difference, pch = 19, ylim = y_range, xlim = c(0.5, length(x) + 0.5),
       xaxt = "n", xlab = "Stage", ylab = "Risk difference (imputation analysis)",
       main = paste("Association:", contrasts$contrast[j]), cex.main = 0.9)
  axis(1, at = x, labels = model_labels)
  segments(x, rows$rd_low, x, rows$rd_high)
  abline(h = 0, lty = 2, col = "grey50")

  y_range <- range(c(0, rows$item_importance_low, rows$item_importance_high,
                     rows$block_importance_low, rows$block_importance_high), na.rm = TRUE)
  plot(x - 0.1, rows$item_importance, pch = 19, ylim = y_range, xlim = c(0.5, length(x) + 0.5),
       xaxt = "n", xlab = "Stage", ylab = "Fall in PR-AUC when shuffled",
       main = "Predictive information: item (filled) and orientation block (open)", cex.main = 0.9)
  axis(1, at = x, labels = landmark_labels)
  segments(x - 0.1, rows$item_importance_low, x - 0.1, rows$item_importance_high)
  points(x + 0.1, rows$block_importance, pch = 1)
  segments(x + 0.1, rows$block_importance_low, x + 0.1, rows$block_importance_high, lty = 3)
  abline(h = 0, lty = 2, col = "grey50")
}
invisible(dev.off())


###############################################################################
# Block 13. Save and write the review tables
#
# Disclosure check (before writing): pooled estimates for whole-population
# models that pass File 12's rule (Block 11), and whole-cohort prediction
# measures from File 18. No counts of people.
###############################################################################

saveRDS(list(rd_values = rd_values, change_results = change_results, prediction = prediction,
             across_stages = across_stages, stage_changes = stage_changes,
             se_ratio = cbind(check_table, se_ratio = se_ratio),
             settings = list(bootstrap_replicates = bootstrap_replicates,
                             bootstrap_seed = bootstrap_seed)),
        file.path(data_dir, "20_bridge.rds"))

numeric_columns <- c("risk_difference", "rd_low", "rd_high", "rd_observed_outcome_sample",
                     "item_importance", "item_importance_low", "item_importance_high",
                     "block_importance", "block_importance_low", "block_importance_high")
across_review <- across_stages[, setdiff(names(across_stages), "shown")]
across_review[, numeric_columns] <- round(across_review[, numeric_columns], 5)
across_review$disclosure_checked <- TRUE
write.csv(across_review, file.path(review_dir, "20_orientation_across_stages.csv"),
          row.names = FALSE, na = "")

numeric_columns <- c("change_in_risk_difference", "change_low", "change_high", "change_fmi",
                     "change_observed_outcome_sample", "pr_auc_gain", "pr_auc_gain_low",
                     "pr_auc_gain_high", "added_block_importance", "added_block_importance_low",
                     "added_block_importance_high", "matched_feature_importance",
                     "matched_feature_importance_low", "matched_feature_importance_high")
changes_review <- stage_changes
changes_review[, numeric_columns] <- round(changes_review[, numeric_columns], 5)
changes_review$disclosure_checked <- TRUE
write.csv(changes_review, file.path(review_dir, "20_stage_changes.csv"), row.names = FALSE, na = "")

disclosure_log <- data.frame(
  table = c("20_orientation_across_stages.csv", "20_stage_changes.csv", "20_bridge_plots.pdf"),
  rule = c("Risk differences only for models with at least 10 events and 10 non-events in every completed dataset; whole-cohort prediction measures",
           "Changes only when both models pass the model rule; whole-cohort prediction measures",
           "Drawn from 20_orientation_across_stages.csv"),
  rows_suppressed = c(sum(!across_stages$shown), sum(!change_shown), 0)
)
disclosure_log$disclosure_checked <- TRUE
write.csv(disclosure_log, file.path(review_dir, "20_disclosure_log.csv"), row.names = FALSE)


###############################################################################
# Block 14. Summary
#
# Checks and run time only. The results are in the review tables.
###############################################################################

cat("\nFile 20 finished.\n")
cat("  Refits reproduce File 12's risk differences:   yes (", nrow(check_table), "checked )\n")
cat("  Bootstrap run time (minutes):                 ", round(bootstrap_minutes, 1), "\n")
cat("  Fits that didn't converge:                    ", fits_not_converged, "of",
    nrow(fit_table) * (bootstrap_replicates + 1), "\n")
cat("  Bootstrap SE / delta-method SE (median ratio):", round(median_se_ratio, 2),
    "(close to 1 is expected)\n")
cat("  Changes estimated:                            ", nrow(change_results), "\n")
cat("  Review: 20_orientation_across_stages.csv, 20_stage_changes.csv, 20_bridge_plots.pdf,\n")
cat("    20_disclosure_log.csv\n")
