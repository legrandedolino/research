###############################################################################
# 12. Pooled survey analysis
#
# Purpose
#   Fit the analysis models of File 5 to every completed dataset from the
#   multiple imputation, with the survey design, and combine the results with
#   Rubin's rules (Rubin 1987). These are the study's population estimates of
#   how adolescent science orientation relates to professional SRET
#   employment at 32, allowing for attrition under missing at random.
#   Missing not at random is examined in File 13.
#
# Analyses (CLS user guide: impute using everyone, then choose the analysis
# sample, and examine how that choice affects the results)
#   Primary: the main imputations (File 8) in the whole target population,
#     weighted by the initial design weight, as Silverwood et al. (2024)
#     weighted their imputation estimates.
#   Sensitivity, analysis sample:
#     - Wave-9 respondents (with a Wave-9 weight), weighted by the Wave-9
#       final weight (W9FINWT);
#     - the observed-outcome sample (the prediction cohort of File 3),
#       weighted by W9FINWT. Only the predictors are imputed here, so this
#       is the closest imputation counterpart of the complete-case analysis.
#   Supplementary imputation (File 9), without the CLS predictors of
#     non-response, in the primary sample: what those predictors add
#     (Silverwood et al. 2024).
#
# Models (prespecified in File 5)
#   M0 (baseline) and the sequential models M1-M4, with the overall SRET
#   outcome, and M0 employed-only (SRET among people in paid work, with
#   employment itself imputed where missing). For each model: odds ratios
#   and the standardized risk differences for the three science contrasts
#   (File 5). For each analysis: weighted SRET prevalence overall and among
#   the employed. The associations are descriptive, not causal effects.
#
# Survey design
#   Each completed dataset is analysed with svydesign(): PSU, the original
#   sampling strata (not the pooled strata used as an imputation predictor)
#   and the weight of the analysis. The Wave-9 samples are analysed as
#   subsets of the Wave-9 respondent design, so their standard errors reflect
#   its full sample design.
#
# Rubin's rules (for each coefficient, risk difference and prevalence)
#   For m completed datasets with estimates q_i and variances u_i:
#     pooled estimate    q = mean of q_i
#     within variance    W = mean of u_i
#     between variance   B = variance of q_i
#     total variance     T = W + (1 + 1/m) B
#     relative increase in variance due to missing data   RVI = (1 + 1/m) B / W
#     proportion of variance due to missing data          lambda = (1 + 1/m) B / T
#     degrees of freedom: Barnard and Rubin (1999), with the survey design's
#       degrees of freedom (PSUs minus strata) as the complete-data value
#     fraction of missing information   FMI = (RVI + 2 / (df + 3)) / (RVI + 1)
#     relative efficiency               RE = 1 / (1 + FMI / m)
#     Monte Carlo error of the estimate  sqrt(B / m); White, Royston and Wood
#       (2011) suggest it should be below 10% of the standard error
#   95% confidence intervals and p-values use the t distribution with df.
#
# Number of imputations (CLS user guide)
#   At least 50, and at least the largest FMI as a percentage. File 7
#   estimated the FMI from a 5-imputation test run; it is checked here on
#   the main imputations. If the rule isn't met for the primary analysis,
#   more imputations should be added (File 8a) and this script rerun.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   Counts of people vary between completed datasets, so every rule uses the
#   smallest count across them:
#     - a model's results, and a prevalence, are shown only if every
#       completed dataset has at least 10 events and 10 non-events;
#     - a coefficient for a category (or any 0/1 term) is shown only if, in
#       every completed dataset, at least 10 people have it, with at least 10
#       events and 10 non-events among them;
#     - sample sizes are shown only where they can't reveal a small group by
#       subtraction (Block 8);
#     - complete-case results in the comparison table follow File 5's rules.
#
# Reads
#   data_dir/05_complete_case_results.rds
#   data_dir/07_imputation_spec.rds         (design table: PSU, strata, weights)
#   data_dir/08_imputations_main.rds
#   data_dir/09_imputations_supplementary.rds
#
# Writes
#   data_dir/12_pooled_estimates.rds        exact pooled results (for Files 13
#                                           and 20)
#   output_dir/private_12_*.csv             exact tables (never shared)
#   review_dir/12_prevalence.csv, 12_coefficients.csv,
#     12_risk_differences.csv, 12_imputation_number_check.csv,
#     12_comparison.csv, 12_disclosure_log.csv
#
# References
#   Barnard J, Rubin DB (1999). Small-sample degrees of freedom with multiple
#   imputation. Biometrika 86(4): 948-955.
#   Rubin DB (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
#   Silverwood RJ, et al. (2024). Longitudinal and Life Course Studies 15(2):
#   227-250.
#   White IR, Royston P, Wood AM (2011). Multiple imputation using chained
#   equations: issues and guidance for practice. Statistics in Medicine
#   30(4): 377-399.
###############################################################################

suppressPackageStartupMessages(library(mice))
suppressPackageStartupMessages(library(survey))


###############################################################################
# Block 1. Paths and settings
#
# Folder locations are set once, in File 0, and read here from
# project_paths.rds. Run this script with the working directory set to the
# folder that holds the scripts.
###############################################################################

paths      <- readRDS("project_paths.rds")
output_dir <- paths$output_dir   # exact (private) tables are written here
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs

disclosure_threshold <- 10   # smallest count that may be shown
minimum_m            <- 50   # CLS user guide: at least 50 imputations

# Standard errors compare PSUs within each stratum. A stratum with only one
# PSU in the data a design is built from has nothing to compare with, and
# the survey package stops by default. "adjust" measures that PSU's
# deviation from the overall mean instead, which is slightly conservative
# (as in File 5).
options(survey.lonely.psu = "adjust")


###############################################################################
# Block 2. Read the models, the design table and the imputations
###############################################################################

file5           <- readRDS(file.path(data_dir, "05_complete_case_results.rds"))
imputation_spec <- readRDS(file.path(data_dir, "07_imputation_spec.rds"))
design_data     <- imputation_spec$design_data     # same row order as the imputations

imputations <- list(
  main          = readRDS(file.path(data_dir, "08_imputations_main.rds")),
  supplementary = readRDS(file.path(data_dir, "09_imputations_supplementary.rds"))
)

# The models: File 5's, except its extra complete-case version of M0
# (fitted to M0's own complete cases), which has no imputation counterpart.
models    <- file5$models[file5$models$model != "M0 baseline (own complete cases)", ]
contrasts <- file5$contrasts


###############################################################################
# Block 3. The analyses
#
# One row per analysis: which imputations, which sample, which weight.
###############################################################################

analyses <- data.frame(
  analysis = c("Primary: target population, design weight",
               "Sensitivity: Wave-9 respondents, Wave-9 weight",
               "Sensitivity: observed-outcome sample, Wave-9 weight",
               "Supplementary imputation: target population, design weight"),
  imputation_set = c("main", "main", "main", "supplementary"),
  sample = c("target_population", "w9_respondents", "observed_outcome",
             "target_population")
)


###############################################################################
# Block 4. Fit every model to every completed dataset
#
# One row of fit_table per imputation set x completed dataset (i) x
# analysis, in that order. For each row:
#   a. the completed dataset, with the design table joined and the variables
#      prepared as in File 5: HE likelihood as unordered categories; SRET
#      among the employed as 0/1;
#   b. the analysis's survey design:
#        target population  all rows, initial design weight
#        Wave-9 respondents respondents with a Wave-9 weight, W9FINWT
#        observed outcome   the respondents with an observed outcome, as a
#                           subset of the respondent design;
#   c. SRET prevalence, overall and among the employed (svymean);
#   d. each model: a survey-weighted logistic regression (quasibinomial, as
#      in File 5); the employed-only model is fitted to the subset in paid
#      work (imputed where missing); and the standardized risk differences of
#      the science contrasts, with delta-method standard errors, as in File 5.
# Each result is stored with its estimate, variance, the design's degrees of
# freedom and, for the disclosure checks, the numbers of people, events and
# non-events. Results are collected in a list and joined at the end.
###############################################################################

fit_table <- data.frame()
for (set in names(imputations)) {
  set_rows <- expand.grid(analysis = analyses$analysis[analyses$imputation_set == set],
                          imputation = seq_len(imputations[[set]]$m), stringsAsFactors = FALSE)
  fit_table <- rbind(fit_table, data.frame(set = set, imputation = set_rows$imputation,
                                           analysis = set_rows$analysis))
}

pieces         <- list()   # one data frame of results per fit
not_converged  <- 0

for (k in seq_len(nrow(fit_table))) {
  set           <- fit_table$set[k]
  i             <- fit_table$imputation[k]
  analysis_name <- fit_table$analysis[k]
  sample_name   <- analyses$sample[analyses$analysis == analysis_name]
  if (i == 1 && analysis_name == analyses$analysis[analyses$imputation_set == set][1]) {
    cat("\nFitting the models to the", imputations[[set]]$m, set, "imputations\n")
  }

  # a. The completed dataset.
  analysis_data <- cbind(design_data, complete(imputations[[set]], i))
  analysis_data$he_likelihood16 <- factor(analysis_data$he_likelihood16, ordered = FALSE)
  analysis_data$sret_among_employed32 <-
    as.numeric(analysis_data$sret_among_employed32 == "SRET professional")

  # b. The analysis's survey design.
  respondents <- analysis_data[analysis_data$responded_w9 %in% TRUE &
                                 !is.na(analysis_data$w9_final_weight), ]
  if (sample_name == "target_population") {
    design <- svydesign(ids = ~psu, strata = ~stratum, weights = ~design_weight, nest = TRUE,
                        data = analysis_data)
  }
  if (sample_name == "w9_respondents") {
    design <- svydesign(ids = ~psu, strata = ~stratum, weights = ~w9_final_weight, nest = TRUE,
                        data = respondents)
  }
  if (sample_name == "observed_outcome") {
    design_respondents <- svydesign(ids = ~psu, strata = ~stratum, weights = ~w9_final_weight,
                                    nest = TRUE, data = respondents)
    design <- subset(design_respondents, outcome_observed32)
  }
  design_employed <- subset(design, employed32 == "In paid work")

  # c. Prevalence, overall and among the employed.
  for (outcome_label in c("SRET overall", "SRET among the employed")) {
    if (outcome_label == "SRET overall") {
      prevalence_design <- design
      y <- design$variables$sret_professional_32
      prevalence <- svymean(~sret_professional_32, design)
    } else {
      prevalence_design <- design_employed
      y <- design_employed$variables$sret_among_employed32
      prevalence <- svymean(~sret_among_employed32, design_employed)
    }
    pieces[[length(pieces) + 1]] <- data.frame(
      analysis = analysis_name, model = outcome_label,
      quantity = "prevalence", term = outcome_label, imputation = i,
      estimate = coef(prevalence)[[1]], variance = SE(prevalence)[[1]]^2,
      df_complete = degf(prevalence_design),
      n = length(y), events = sum(y == 1), non_events = sum(y == 0),
      term_is_binary = FALSE, people_with_term = NA,
      events_with_term = NA, non_events_with_term = NA
    )
  }

  # d. The models.
  for (m in seq_len(nrow(models))) {
    model_outcome <- models$outcome[m]
    model_design  <- design
    if (model_outcome == "sret_among_employed32") model_design <- design_employed
    model_formula <- as.formula(paste(model_outcome, "~", models$predictors[m]))

    fit <- svyglm(model_formula, design = model_design, family = quasibinomial(),
                  control = glm.control(maxit = 100))
    if (!fit$converged) not_converged <- not_converged + 1
    if (any(is.na(coef(fit)))) {
      stop(models$model[m], " (", analysis_name, ", imputation ", i, "): some ",
           "terms can't be estimated: ",
           paste(names(coef(fit))[is.na(coef(fit))], collapse = ", "))
    }

    y <- model_design$variables[[model_outcome]]
    w <- weights(model_design)
    X <- model.matrix(fit)
    df_complete <- degf(model_design)

    # Coefficients (log odds), with the counts for the disclosure check.
    pieces[[length(pieces) + 1]] <- data.frame(
      analysis = analysis_name, model = models$model[m],
      quantity = "coefficient", term = names(coef(fit)), imputation = i,
      estimate = unname(coef(fit)), variance = unname(diag(vcov(fit))),
      df_complete = df_complete,
      n = length(y), events = sum(y == 1), non_events = sum(y == 0),
      term_is_binary       = unname(colSums(X != 0 & X != 1) == 0),
      people_with_term     = unname(colSums(X != 0)),
      events_with_term     = unname(colSums((X != 0) * (y == 1))),
      non_events_with_term = unname(colSums((X != 0) * (y == 0)))
    )

    # Standardized risk differences for the science contrasts: each person's
    # predicted probability with the item set to "high" and then to "low",
    # and the weighted averages compared.
    for (j in seq_len(nrow(contrasts))) {
      variable  <- contrasts$variable[j]
      data_high <- model_design$variables
      data_low  <- model_design$variables
      numeric_item <- is.numeric(data_high[[variable]])
      if (numeric_item)  data_high[[variable]] <- as.numeric(contrasts$high[j])
      if (numeric_item)  data_low[[variable]]  <- as.numeric(contrasts$low[j])
      if (!numeric_item) data_high[[variable]] <- factor(contrasts$high[j], levels = levels(data_high[[variable]]))
      if (!numeric_item) data_low[[variable]]  <- factor(contrasts$low[j],  levels = levels(data_low[[variable]]))
      predictor_terms <- delete.response(terms(fit))
      X_high <- model.matrix(predictor_terms, data = data_high, xlev = fit$xlevels)
      X_low  <- model.matrix(predictor_terms, data = data_low,  xlev = fit$xlevels)
      p_high <- as.vector(plogis(X_high %*% coef(fit)))
      p_low  <- as.vector(plogis(X_low  %*% coef(fit)))

      rd <- sum(w * (p_high - p_low)) / sum(w)
      gradient <- colSums(w * (p_high * (1 - p_high) * X_high -
                               p_low  * (1 - p_low)  * X_low)) / sum(w)
      rd_variance <- as.numeric(t(gradient) %*% vcov(fit) %*% gradient)

      pieces[[length(pieces) + 1]] <- data.frame(
        analysis = analysis_name, model = models$model[m],
        quantity = "risk difference", term = contrasts$contrast[j],
        imputation = i, estimate = rd, variance = rd_variance,
        df_complete = df_complete,
        n = length(y), events = sum(y == 1), non_events = sum(y == 0),
        term_is_binary = FALSE, people_with_term = NA,
        events_with_term = NA, non_events_with_term = NA
      )
    }
  }

  last_analysis <- tail(analyses$analysis[analyses$imputation_set == set], 1)
  if (i %% 10 == 0 && analysis_name == last_analysis) cat("  ", i, "completed datasets done\n")
}

# Join all the pieces into one table, in the order they were made
# (do.call passes every piece in the list to rbind at once).
results_by_imputation <- do.call(rbind, pieces)


###############################################################################
# Block 5. Combine with Rubin's rules
#
# One pooled row for each analysis x model x quantity x term. The smallest
# counts across completed datasets are kept for the disclosure checks.
###############################################################################

# One key per estimate (analysis | model | quantity | term); each key's rows
# are its values in the completed datasets.
key <- paste(results_by_imputation$analysis, results_by_imputation$model,
             results_by_imputation$quantity, results_by_imputation$term, sep = " | ")

pooled <- data.frame()
for (g in unique(key)) {
  rows <- results_by_imputation[key == g, ]
  m    <- nrow(rows)

  q_bar  <- mean(rows$estimate)
  within <- mean(rows$variance)
  between <- var(rows$estimate)
  total  <- within + (1 + 1 / m) * between
  rvi    <- (1 + 1 / m) * between / within
  lambda <- (1 + 1 / m) * between / total

  # Barnard-Rubin degrees of freedom.
  df_complete <- rows$df_complete[1]
  df_observed <- (df_complete + 1) / (df_complete + 3) * df_complete * (1 - lambda)
  if (lambda > 0) {
    df_old <- (m - 1) / lambda^2
    df     <- df_old * df_observed / (df_old + df_observed)
  } else {
    df <- df_observed
  }
  fmi <- (rvi + 2 / (df + 3)) / (rvi + 1)

  se       <- sqrt(total)
  critical <- qt(0.975, df = df)

  pooled <- rbind(pooled, data.frame(
    analysis = rows$analysis[1], model = rows$model[1],
    quantity = rows$quantity[1], term = rows$term[1],
    imputations = m,
    estimate = q_bar, se = se,
    low = q_bar - critical * se, high = q_bar + critical * se,
    p_value = 2 * pt(-abs(q_bar / se), df = df),
    df = df, rvi = rvi, lambda = lambda, fmi = fmi,
    relative_efficiency = 1 / (1 + fmi / m),
    mc_error = sqrt(between / m),
    mc_error_vs_se = sqrt(between / m) / se,
    # For the disclosure checks (not shown in the checked tables).
    n_mean = mean(rows$n),
    smallest_events = min(rows$events), smallest_non_events = min(rows$non_events),
    term_is_binary = rows$term_is_binary[1],
    smallest_people_with_term = min(rows$people_with_term),
    smallest_events_with_term = min(rows$events_with_term),
    smallest_non_events_with_term = min(rows$non_events_with_term)
  ))
}


###############################################################################
# Block 6. Result tables
#
#   prevalence         weighted SRET %, overall and among the employed
#   coefficients       log odds, odds ratios with 95% CIs, p-values, and
#                      the missing-information statistics
#   risk differences   the science contrasts
#   imputation number  the CLS rule: at least 50 imputations and at least
#                      the largest FMI (%) of each model's coefficients
###############################################################################

prevalence <- pooled[pooled$quantity == "prevalence", ]
prevalence$weighted_percent      <- round(100 * prevalence$estimate, 2)
prevalence$weighted_percent_low  <- round(100 * prevalence$low, 2)
prevalence$weighted_percent_high <- round(100 * prevalence$high, 2)

coefficients <- pooled[pooled$quantity == "coefficient", ]
coefficients$odds_ratio      <- exp(coefficients$estimate)
coefficients$odds_ratio_low  <- exp(coefficients$low)
coefficients$odds_ratio_high <- exp(coefficients$high)

risk_differences <- pooled[pooled$quantity == "risk difference", ]

imputation_check <- data.frame()
for (a in unique(coefficients$analysis)) {
  for (model_name in models$model) {
    rows <- coefficients[coefficients$analysis == a & coefficients$model == model_name, ]
    largest <- max(rows$fmi)
    imputation_check <- rbind(imputation_check, data.frame(
      analysis = a, model = model_name,
      largest_fmi = largest,
      term_with_largest_fmi = rows$term[which.max(rows$fmi)],
      imputations = rows$imputations[1],
      imputations_needed = max(minimum_m, ceiling(100 * largest)),
      smallest_events = min(rows$smallest_events),
      smallest_non_events = min(rows$smallest_non_events)
    ))
  }
}
imputation_check$rule_met <- imputation_check$imputations >= imputation_check$imputations_needed

write.csv(prevalence, file.path(output_dir, "private_12_prevalence.csv"), row.names = FALSE)
write.csv(coefficients, file.path(output_dir, "private_12_coefficients.csv"), row.names = FALSE)
write.csv(risk_differences, file.path(output_dir, "private_12_risk_differences.csv"),
          row.names = FALSE)
write.csv(imputation_check, file.path(output_dir, "private_12_imputation_number_check.csv"),
          row.names = FALSE)


###############################################################################
# Block 7. Disclosure checks and review tables
###############################################################################

disclosure_log <- data.frame(table = character(), primary_suppressed = integer(),
                             secondary_suppressed = integer())

# Which analysis x model combinations fail the model rule: fewer than 10
# events or non-events in some completed dataset.
model_rows <- pooled[pooled$quantity == "coefficient", ]
small_model_key <- unique(paste(model_rows$analysis, model_rows$model)[
  model_rows$smallest_events < disclosure_threshold |
  model_rows$smallest_non_events < disclosure_threshold])

# --- Disclosure check: prevalence ---------------------------------------------
# A prevalence is shown only if every completed dataset has at least 10
# events and 10 non-events. Sample sizes (average over completed datasets,
# which is the exact count where nothing in the sample definition is
# imputed) are shown only if at least 10, and hidden where subtraction
# would reveal a group of 1-9 people:
#   - the Wave-9 respondent sample minus the observed-outcome sample (people
#     with an unresolved outcome), both of which are exact counts;
#   - within an analysis, the whole sample minus the employed (people not in
#     paid work), where both are exact counts.
checked <- prevalence
too_small <- checked$smallest_events < disclosure_threshold |
             checked$smallest_non_events < disclosure_threshold
checked[too_small, c("weighted_percent", "weighted_percent_low",
                     "weighted_percent_high", "fmi")] <- NA
checked$n <- round(checked$n_mean)
hide_n <- checked$n < disclosure_threshold

overall <- checked$model == "SRET overall"
n_respondents <- checked$n_mean[overall & checked$analysis == analyses$analysis[2]]
n_observed    <- checked$n_mean[overall & checked$analysis == analyses$analysis[3]]
difference    <- n_respondents - n_observed
if (difference > 0 && difference < disclosure_threshold) {
  hide_n[overall & checked$analysis == analyses$analysis[2]] <- TRUE
}
for (a in unique(checked$analysis)) {
  n_all      <- checked$n_mean[checked$analysis == a & overall]
  n_employed <- checked$n_mean[checked$analysis == a & !overall]
  exact <- n_all == round(n_all) && n_employed == round(n_employed)
  if (exact && (n_all - n_employed) > 0 && (n_all - n_employed) < disclosure_threshold) {
    hide_n[checked$analysis == a & !overall] <- TRUE
  }
}
checked$n[hide_n] <- NA
checked$fmi <- round(checked$fmi, 3)
checked <- checked[, c("analysis", "model", "n", "weighted_percent",
                       "weighted_percent_low", "weighted_percent_high", "fmi",
                       "imputations")]
names(checked)[2] <- "outcome"
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "12_prevalence.csv"), row.names = FALSE, na = "")
checked_prevalence <- checked
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "12_prevalence", primary_suppressed = sum(too_small) + sum(hide_n),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------

# --- Disclosure check: coefficients -------------------------------------------
# All coefficients of a model that fails the model rule are hidden. A 0/1
# term (a category, or science as a favourite subject) is hidden unless, in
# every completed dataset, at least 10 people have it, with at least 10
# events and 10 non-events among them. The counts are removed from the
# shared table.
checked <- coefficients
hidden_model <- paste(checked$analysis, checked$model) %in% small_model_key
hidden_term  <- checked$term_is_binary &
  (checked$smallest_people_with_term < disclosure_threshold |
   checked$smallest_events_with_term < disclosure_threshold |
   checked$smallest_non_events_with_term < disclosure_threshold)
hide <- hidden_model | hidden_term
shown_columns <- c("estimate", "se", "odds_ratio", "odds_ratio_low", "odds_ratio_high",
                   "p_value", "df", "rvi", "lambda", "fmi", "relative_efficiency",
                   "mc_error", "mc_error_vs_se")
checked[hide, shown_columns] <- NA
checked <- data.frame(
  analysis = checked$analysis, model = checked$model, term = checked$term,
  log_odds = round(checked$estimate, 4), se = round(checked$se, 4),
  odds_ratio = round(checked$odds_ratio, 3),
  odds_ratio_low = round(checked$odds_ratio_low, 3),
  odds_ratio_high = round(checked$odds_ratio_high, 3),
  p_value = signif(checked$p_value, 3), df = round(checked$df, 1),
  rvi = round(checked$rvi, 3), lambda = round(checked$lambda, 3),
  fmi = round(checked$fmi, 3),
  relative_efficiency = round(checked$relative_efficiency, 4),
  mc_error = signif(checked$mc_error, 3),
  mc_error_vs_se = round(checked$mc_error_vs_se, 3),
  imputations = checked$imputations
)
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "12_coefficients.csv"), row.names = FALSE, na = "")
checked_coefficients <- checked
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "12_coefficients", primary_suppressed = sum(hide),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------

# --- Disclosure check: risk differences ---------------------------------------
# Hidden for any model that fails the model rule.
checked <- risk_differences
hide <- paste(checked$analysis, checked$model) %in% small_model_key
checked <- data.frame(
  analysis = checked$analysis, model = checked$model, contrast = checked$term,
  risk_difference = round(checked$estimate, 4), se = round(checked$se, 4),
  low = round(checked$low, 4), high = round(checked$high, 4),
  p_value = signif(checked$p_value, 3), df = round(checked$df, 1),
  fmi = round(checked$fmi, 3), mc_error = signif(checked$mc_error, 3),
  mc_error_vs_se = round(checked$mc_error_vs_se, 3),
  imputations = checked$imputations
)
checked[hide, c("risk_difference", "se", "low", "high", "p_value", "df", "fmi",
                "mc_error", "mc_error_vs_se")] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "12_risk_differences.csv"), row.names = FALSE, na = "")
checked_risk_differences <- checked
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "12_risk_differences", primary_suppressed = sum(hide),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------

# --- Disclosure check: number of imputations -----------------------------------
# A model's row is shown only if it passes the model rule; the event counts
# are removed.
checked <- imputation_check
hide <- checked$smallest_events < disclosure_threshold |
        checked$smallest_non_events < disclosure_threshold
checked[hide, c("largest_fmi", "term_with_largest_fmi", "imputations_needed",
                "rule_met")] <- NA
checked$largest_fmi <- round(checked$largest_fmi, 3)
checked <- checked[, c("analysis", "model", "largest_fmi", "term_with_largest_fmi",
                       "imputations", "imputations_needed", "rule_met")]
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "12_imputation_number_check.csv"),
          row.names = FALSE, na = "")
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "12_imputation_number_check", primary_suppressed = sum(hide),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------


###############################################################################
# Block 8. Comparison with the complete-case analysis (File 5)
#
# For the science terms, side by side: complete cases (File 5), then each
# imputation analysis. Quantities: the three standardized risk differences,
# and the odds ratios of the four age-14 orientation terms. The imputation
# rows come from the checked tables above. The complete-case rows get File
# 5's disclosure rules again, written out here.
# File 5's models are labelled as follows:
#   "M0 baseline (own complete cases)"  M0, complete cases: M0's own sample
#   M0-M4 (common sample)               complete cases: common M0-M4 sample
#   "M0 baseline, employed only"        complete cases
###############################################################################

orientation_variables <- c("science_liking14", "science_competence14",
                           "science_favourite14", "math_liking14")

# --- Disclosure check: complete-case rows (File 5's rules) --------------------
cc_samples <- file5$samples
cc_small_models <- cc_samples$model[cc_samples$events < disclosure_threshold |
                                      cc_samples$non_events < disclosure_threshold]
cc_coefficients <- file5$coefficients
cc_hide_coefficient <- cc_coefficients$model %in% cc_small_models |
  (cc_coefficients$term_is_binary &
     (cc_coefficients$people_with_term < disclosure_threshold |
      cc_coefficients$events_with_term < disclosure_threshold |
      cc_coefficients$non_events_with_term < disclosure_threshold))
cc_risk_differences <- file5$risk_differences
cc_hide_rd <- cc_risk_differences$model %in% cc_small_models
# -----------------------------------------------------------------------------

cc_label <- c("M0 baseline (own complete cases)" = "Complete cases: M0's own sample",
              "M0 baseline"                      = "Complete cases: common M0-M4 sample",
              "M1 + course route (15)"           = "Complete cases: common M0-M4 sample",
              "M2 + intentions (16)"             = "Complete cases: common M0-M4 sample",
              "M3 + HE pathway (20)"             = "Complete cases: common M0-M4 sample",
              "M4 + degree pathway (25)"         = "Complete cases: common M0-M4 sample",
              "M0 baseline, employed only"       = "Complete cases")
# The model each complete-case model corresponds to (M0's own sample is
# still M0).
cc_model_name <- names(cc_label)
cc_model_name[cc_model_name == "M0 baseline (own complete cases)"] <- "M0 baseline"
names(cc_model_name) <- names(cc_label)

comparison <- data.frame()

# Risk differences: complete cases (hidden values set to missing).
cc_rd_rows <- data.frame(
  model    = cc_model_name[cc_risk_differences$model],
  quantity = paste("Risk difference:", cc_risk_differences$contrast),
  analysis = cc_label[cc_risk_differences$model],
  estimate = cc_risk_differences$risk_difference,
  low      = cc_risk_differences$low,
  high     = cc_risk_differences$high
)
cc_rd_rows[cc_hide_rd, c("estimate", "low", "high")] <- NA
comparison <- rbind(comparison, cc_rd_rows)

# Risk differences: imputation analyses.
comparison <- rbind(comparison, data.frame(
  model    = checked_risk_differences$model,
  quantity = paste("Risk difference:", checked_risk_differences$contrast),
  analysis = checked_risk_differences$analysis,
  estimate = checked_risk_differences$risk_difference,
  low      = checked_risk_differences$low,
  high     = checked_risk_differences$high
))

# Odds ratios of the orientation terms: complete cases.
is_orientation <- rep(FALSE, nrow(cc_coefficients))
for (v in orientation_variables) {
  is_orientation <- is_orientation | startsWith(cc_coefficients$term, v)
}
cc_or <- cc_coefficients[is_orientation, ]
cc_or_hidden <- cc_hide_coefficient[is_orientation]
cc_or_rows <- data.frame(
  model    = cc_model_name[cc_or$model],
  quantity = paste("Odds ratio:", cc_or$term),
  analysis = cc_label[cc_or$model],
  estimate = cc_or$odds_ratio,
  low      = cc_or$odds_ratio_low,
  high     = cc_or$odds_ratio_high
)
cc_or_rows[cc_or_hidden, c("estimate", "low", "high")] <- NA
comparison <- rbind(comparison, cc_or_rows)

# Odds ratios of the orientation terms: imputation analyses.
is_orientation <- rep(FALSE, nrow(checked_coefficients))
for (v in orientation_variables) {
  is_orientation <- is_orientation | startsWith(checked_coefficients$term, v)
}
mi_or <- checked_coefficients[is_orientation, ]
comparison <- rbind(comparison, data.frame(
  model    = mi_or$model,
  quantity = paste("Odds ratio:", mi_or$term),
  analysis = mi_or$analysis,
  estimate = mi_or$odds_ratio,
  low      = mi_or$odds_ratio_low,
  high     = mi_or$odds_ratio_high
))

# Order: by model, then quantity, then analysis (complete cases first).
analysis_order <- c("Complete cases: M0's own sample", "Complete cases: common M0-M4 sample",
                    "Complete cases", analyses$analysis)
comparison <- comparison[order(match(comparison$model, models$model),
                               comparison$quantity,
                               match(comparison$analysis, analysis_order)), ]
rownames(comparison) <- NULL
comparison$disclosure_checked <- TRUE
write.csv(comparison, file.path(review_dir, "12_comparison.csv"), row.names = FALSE, na = "")
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "12_comparison (complete-case rows)",
  primary_suppressed = sum(cc_hide_rd) + sum(cc_or_hidden),
  secondary_suppressed = 0L))

disclosure_log$threshold <- disclosure_threshold
write.csv(disclosure_log, file.path(review_dir, "12_disclosure_log.csv"), row.names = FALSE)


###############################################################################
# Block 9. Save the exact results
###############################################################################

saveRDS(
  list(analyses = analyses, models = models, contrasts = contrasts,
       pooled = pooled, prevalence = prevalence, coefficients = coefficients,
       risk_differences = risk_differences, imputation_check = imputation_check,
       results_by_imputation = results_by_imputation),
  file.path(data_dir, "12_pooled_estimates.rds")
)


###############################################################################
# Block 10. Summary
#
# Counts, the imputation-number rule and missing-information summaries
# only. No estimates, which are read from the checked tables.
###############################################################################

primary_check <- imputation_check[imputation_check$analysis == analyses$analysis[1], ]
primary_rows  <- coefficients$analysis == analyses$analysis[1]

cat("\nFile 12 finished.\n")
cat("  Analyses:                              ", nrow(analyses), "\n")
cat("  Models per analysis:                   ", nrow(models), "\n")
cat("  Model fits that did not converge:      ", not_converged, "\n")
cat("  Primary analysis, imputation-number rule (m >= 50 and >= largest FMI %):",
    if (all(primary_check$rule_met)) "met" else "NOT met", "\n")
cat("    largest FMI:", round(max(primary_check$largest_fmi), 3),
    " imputations:", primary_check$imputations[1],
    " needed:", max(primary_check$imputations_needed), "\n")
cat("    Monte Carlo error above 10% of the SE:",
    sum(coefficients$mc_error_vs_se[primary_rows] > 0.1), "of",
    sum(primary_rows), "coefficients\n")
if (!all(primary_check$rule_met)) {
  cat("  More imputations are needed for the primary analysis: add tasks in\n")
  cat("  File 8a, rerun 8b, 10 and this file.\n")
}
cat("  Review: 12_prevalence.csv, 12_coefficients.csv, 12_risk_differences.csv,\n")
cat("    12_imputation_number_check.csv, 12_comparison.csv\n")
