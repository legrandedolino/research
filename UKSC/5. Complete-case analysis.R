###############################################################################
# 5. Complete-case analysis
#
# Purpose
#   The first step of the missing-data analysis, as recommended in the CLS
#   missing-data user guide:
#     Table 1  compare complete and incomplete cases on every variable in the
#              baseline analysis model
#     Table 2  analysis samples: size, events and weighted SRET prevalence
#     Table 3  survey-weighted logistic models fitted to complete cases:
#              the baseline model, the sequential models that add each
#              later educational stage, and the employed-only version
#     Table 4  standardized risk differences for the age-14 science
#              contrasts in every model
#   These complete-case results are the reference point for the multiple
#   imputation results (File 12). A complete-case analysis is valid if data
#   are missing completely at random, and in some situations beyond that;
#   differences from the imputation results suggest bias from missing data.
#
# Analysis models (prespecified; the same models are used in File 12)
#   Outcome: professional SRET employment at 32 (sret_professional_32).
#   Baseline model (M0): age-14 science liking, perceived science competence
#     and science as a favourite subject, with mathematics liking as a
#     comparator, adjusted for sex, ethnic group, family NS-SEC, highest
#     resident-parent qualification and household income band.
#   Sequential models, each adding one later stage to the one before:
#     M1 + science course route (15)
#     M2 + continuation in full-time education, likelihood of applying to
#          university, intended HE subject route (16)
#     M3 + HE pathway by 20
#     M4 + degree pathway by 25
#   M0-M4 are fitted to one common set of complete cases (complete on every
#   M4 variable), so changes in the science associations across models
#   reflect the added stages, not a changing sample. M0 is also fitted to its
#   own, larger set of complete cases, which is the main complete-case
#   reference. Employed-only sensitivity: M0 with SRET among the employed
#   (sret_among_employed32) as the outcome.
#   The associations are descriptive, not causal effects.
#
# Survey design
#   All models and prevalences use the Wave-9 final weight (W9FINWT), PSU and
#   stratum. Complete cases are analysed as a subset of the full design, so
#   the standard errors reflect the whole sample design.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   Exact tables go to private files. Checked copies go to review/:
#     - counts 1-9 suppressed, with secondary suppression within each group;
#     - means and SDs shown only for groups of at least 10;
#     - a model is shown only if it has at least 10 events and 10
#       non-events;
#     - a coefficient for a category (or any 0/1 term) is shown only if at
#       least 10 people in the model sample have it, with at least 10 events
#       and 10 non-events among them. This prevents disclosures such as "no
#       one in this group is in SRET".
#
# Reads
#   data_dir/03_missingness_and_landmarks.rds
#
# Writes
#   data_dir/05_complete_case_results.rds   exact result tables (for File 12)
#   output_dir/private_05_*.csv             exact tables (never shared)
#   review_dir/05_*.csv                     disclosure-checked tables
#   review_dir/05_disclosure_log.csv        what was suppressed
#
# References
#   Silverwood R, et al. (2026). Handling missing data in the CLS cohort
#   studies: User guide (Version 2). UCL Centre for Longitudinal Studies.
#   White IR, Carlin JB (2010). Bias and efficiency of multiple imputation
#   compared with complete-case analysis for missing covariate values.
#   Statistics in Medicine 29(28): 2920-2931.
###############################################################################

library(haven)   # needed to read the deposited (SPSS-labelled) variables
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(survey))


###############################################################################
# Block 1. Paths and disclosure threshold
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

# Standard errors compare PSUs within each stratum. A stratum with only one
# PSU in the data the design is built from (e.g. after attrition) has nothing
# to compare with, and the survey package stops with an error by default.
# "adjust" instead measures that PSU's deviation from the overall mean, which
# is slightly conservative. Subsets made with subset() (e.g. complete cases)
# keep the full design's PSUs and use the standard subpopulation formulas.
options(survey.lonely.psu = "adjust")


###############################################################################
# Block 2. Read the cohort from File 3
###############################################################################

cohort <- readRDS(file.path(data_dir, "03_missingness_and_landmarks.rds"))$cohort

# The HE-likelihood item is an ordered factor. It enters the models as
# ordinary categories (each compared with the first), so its order doesn't
# change the coding.
cohort$he_likelihood16 <- factor(cohort$he_likelihood16, ordered = FALSE)


###############################################################################
# Block 3. Model variables and complete cases
#
# Complete cases are people in the prediction cohort (productive Wave-9
# interview, observed outcome, Wave-9 weight, PSU and stratum; defined in
# File 3) who are observed on every variable of the model concerned.
###############################################################################

outcome <- "sret_professional_32"

orientation_terms <- c("science_liking14", "science_competence14",
                       "science_favourite14", "math_liking14")
background_terms  <- c("sex14", "ethnicity14", "nssec_family14",
                       "highest_parent_education14", "household_income14")
route15_terms     <- "science_course_route15"
intention16_terms <- c("w3_continue_fte16", "he_likelihood16",
                       "w3_intention_route16")
he20_terms        <- "he_pathway20"
degree25_terms    <- "degree_pathway25"

baseline_terms   <- c(orientation_terms, background_terms)
sequential_terms <- c(baseline_terms, route15_terms, intention16_terms,
                      he20_terms, degree25_terms)

cohort$complete_baseline <- cohort$in_prediction_cohort &
  complete.cases(cohort[c(outcome, baseline_terms)])
cohort$complete_sequential <- cohort$in_prediction_cohort &
  complete.cases(cohort[c(outcome, sequential_terms)])
cohort$complete_employed <- cohort$in_prediction_cohort &
  complete.cases(cohort[c("sret_among_employed32", baseline_terms)])


###############################################################################
# Block 4. Table 1: complete versus incomplete cases
#
# Within the target population (the imputation population; File 3), the
# distribution of each baseline-model variable among complete cases
# (complete_baseline) and among everyone else, each among those observed on
# that variable. Categorical variables: n and % per category. Numeric
# variables: n, mean and SD. Unweighted, as in the CLS user guide.
# Systematic differences suggest the data aren't missing completely at
# random.
###############################################################################

population <- cohort[cohort$in_target_population, ]
population$case_group <- ifelse(population$complete_baseline,
                                "Complete case", "Incomplete case")
population$sret_outcome_32 <- factor(population[[outcome]], levels = 0:1,
                                     labels = c("Not SRET", "SRET"))

comparison_variables <- c("sret_outcome_32", baseline_terms)
comparison <- data.frame()

# The rows to write, in order: for each variable and group, one row per
# category of a categorical variable, or one "(numeric)" row.
row_table <- data.frame()
for (variable in comparison_variables) {
  categories <- "(numeric)"
  if (!is.numeric(population[[variable]])) categories <- levels(population[[variable]])
  row_table <- rbind(row_table, expand.grid(category = categories,
                                            group = c("Complete case", "Incomplete case"),
                                            variable = variable, stringsAsFactors = FALSE))
}

for (k in seq_len(nrow(row_table))) {
  variable <- row_table$variable[k]
  group    <- row_table$group[k]
  category <- row_table$category[k]
  values <- population[[variable]][population$case_group == group &
                                     !is.na(population[[variable]])]
  if (category == "(numeric)") {
    comparison <- rbind(comparison, data.frame(
      variable = variable, group = group, category = "(numeric)",
      n = length(values), percent = NA,
      mean = round(mean(values), 2), sd = round(sd(values), 2)))
  } else {
    comparison <- rbind(comparison, data.frame(
      variable = variable, group = group, category = category,
      n = sum(values == category),
      percent = round(100 * sum(values == category) / length(values), 1),
      mean = NA, sd = NA))
  }
}

write.csv(comparison, file.path(output_dir, "private_05_complete_vs_incomplete.csv"),
          row.names = FALSE)

# --- Disclosure check: Table 1 -----------------------------------------------
# Categories: counts 1-9 suppressed, with secondary suppression within each
# variable x group. Numeric rows: n, mean and SD shown only if n >= 10.
checked <- comparison
is_numeric_row <- checked$category == "(numeric)"
primary <- !is_numeric_row & checked$n > 0 & checked$n < disclosure_threshold
secondary <- rep(FALSE, nrow(checked))
group_id <- paste(checked$variable, checked$group)
for (g in unique(group_id[!is_numeric_row])) {
  rows <- which(group_id == g)
  if (sum(primary[rows]) == 1) {
    candidates <- rows[!primary[rows] & checked$n[rows] > 0]
    if (length(candidates) > 0) {
      secondary[candidates[which.min(checked$n[candidates])]] <- TRUE
    }
  }
}
numeric_small <- is_numeric_row & checked$n < disclosure_threshold
checked$n[primary | secondary]       <- NA
checked$percent[primary | secondary] <- NA
checked[numeric_small, c("n", "mean", "sd")] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "05_complete_vs_incomplete.csv"),
          row.names = FALSE)
disclosure_log <- data.frame(
  table = "05_complete_vs_incomplete",
  primary_suppressed = sum(primary) + sum(numeric_small),
  secondary_suppressed = sum(secondary))
# -----------------------------------------------------------------------------


###############################################################################
# Block 5. Survey design
#
# One design for the prediction cohort. Each model uses a subset of it
# (subset() keeps the full design information for the standard errors).
###############################################################################

prediction_cohort <- cohort[cohort$in_prediction_cohort, ]
design <- svydesign(ids = ~psu, strata = ~stratum, weights = ~w9_final_weight,
                    nest = TRUE, data = prediction_cohort)


###############################################################################
# Block 6. The models
#
# One row per model: its name, outcome, complete-case sample and
# predictors. Block 7 fits each one in turn.
###############################################################################

models <- data.frame(
  model = c("M0 baseline (own complete cases)",
            "M0 baseline", "M1 + course route (15)", "M2 + intentions (16)",
            "M3 + HE pathway (20)", "M4 + degree pathway (25)",
            "M0 baseline, employed only"),
  outcome = c(rep(outcome, 6), "sret_among_employed32"),
  sample  = c("complete_baseline", rep("complete_sequential", 5),
              "complete_employed"),
  predictors = c(
    paste(baseline_terms, collapse = " + "),
    paste(baseline_terms, collapse = " + "),
    paste(c(baseline_terms, route15_terms), collapse = " + "),
    paste(c(baseline_terms, route15_terms, intention16_terms), collapse = " + "),
    paste(c(baseline_terms, route15_terms, intention16_terms, he20_terms),
          collapse = " + "),
    paste(sequential_terms, collapse = " + "),
    paste(baseline_terms, collapse = " + ")
  )
)

# The science contrasts for the standardized risk differences: each person's
# predicted probability is computed with the variable set to "high" and then
# to "low", and the weighted averages are compared.
contrasts <- data.frame(
  contrast = c("Science liking: like it a lot vs don't like it at all",
               "Perceived science competence: very good vs no good at all",
               "Science is a favourite subject: yes vs no"),
  variable = c("science_liking14", "science_competence14", "science_favourite14"),
  high     = c("4", "4", "Science is a favourite subject"),
  low      = c("1", "1", "Science is not a favourite subject")
)


###############################################################################
# Block 7. Fit each model and compute its results
#
# For each model:
#   a. fit a survey-weighted logistic regression (quasibinomial, which is
#      the usual choice with survey weights) to its complete cases, allowing
#      up to 100 iterations, and record whether it converged. A model that
#      doesn't converge is reported in the console and flagged in Table 2,
#      because its estimates aren't reliable;
#   b. record sample size, events, non-events and weighted prevalence;
#   c. record every coefficient with its odds ratio and 95% CI (t
#      distribution with the design's degrees of freedom: number of PSUs
#      minus number of strata), p-value, and, for the
#      disclosure check, how many people have each term and how many of them
#      are events;
#   d. compute the standardized risk difference for each science contrast:
#        RD = weighted mean of p(high) - weighted mean of p(low)
#      with a delta-method standard error from the model's variance matrix.
#      This treats the covariate distribution of the sample as fixed.
###############################################################################

samples      <- data.frame()
coefficients <- data.frame()
risk_differences <- data.frame()

for (m in seq_len(nrow(models))) {
  model_name <- models$model[m]
  model_outcome <- models$outcome[m]

  # a. Fit.
  sub_design <- subset(design, prediction_cohort[[models$sample[m]]])
  model_formula <- as.formula(paste(model_outcome, "~", models$predictors[m]))
  fit <- svyglm(model_formula, design = sub_design, family = quasibinomial(),
                control = glm.control(maxit = 100))
  if (!fit$converged) {
    warning(model_name, ": the model did not converge; its estimates are ",
            "not reliable.")
  }

  if (any(is.na(coef(fit)))) {
    stop(model_name, ": some terms can't be estimated (perfectly collinear): ",
         paste(names(coef(fit))[is.na(coef(fit))], collapse = ", "))
  }

  y <- sub_design$variables[[model_outcome]]
  w <- weights(sub_design)

  # b. Sample and weighted prevalence.
  prevalence <- svymean(as.formula(paste("~", model_outcome)), sub_design)
  samples <- rbind(samples, data.frame(
    model = model_name, outcome = model_outcome,
    n = length(y), events = sum(y == 1), non_events = sum(y == 0),
    weighted_percent = round(100 * coef(prevalence)[[1]], 2),
    weighted_percent_low  = round(100 * (coef(prevalence)[[1]] - 1.96 * SE(prevalence)[[1]]), 2),
    weighted_percent_high = round(100 * (coef(prevalence)[[1]] + 1.96 * SE(prevalence)[[1]]), 2),
    design_df = degf(sub_design),
    converged = fit$converged
  ))

  # c. Coefficients.
  estimates <- summary(fit)$coefficients
  critical  <- qt(0.975, df = degf(sub_design))
  X <- model.matrix(fit)
  coefficients <- rbind(coefficients, data.frame(
    model     = model_name,
    term      = rownames(estimates),
    log_odds  = round(estimates[, "Estimate"], 4),
    se        = round(estimates[, "Std. Error"], 4),
    odds_ratio      = round(exp(estimates[, "Estimate"]), 3),
    odds_ratio_low  = round(exp(estimates[, "Estimate"] - critical * estimates[, "Std. Error"]), 3),
    odds_ratio_high = round(exp(estimates[, "Estimate"] + critical * estimates[, "Std. Error"]), 3),
    p_value   = signif(estimates[, "Pr(>|t|)"], 3),
    # Disclosure information (not shown in the checked table): for 0/1
    # terms, the people who have the term, and the events and non-events
    # among them.
    term_is_binary = colSums(X != 0 & X != 1) == 0,
    people_with_term     = colSums(X != 0),
    events_with_term     = colSums((X != 0) * (y == 1)),
    non_events_with_term = colSums((X != 0) * (y == 0)),
    row.names = NULL
  ))

  # d. Standardized risk differences.
  for (k in seq_len(nrow(contrasts))) {
    variable <- contrasts$variable[k]
    data_high <- sub_design$variables
    data_low  <- sub_design$variables
    if (is.numeric(data_high[[variable]])) {
      data_high[[variable]] <- as.numeric(contrasts$high[k])
      data_low[[variable]]  <- as.numeric(contrasts$low[k])
    } else {
      data_high[[variable]] <- factor(contrasts$high[k], levels = levels(data_high[[variable]]))
      data_low[[variable]]  <- factor(contrasts$low[k],  levels = levels(data_low[[variable]]))
    }

    predictor_terms <- delete.response(terms(fit))
    X_high <- model.matrix(predictor_terms, data = data_high, xlev = fit$xlevels)
    X_low  <- model.matrix(predictor_terms, data = data_low,  xlev = fit$xlevels)
    p_high <- as.vector(plogis(X_high %*% coef(fit)))
    p_low  <- as.vector(plogis(X_low  %*% coef(fit)))

    rd <- sum(w * (p_high - p_low)) / sum(w)
    gradient <- colSums(w * (p_high * (1 - p_high) * X_high -
                             p_low  * (1 - p_low)  * X_low)) / sum(w)
    rd_se <- sqrt(as.numeric(t(gradient) %*% vcov(fit) %*% gradient))

    risk_differences <- rbind(risk_differences, data.frame(
      model = model_name, contrast = contrasts$contrast[k],
      risk_difference = round(rd, 4), se = round(rd_se, 4),
      low  = round(rd - critical * rd_se, 4),
      high = round(rd + critical * rd_se, 4)
    ))
  }
}

write.csv(samples, file.path(output_dir, "private_05_samples.csv"), row.names = FALSE)
write.csv(coefficients, file.path(output_dir, "private_05_coefficients.csv"),
          row.names = FALSE)
write.csv(risk_differences, file.path(output_dir, "private_05_risk_differences.csv"),
          row.names = FALSE)


###############################################################################
# Block 8. Disclosure checks: Tables 2-4
###############################################################################

# --- Disclosure check: Table 2 (samples) -------------------------------------
# A model's row is shown only if it has at least 10 events and 10 non-events.
checked <- samples
model_too_small <- checked$events < disclosure_threshold |
                   checked$non_events < disclosure_threshold
checked[model_too_small, setdiff(names(checked), c("model", "outcome"))] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "05_samples.csv"), row.names = FALSE)
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "05_samples", primary_suppressed = sum(model_too_small),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------

# --- Disclosure check: Table 3 (coefficients) --------------------------------
# All coefficients of a model that fails the sample check are hidden. A 0/1
# term (a category, or science as favourite subject) is hidden unless at
# least 10 people have it, with at least 10 events and 10 non-events among
# them. The counts used for this check are removed from the shared table.
checked <- coefficients
hidden_model <- checked$model %in% samples$model[model_too_small]
hidden_term  <- checked$term_is_binary &
  (checked$people_with_term < disclosure_threshold |
   checked$events_with_term < disclosure_threshold |
   checked$non_events_with_term < disclosure_threshold)
hide <- hidden_model | hidden_term
checked[hide, c("log_odds", "se", "odds_ratio", "odds_ratio_low",
                "odds_ratio_high", "p_value")] <- NA
checked <- checked[, c("model", "term", "log_odds", "se", "odds_ratio",
                       "odds_ratio_low", "odds_ratio_high", "p_value")]
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "05_coefficients.csv"), row.names = FALSE)
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "05_coefficients", primary_suppressed = sum(hide),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------

# --- Disclosure check: Table 4 (risk differences) ----------------------------
# Hidden for any model that fails the sample check.
checked <- risk_differences
hide <- checked$model %in% samples$model[model_too_small]
checked[hide, c("risk_difference", "se", "low", "high")] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "05_risk_differences.csv"), row.names = FALSE)
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "05_risk_differences", primary_suppressed = sum(hide),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------

disclosure_log$threshold <- disclosure_threshold
write.csv(disclosure_log, file.path(review_dir, "05_disclosure_log.csv"),
          row.names = FALSE)


###############################################################################
# Block 9. Save the exact results for comparison with multiple imputation
###############################################################################

saveRDS(
  list(models = models, contrasts = contrasts, samples = samples,
       coefficients = coefficients, risk_differences = risk_differences,
       complete_vs_incomplete = comparison),
  file.path(data_dir, "05_complete_case_results.rds")
)


###############################################################################
# Block 10. Summary
#
# Counts of models and suppressions only.
###############################################################################

cat("\nFile 5 finished.\n")
cat("  Models fitted:                          ", nrow(models), "\n")
cat("  Tables written (checked, in review/):   ", nrow(disclosure_log), "\n")
cat("  Cells suppressed (primary / secondary): ",
    sum(disclosure_log$primary_suppressed), "/",
    sum(disclosure_log$secondary_suppressed), "\n")
cat("  Saved:", file.path(data_dir, "05_complete_case_results.rds"), "\n")
