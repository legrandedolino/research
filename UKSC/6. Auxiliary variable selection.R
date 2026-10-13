###############################################################################
# 6. Auxiliary variable selection
#
# Purpose
#   Choose the auxiliary variables for the multiple imputation model: variables
#   that aren't in the analysis model but make the missing-at-random
#   assumption more plausible, and the imputations more precise.
#
# Categories of auxiliary variables
#   Let Y be the incomplete variable of interest. Candidates are sorted into
#   three categories, similar to those described by Collins, Schafer and Kam
#   (2001):
#     Category a  associated with Y and with missingness of Y
#     Category b  associated with Y but not with missingness of Y
#     Category c  not associated with Y (whether or not it relates to
#                 missingness)
#   Collins et al. report that omitting variables like Category a can bias
#   results, that including variables like Category b can improve precision,
#   and that variables like Category c bring little benefit. The CLS
#   missing-data user guide gives similar advice: include variables that
#   predict both missingness and the missing values, and variables that
#   predict the missing values only, and leave out variables that predict
#   missingness only.
#
# Procedure (CLS user guide)
#   Y is the incomplete variable with the most missing data: the age-32 SRET
#   outcome. Candidates are the CLS predictors of non-response (user guide
#   Appendix C, Table C1; Silverwood et al. 2024), which are known correlates
#   of missingness. They are used as published, with one logical update:
#   Appendix C's "previous non-response at Waves 2-7" covers the sweeps
#   before its Wave-8 target, so for a Wave-9 outcome it becomes "previous
#   non-response at Waves 2-8". Of the 21 published predictors, sex is
#   already in the analysis model, and consent to linkage with the
#   Department for Work and Pensions (Wave 4) isn't in the End User Licence
#   data, which leaves 19 candidates.
#   Two tests are run for each candidate:
#     1. Association with Y: a logistic regression of the observed age-32
#        outcome on the candidate, with a joint Wald test across all its
#        categories. This is the selection test: keep if p < 0.05. The guide
#        recommends a more lenient cutoff for a small pool of candidates,
#        and 19 candidates are screened here.
#     2. Association with missingness of Y: the same test for whether the
#        age-32 outcome is missing, in the whole imputation population. This
#        is descriptive. Together with test 1, it assigns each candidate to
#        Category a, b or c.
#   Decision: keep Category a and Category b; drop Category c.
#   How the categories differ from Collins et al.'s: Collins et al. define
#   their categories by population correlations and give no procedure for
#   classifying real variables. Here the categories are assigned with the
#   CLS guide's association screen, i.e. by statistical significance, not
#   effect size. With samples of this size, small associations reach
#   p < 0.05, so "associated" means a detectable association, not
#   necessarily a strong one. The categories are therefore similar to,
#   but not the same as, Collins et al.'s. This fits the CLS guide's view
#   that weakly associated auxiliaries do little harm, while those
#   unrelated to Y are left out.
#   Both tests are unweighted, as in the guide's worked example; the aim is
#   to detect association, not to estimate it. Only joint tests are
#   reported, so no category-level estimates are shown.
#
#   Not screened:
#   - sex, because it's already in the analysis model;
#   - Wave-9 response, because it's constant among people with an observed
#     outcome and defines missingness itself;
#   - the design weight and stratum, which enter every imputation equation
#     as predictors (File 7).
#
# Prespecified Category b auxiliaries (always included)
#   SRET employment at 25 and main activity at 25 are included for their
#   correlation with Y, whether or not they relate to missingness. They are
#   the most direct earlier measures of the age-32 outcome and of employment
#   at 32. As in the CLS guide's example (earlier income as an auxiliary for
#   later income), they improve the imputed values even though they come
#   after the age-14 exposures. Their test results, and so their empirical
#   category, are still reported.
#
# Completeness
#   The share of the imputation population observed on each candidate is
#   reported, because the guide prefers more complete auxiliaries.
#   Incomplete auxiliaries are themselves imputed.
#
# Reads
#   data_dir/03_missingness_and_landmarks.rds
#
# Writes
#   data_dir/06_auxiliary_list.rds         the frozen auxiliary list and the
#                                          full screening table
#   output_dir/private_06_auxiliary_screen.csv   exact screening table
#   review_dir/06_auxiliary_screen.csv     disclosure-checked screening table
#   review_dir/06_auxiliary_list.csv       the frozen list (names and
#                                          categories only)
#   review_dir/06_disclosure_log.csv       what was suppressed
#
# References
#   Collins LM, Schafer JL, Kam C-M (2001). A comparison of inclusive and
#   restrictive strategies in modern missing data procedures. Psychological
#   Methods 6(4): 330-351.
#   Silverwood R, et al. (2026). Handling missing data in the CLS cohort
#   studies: User guide (Version 2). UCL Centre for Longitudinal Studies.
#   Silverwood RJ, et al. (2024). Longitudinal and Life Course Studies 15(2):
#   227-250.
###############################################################################

library(haven)   # needed to read the deposited (SPSS-labelled) variables
suppressPackageStartupMessages(library(dplyr))


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

disclosure_threshold <- 10     # smallest count that may be shown
screen_p_cutoff      <- 0.05   # p-value cutoff for both association tests


###############################################################################
# Block 2. Read the cohort from File 3
###############################################################################

cohort <- readRDS(file.path(data_dir, "03_missingness_and_landmarks.rds"))$cohort

# The imputation population (the target population; File 3), whether each
# person's outcome (Y) is missing, and the people with an observed outcome.
population <- cohort[cohort$in_target_population, ]
population$outcome_missing <- as.integer(!population$outcome_observed32)
screening  <- population[population$outcome_observed32, ]
outcome    <- "sret_professional_32"


###############################################################################
# Block 3. Candidates
#
# One row per candidate, with its group. SRET at 25 (0/1) is turned into a
# factor, so every categorical candidate is handled the same way.
###############################################################################

candidates <- data.frame(
  variable = c(
    # CLS predictors of non-response (Appendix C, Table C1), except sex
    "nr_parents_know_where14", "nr_upset_by_names14", "nr_computer_games14",
    "nr_played_instrument14", "nr_housing_tenure14", "nr_internet_at_home14",
    "nr_school_contacted_parents15", "nr_under_strain15",
    "nr_smokes16", "nr_parent_age16",
    "nr_nightclubs17",
    "nr_same_address18", "nr_groups_treated_better18",
    "nr_teacher_expectations18", "nr_main_activity18",
    "nr_spoke_to_teacher19", "nr_willing_sex_questions19",
    "nr_willing_sex_questions20",
    "nr_previous_nonresponse_w2_w8",   # Appendix C item, updated to Waves 2-8
    # Prespecified Category b auxiliaries
    "sret_professional_25", "activity25"
  ),
  group = c(
    rep("CLS predictor of non-response", 19),
    rep("Prespecified Category b", 2)
  )
)

population$sret_professional_25 <- factor(population$sret_professional_25, levels = 0:1)
screening <- population[population$outcome_observed32, ]


###############################################################################
# Block 4. Completeness in the imputation population
###############################################################################

candidates$n_observed <- NA_integer_
candidates$n_missing  <- NA_integer_
for (i in seq_len(nrow(candidates))) {
  values <- population[[candidates$variable[i]]]
  candidates$n_observed[i] <- sum(!is.na(values))
  candidates$n_missing[i]  <- sum(is.na(values))
}
candidates$percent_observed <- round(100 * candidates$n_observed / nrow(population), 1)


###############################################################################
# Block 5. Test 1: association with Y (the selection test)
#
# For each candidate:
#   fit    outcome ~ candidate, unweighted logistic regression, among people
#          with an observed outcome who are observed on the candidate;
#   test   joint Wald test of all the candidate's coefficients:
#            W = b' V^-1 b, compared with chi-square on (number of
#            coefficients) degrees of freedom;
#   record people used, events and non-events (for the disclosure check),
#          W, degrees of freedom and p.
# If a coefficient can't be estimated, or the variance matrix is singular
# (e.g. a category with no events), the candidate is "not testable".
###############################################################################

candidates$y_n_used          <- NA_integer_
candidates$y_events_used     <- NA_integer_
candidates$y_non_events_used <- NA_integer_
candidates$y_wald_chi2       <- NA_real_
candidates$y_wald_df         <- NA_integer_
candidates$y_p_value         <- NA_real_
candidates$y_testable        <- NA

for (i in seq_len(nrow(candidates))) {
  variable <- candidates$variable[i]
  fit_data <- screening[!is.na(screening[[variable]]), c(outcome, variable)]
  if (is.factor(fit_data[[variable]])) fit_data[[variable]] <- droplevels(fit_data[[variable]])

  fit <- glm(as.formula(paste(outcome, "~", variable)), data = fit_data,
             family = binomial())
  b <- coef(fit)[-1]                       # the candidate's coefficients
  V <- vcov(fit)[-1, -1, drop = FALSE]     # and their variance matrix

  candidates$y_n_used[i]          <- nrow(fit_data)
  candidates$y_events_used[i]     <- sum(fit_data[[outcome]] == 1)
  candidates$y_non_events_used[i] <- sum(fit_data[[outcome]] == 0)

  testable <- length(b) > 0 && !any(is.na(b)) && rcond(V) > 1e-10 && fit$converged
  candidates$y_testable[i] <- testable
  if (testable) {
    W <- as.numeric(t(b) %*% solve(V) %*% b)
    candidates$y_wald_chi2[i] <- round(W, 2)
    candidates$y_wald_df[i]   <- length(b)
    candidates$y_p_value[i]   <- signif(pchisq(W, df = length(b), lower.tail = FALSE), 3)
  }
}


###############################################################################
# Block 6. Test 2: association with missingness of Y (descriptive)
#
# The same test for outcome_missing (1 = age-32 outcome missing), in the
# whole imputation population, among people observed on the candidate. It
# doesn't change which candidates are kept; it separates Category a from
# Category b. The "events" here are people with a missing outcome.
###############################################################################

candidates$r_n_used          <- NA_integer_
candidates$r_missing_used    <- NA_integer_
candidates$r_observed_used   <- NA_integer_
candidates$r_wald_chi2       <- NA_real_
candidates$r_wald_df         <- NA_integer_
candidates$r_p_value         <- NA_real_
candidates$r_testable        <- NA

for (i in seq_len(nrow(candidates))) {
  variable <- candidates$variable[i]
  fit_data <- population[!is.na(population[[variable]]), c("outcome_missing", variable)]
  if (is.factor(fit_data[[variable]])) fit_data[[variable]] <- droplevels(fit_data[[variable]])

  fit <- glm(as.formula(paste("outcome_missing ~", variable)), data = fit_data,
             family = binomial())
  b <- coef(fit)[-1]
  V <- vcov(fit)[-1, -1, drop = FALSE]

  candidates$r_n_used[i]        <- nrow(fit_data)
  candidates$r_missing_used[i]  <- sum(fit_data$outcome_missing == 1)
  candidates$r_observed_used[i] <- sum(fit_data$outcome_missing == 0)

  testable <- length(b) > 0 && !any(is.na(b)) && rcond(V) > 1e-10 && fit$converged
  candidates$r_testable[i] <- testable
  if (testable) {
    W <- as.numeric(t(b) %*% solve(V) %*% b)
    candidates$r_wald_chi2[i] <- round(W, 2)
    candidates$r_wald_df[i]   <- length(b)
    candidates$r_p_value[i]   <- signif(pchisq(W, df = length(b), lower.tail = FALSE), 3)
  }
}


###############################################################################
# Block 7. Categories, decisions and the frozen auxiliary list
#
#   Category a  associated with Y (p < 0.05) and with missingness of Y
#   Category b  associated with Y, not with missingness of Y
#   Category c  not associated with Y
#   Decision    keep Category a and b; drop Category c and anything not
#               testable. The prespecified Category b auxiliaries are always
#               kept; their empirical category is still recorded.
###############################################################################

associated_with_y           <- candidates$y_testable & candidates$y_p_value < screen_p_cutoff
associated_with_missingness <- candidates$r_testable & candidates$r_p_value < screen_p_cutoff

candidates$auxiliary_category <- case_when(
  !candidates$y_testable                        ~ "Not testable",
  associated_with_y & associated_with_missingness ~ "a",
  associated_with_y                             ~ "b",
  TRUE                                          ~ "c"
)

is_prespecified <- candidates$group == "Prespecified Category b"
candidates$decision <- case_when(
  is_prespecified                               ~ "Keep (prespecified Category b)",
  candidates$auxiliary_category == "a"            ~ "Keep (Category a)",
  candidates$auxiliary_category == "b"            ~ "Keep (Category b)",
  candidates$auxiliary_category == "c"            ~ "Drop (Category c)",
  TRUE                                          ~ "Drop (not testable)"
)

auxiliary_list <- candidates$variable[startsWith(candidates$decision, "Keep")]

write.csv(candidates, file.path(output_dir, "private_06_auxiliary_screen.csv"),
          row.names = FALSE)


###############################################################################
# Block 8. Disclosure check and review tables
###############################################################################

# --- Disclosure check: screening table ---------------------------------------
# Completeness: if fewer than 10 people are missing (or observed), both
# counts and the percentage are hidden, since either count gives the other
# from the population total.
# Each test: its counts and statistics are hidden unless both of its groups
# (events and non-events for test 1; missing and observed outcome for
# test 2) have at least 10 people. Categories and decisions are always
# shown; they reveal nothing about small groups.
checked <- candidates
small_completeness <- (checked$n_missing  > 0 & checked$n_missing  < disclosure_threshold) |
                      (checked$n_observed > 0 & checked$n_observed < disclosure_threshold)
small_y <- checked$y_events_used < disclosure_threshold |
           checked$y_non_events_used < disclosure_threshold
small_r <- checked$r_missing_used < disclosure_threshold |
           checked$r_observed_used < disclosure_threshold
checked[small_completeness, c("n_observed", "n_missing", "percent_observed")] <- NA
checked[small_y, c("y_n_used", "y_events_used", "y_non_events_used",
                   "y_wald_chi2", "y_wald_df", "y_p_value")] <- NA
checked[small_r, c("r_n_used", "r_missing_used", "r_observed_used",
                   "r_wald_chi2", "r_wald_df", "r_p_value")] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "06_auxiliary_screen.csv"), row.names = FALSE)

# The frozen list holds names and categories only; nothing to suppress.
write.csv(data.frame(auxiliary        = auxiliary_list,
                     group            = candidates$group[match(auxiliary_list, candidates$variable)],
                     auxiliary_category = candidates$auxiliary_category[match(auxiliary_list, candidates$variable)],
                     disclosure_checked = TRUE),
          file.path(review_dir, "06_auxiliary_list.csv"), row.names = FALSE)

disclosure_log <- data.frame(
  table = "06_auxiliary_screen",
  primary_suppressed = sum(small_completeness) + sum(small_y) + sum(small_r),
  secondary_suppressed = 0L,
  threshold = disclosure_threshold
)
write.csv(disclosure_log, file.path(review_dir, "06_disclosure_log.csv"), row.names = FALSE)
# -----------------------------------------------------------------------------


###############################################################################
# Block 9. Save
###############################################################################

saveRDS(list(auxiliary_list = auxiliary_list, screen = candidates,
             p_cutoff = screen_p_cutoff),
        file.path(data_dir, "06_auxiliary_list.rds"))


###############################################################################
# Block 10. Summary
#
# Numbers of candidates by category and decision only.
###############################################################################

screened <- !is_prespecified
cat("\nFile 6 finished.\n")
cat("  Candidates screened:                 ", sum(screened), "\n")
cat("    Category a (kept):                 ", sum(screened & candidates$auxiliary_category == "a"), "\n")
cat("    Category b (kept):                 ", sum(screened & candidates$auxiliary_category == "b"), "\n")
cat("    Category c (dropped):              ", sum(screened & candidates$auxiliary_category == "c"), "\n")
cat("    Not testable (dropped):            ", sum(screened & candidates$auxiliary_category == "Not testable"), "\n")
cat("  Prespecified Category b (kept):      ", sum(is_prespecified), "\n")
cat("  Auxiliaries in the frozen list:      ", length(auxiliary_list), "\n")
cat("  Saved:", file.path(data_dir, "06_auxiliary_list.rds"), "\n")
