###############################################################################
# 7. Imputation specification and test run
#
# Purpose
#   Build the multiple imputation model for the target population (the
#   original cohort, less those out of scope at Wave 9; File 3), check it,
#   and run the CLS user guide's test imputation before the main run
#   (File 8a). The model is fitted to everyone in the target population (CLS
#   user guide: impute using everyone; the analysis samples are chosen
#   afterwards, in File 12).
#
# Survey weights and the imputation model (a deliberate departure)
#   Silverwood et al. (2024) weighted their imputation models by the initial
#   design weights. That is not done here. The imputation uses R's mice
#   package, which has no weighted imputation methods for categorical
#   variables, and most of the variables imputed here are categorical.
#   Instead:
#     - the log of the initial design weight and the sampling stratum are
#       included as predictors in every imputation equation (Block 7). This
#       is a standard way of carrying survey design information into an
#       imputation model (Reiter, Raghunathan and Kinney 2006). As a
#       predictor, strata too small for the imputation models are pooled
#       into one group (Block 4); the analysis uses the original strata;
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
# What the imputation model contains (CLS user guide)
#   - Every variable in the analysis models (File 5), including the outcome,
#     in the form it enters the analysis. There are no interactions in the
#     analysis models, so none are needed here.
#   - The auxiliary variables chosen in File 6.
#   - The log design weight and the sampling stratum.
#   Variables used only in the prediction analysis (e.g. home language,
#   school type, activity histories) are not included. The prediction
#   analysis imputes within its own training folds (File 17).
#
# Routed questions
#   Some questions apply only given an earlier answer:
#     the HE route at 16              only for those staying on in full-time
#                                     education
#     the intended HE subject at 16   only if a subject area was given
#     the HE subject at 20            only if in HE by 20
#     the degree subject at 25        only for degree holders
#     SRET occupation at 32           only for those in paid work
#   Each is imputed in parts that follow the questionnaire's routing:
#     1. the earlier answer (e.g. in paid work at 32 or not), imputed for
#        everyone for whom it is missing;
#     2. the routed answer (e.g. SRET or other occupation), imputed only
#        where it could apply. It is never imputed where the observed
#        earlier answer makes it not applicable (the "structural" state
#        from File 3).
#   At 16 there are three parts: staying on in full-time education, then
#   the HE route for those staying on, then the subject domain where a
#   subject area was given (File 2).
#   The variable used in the analysis is then calculated from its parts at
#   every iteration ("passive imputation"; van Buuren 2018), so the
#   completed data always respect the routing. For the outcome:
#     sret_professional_32 = SRET among the employed if in paid work, else 0.
#   Imputing a routed answer only among the people it applies to also
#   removes near-certain links between the parts (e.g. leaving education and
#   having no university route). Those links would otherwise form a
#   feedback loop that slows convergence (van Buuren 2018, section 6.5.2).
#   A routed part imputed for someone whose imputed earlier answer doesn't
#   lead to it has no effect, because the calculation ignores it. Analyses
#   use the calculated variables.
#
# Imputation methods (one per variable, chosen by its type; Block 6)
#   numeric                       predictive mean matching (pmm)
#   two categories                logistic regression, bootstrap version
#                                 (logreg.boot)
#   more than two categories      classification tree (cart)
#   calculated from other parts   passive (see above)
#   Why not mice's default logistic, proportional-odds and multinomial
#   models: to prevent perfect prediction, they add pseudo-observations to
#   the data (White, Daniel and Royston 2010) with a total weight equal to
#   the number of predictor columns plus one, shared equally between the
#   categories. This model has over 100 predictor columns (the strata alone
#   give over 30). Against a rare outcome, or a rare category, that added
#   weight is large: it pulls imputations toward equal shares of the
#   categories and weakens their associations. Rare categories and the rare
#   SRET outcomes are then imputed too often. In simulations of this
#   setting (about 100 predictor columns, a few hundred events), the
#   default logistic model imputed a rare outcome about 40% too often among
#   the missing, while the methods below were close to the truth. The
#   methods used instead add nothing:
#     - The bootstrap logistic model refits the model to a bootstrap sample
#       of the observed rows at each step, which carries the uncertainty in
#       the model's coefficients into the imputations (van Buuren 2018). If
#       a category of a predictor were missing from a bootstrap sample, its
#       coefficient couldn't be estimated and the imputations would fail, so
#       every category of every categorical predictor must have at least 20
#       people among the rows each such model is fitted to (checked before
#       the test run, Block 8). Sparse categories are collapsed where
#       needed (Block 4).
#     - Classification trees (Burgette and Reiter 2010; Doove, van Buuren
#       and Dusseldorp 2014) impute from the observed values of similar
#       people, so rare categories and perfect prediction need no special
#       handling. The CLS user guide names them as an alternative to
#       collapsing sparse categories. mice's tree method doesn't draw its
#       own parameters, so the variation between imputations may be slightly
#       understated for these variables.
#
# Test run (CLS user guide)
#   5 imputations of 30 iterations each, with a fixed seed. It checks that
#   the model runs and that the completed data respect the routing, and it
#   produces the diagnostics the guide asks for:
#     - trace plots of the mean and SD of imputed values over iterations,
#       with a numeric convergence table (split R-hat and trend);
#     - observed-versus-imputed distributions;
#     - mice's logged events (e.g. predictors removed as collinear) and any
#       warnings.
#   It also gives a first estimate of the fraction of missing information
#   (FMI) in the analysis models. The main run then uses at least 50
#   imputations, and at least the largest FMI as a percentage (CLS user
#   guide). File 12 checks this rule again on the main imputations.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   Exact tables and plots go to private files. Checked copies go to review/:
#     - category percentages are hidden where fewer than 10 people are
#       observed in the category, or fewer than 10 values are imputed into
#       it, with secondary suppression within each variable. Counts are not
#       shown, so they can't be combined with File 5's complete/incomplete
#       counts;
#     - means, SDs and quartiles are shown only for groups of at least 10;
#     - trace plots are shown only for variables with at least 10 imputed
#       values;
#     - FMI is shown only for models with at least 10 events and 10
#       non-events in every completed dataset.
#
# Reads
#   data_dir/03_missingness_and_landmarks.rds
#   data_dir/05_complete_case_results.rds   (the analysis models)
#   data_dir/06_auxiliary_list.rds
#
# Writes
#   data_dir/07_imputation_spec.rds         the imputation model for Files 8a,
#                                           9a, 11a and 13a (participant level)
#   data_dir/07_test_imputation.rds         the test imputation (participant
#                                           level)
#   output_dir/private_07_*                 exact tables and all trace plots
#   review_dir/07_imputation_model.csv      variables, roles and methods
#   review_dir/07_predictor_matrix.csv      which variables predict which
#   review_dir/07_routing_checks.csv        pass/fail checks on the test run
#   review_dir/07_trace_plots.pdf           trace plots (checked)
#   review_dir/07_convergence.csv           split R-hat and trend per variable
#                                           (checked)
#   review_dir/07_observed_vs_imputed.csv   distributions (checked)
#   review_dir/07_logged_events.csv         mice's logged events
#   review_dir/07_warnings.csv              warnings from the test run
#   review_dir/07_fmi.csv                   preliminary FMI (checked)
#   review_dir/07_disclosure_log.csv        what was suppressed
#
# References
#   Reiter JP, Raghunathan TE, Kinney SK (2006). The importance of modeling
#   the sampling design in multiple imputation for missing data. Survey
#   Methodology 32(2): 143-149.
#   Seaman SR, White IR, Copas AJ, Li L (2012). Combining multiple
#   imputation and inverse-probability weighting. Biometrics 68(1): 129-137.
#   Oberman HI, van Buuren S, Vink G (2021). Missing the point: non-convergence
#   in iterative imputation algorithms. arXiv:2110.11951.
#   Vehtari A, Gelman A, Simpson D, Carpenter B, Burkner P-C (2021).
#   Rank-normalization, folding, and localization: an improved R-hat for
#   assessing convergence of MCMC. Bayesian Analysis 16(2): 667-718.
#   Burgette LF, Reiter JP (2010). Multiple imputation for missing data via
#   sequential regression trees. American Journal of Epidemiology 172(9):
#   1070-1076.
#   Doove LL, van Buuren S, Dusseldorp E (2014). Recursive partitioning for
#   missing data imputation in the presence of interaction effects.
#   Computational Statistics & Data Analysis 72: 92-104.
#   Silverwood R, et al. (2026). Handling missing data in the CLS cohort
#   studies: User guide (Version 2). UCL Centre for Longitudinal Studies.
#   Silverwood RJ, et al. (2024). Longitudinal and Life Course Studies 15(2):
#   227-250.
#   van Buuren S (2018). Flexible Imputation of Missing Data (2nd ed.).
#   Chapman & Hall/CRC.
#   White IR, Daniel R, Royston P (2010). Avoiding bias due to perfect
#   prediction in multiple imputation of incomplete categorical variables.
#   Computational Statistics & Data Analysis 54(10): 2267-2275.
###############################################################################

library(haven)   # needed to read the deposited (SPSS-labelled) variables
suppressPackageStartupMessages(library(dplyr))
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

disclosure_threshold <- 10     # smallest count that may be shown

test_m          <- 5       # imputations in the test run
iterations      <- 30      # iterations per imputation (test run and main run).
                           # The CLS guide suggests 10-20 or more. The age-16
                           # questions are missing together and closely
                           # related, so their chains take about 10-15
                           # iterations to settle; 30 leave a margin (van
                           # Buuren 2018, section 6.5.2).
test_seed       <- 7007    # fixed seed for the test run
minimum_m       <- 50      # the main run uses at least this many imputations
minimum_category_size <- 20  # people per predictor category among the rows
                             # each bootstrap logistic model is fitted to

# Standard errors compare PSUs within each stratum. A stratum with only one
# PSU in the data the design is built from (e.g. after attrition) has nothing
# to compare with, and the survey package stops with an error by default.
# "adjust" instead measures that PSU's deviation from the overall mean, which
# is slightly conservative. Subsets made with subset() (e.g. the employed
# only) keep the full design's PSUs and use the standard subpopulation
# formulas.
options(survey.lonely.psu = "adjust")


###############################################################################
# Block 2. Read the cohort, the analysis models and the auxiliary list
###############################################################################

cohort          <- readRDS(file.path(data_dir, "03_missingness_and_landmarks.rds"))$cohort
analysis_models <- readRDS(file.path(data_dir, "05_complete_case_results.rds"))$models
auxiliary_list  <- readRDS(file.path(data_dir, "06_auxiliary_list.rds"))$auxiliary_list

# The imputation population: the target population (File 3).
population <- cohort[cohort$in_target_population, ]


###############################################################################
# Block 3. Variables in the imputation model
#
# One row per variable, with its role. The order is the order in which mice
# visits the variables: design, auxiliaries, then the analysis variables by
# age. Each calculated pathway comes straight after its parts.
###############################################################################

analysis_variables <- data.frame(
  variable = c(
    # Age 14: background
    "sex14", "ethnicity14", "nssec_family14", "highest_parent_education14",
    "household_income14",
    # Age 14: science and mathematics orientation
    "science_liking14", "science_competence14", "science_favourite14",
    "math_liking14",
    # Age 15: science course route
    "science_course_route15",
    # Age 16: continuation, HE likelihood, intended HE subject route
    "w3_continue_fte16", "he_likelihood16",
    "w3_route_stage16", "w3_intended_domain16", "w3_intention_route16",
    # Age 20: HE pathway
    "he_participation20", "he_domain20", "he_pathway20",
    # Age 25: degree pathway
    "degree_attained25", "degree_domain25", "degree_pathway25",
    # Age 32: outcome
    "employed32", "sret_among_employed32", "sret_professional_32"
  ),
  role = c(
    rep("Analysis variable", 10),
    "Analysis variable and routed: first part", "Analysis variable",
    "Routed: second part", "Routed: third part", "Analysis variable (calculated from its parts)",
    "Routed: first part", "Routed: second part", "Analysis variable (calculated from its parts)",
    "Routed: first part", "Routed: second part", "Analysis variable (calculated from its parts)",
    "Outcome: first part", "Outcome: second part (employed-only outcome)",
    "Outcome (calculated from its parts)"
  )
)

variable_table <- rbind(
  data.frame(variable = c("log_design_weight", "stratum_group"),
             role     = "Survey design (predictor only)"),
  data.frame(variable = auxiliary_list,
             role     = "Auxiliary (File 6)"),
  analysis_variables
)

# Check: every variable in every analysis model (File 5), outcome included,
# is in the imputation model (CLS user guide).
analysis_model_variables <- c()
for (m in seq_len(nrow(analysis_models))) {
  model_formula <- as.formula(paste(analysis_models$outcome[m], "~",
                                    analysis_models$predictors[m]))
  analysis_model_variables <- c(analysis_model_variables, all.vars(model_formula))
}
not_in_model <- setdiff(unique(analysis_model_variables), variable_table$variable)
if (length(not_in_model) > 0) {
  stop("Analysis-model variables missing from the imputation model: ",
       paste(not_in_model, collapse = ", "))
}

# Check: no auxiliary is also an analysis variable.
if (any(auxiliary_list %in% analysis_variables$variable)) {
  stop("An auxiliary variable is also an analysis variable.")
}


###############################################################################
# Block 4. The imputation dataset
#
# One row per person in the imputation population, one column per variable
# in Block 3. The survey design variables used at the analysis stage (PSU,
# the original stratum, weights) and the sample flags are kept in a separate
# table, in the same row order.
###############################################################################

# The design information must be complete: it is a predictor in every
# equation, and every completed dataset is analysed with it.
if (anyNA(population$design_weight) || anyNA(population$stratum) ||
    anyNA(population$psu)) {
  stop("Some people in the imputation population lack a design weight, ",
       "stratum or PSU.")
}

population$log_design_weight <- log(population$design_weight)

# The sampling stratum as an imputation predictor. A few strata have very
# few people among the later-sweep rows that some models are fitted to, too
# few for the bootstrap logistic models (Block 8e). For the predictor only,
# the smallest strata are pooled into one group. Sizes are counted among the
# smallest of those sets of rows: people in paid work at 32 with an
# observed occupation (the rows the model for SRET among the employed is
# fitted to). The pooled strata are every stratum with fewer than 20 people
# there, plus, if needed, the next smallest strata until the pooled group
# itself has at least 20. The pooling uses group sizes only, never outcome
# values. The survey analysis keeps the original strata (design_data,
# File 12).
smallest_fitted_set <- !is.na(population$sret_among_employed32)
stratum_sizes <- sort(table(population$stratum[smallest_fitted_set]))   # smallest first
n_below_minimum <- sum(stratum_sizes < minimum_category_size)
n_pooled <- 0
if (n_below_minimum > 0) {
  n_pooled <- max(n_below_minimum,
                  which(cumsum(stratum_sizes) >= minimum_category_size)[1])
}
small_strata <- names(stratum_sizes)[seq_len(n_pooled)]
population$stratum_group <- as.character(population$stratum)
population$stratum_group[population$stratum_group %in% small_strata] <- "Pooled small strata"
population$stratum_group <- factor(population$stratum_group)

# The two 0/1 variables imputed by logistic regression are stored as
# two-category factors, the form mice's logistic method expects. The overall
# outcome stays 0/1, the form it takes in the analysis. For the employed-only
# model, SRET among the employed is turned back into 0/1 after imputation
# (Block 15 and File 12).
population$sret_professional_25 <- factor(population$sret_professional_25, levels = 0:1,
                                          labels = c("Not SRET professional",
                                                     "SRET professional"))
population$sret_among_employed32 <- factor(population$sret_among_employed32, levels = 0:1,
                                           labels = c("Other occupation",
                                                      "SRET professional"))
population$sret_professional_32 <- as.numeric(population$sret_professional_32)

# Main activity at 25 (an auxiliary) has several very small categories. They
# are grouped by meaning, without reference to the outcome (CLS user guide:
# collapse sparse categories). The script stops if a category has no group.
activity25_groups <- c(
  "Employee - in paid work"                 = "Employee",
  "Self employed"                           = "Self-employed",
  "Unemployed"                              = "Unemployed",
  "Education: School/college/university"    = "Education or training",
  "Apprenticeship"                          = "Education or training",
  "On gov't scheme for employment training" = "Education or training",
  "Looking after home or family"            = "Looking after home or family",
  "Sick or disabled"                        = "Sick, disabled, voluntary work or other",
  "In unpaid/voluntary work"                = "Sick, disabled, voluntary work or other",
  "Something else"                          = "Sick, disabled, voluntary work or other"
)
not_grouped <- setdiff(levels(population$activity25), names(activity25_groups))
if (length(not_grouped) > 0) {
  stop("Main activity at 25 has categories without a group: ",
       paste(not_grouped, collapse = ", "))
}
population$activity25 <- factor(unname(activity25_groups[as.character(population$activity25)]),
                                levels = unique(activity25_groups))

imputation_data <- as.data.frame(population[, variable_table$variable])

# Drop categories that nobody in the imputation population has, so every
# category can be estimated.
imputation_data <- droplevels(imputation_data)

design_data <- as.data.frame(population[, c("NSID", "psu", "stratum", "design_weight",
                                            "w9_final_weight", "responded_w9",
                                            "outcome_observed32",
                                            "in_prediction_cohort")])

# Every column must be numeric or a factor (mice can't use text columns).
for (variable in names(imputation_data)) {
  values <- imputation_data[[variable]]
  if (!(is.numeric(values) || is.factor(values))) {
    stop(variable, " is neither numeric nor a factor.")
  }
}


###############################################################################
# Block 5. Which values are imputed (the "where" matrix)
#
# TRUE = impute this value. Imputed: item and wave non-response (File 3).
# Never imputed: observed values and structural (not applicable) values.
# Variables with no missingness state (the design variables, previous
# non-response) are complete.
###############################################################################

where <- matrix(FALSE, nrow = nrow(imputation_data), ncol = ncol(imputation_data),
                dimnames = list(NULL, names(imputation_data)))

for (variable in names(imputation_data)) {
  state_name <- paste0("state_", variable)
  if (state_name %in% names(population)) {
    where[, variable] <- population[[state_name]] %in% c("Item non-response",
                                                         "Wave non-response")
  } else {
    where[, variable] <- is.na(imputation_data[[variable]])
  }
}

# Checks: every missing value is either imputed or structural, and no
# observed value is marked for imputation.
for (variable in names(imputation_data)) {
  state_name <- paste0("state_", variable)
  if (state_name %in% names(population)) {
    structural <- population[[state_name]] == "Structural (not applicable)"
  } else {
    structural <- rep(FALSE, nrow(population))
  }
  missing_value <- is.na(imputation_data[[variable]])
  if (any(missing_value & !where[, variable] & !structural)) {
    stop(variable, " has missing values that are neither imputed nor structural.")
  }
  if (any(where[, variable] & !missing_value)) {
    stop(variable, " has observed values marked for imputation.")
  }
}


###############################################################################
# Block 6. Imputation method for each variable
#
# The method follows the variable's type (see the header). Variables with
# nothing to impute get "" (no method). The four calculated variables get a
# formula ("passive imputation"). For the pathways, the categories are
# written into the formula, so the calculated variable keeps exactly the
# categories used in the analysis.
###############################################################################

method <- setNames(rep("", ncol(imputation_data)), names(imputation_data))

for (variable in names(imputation_data)) {
  values <- imputation_data[[variable]]
  if (!any(where[, variable])) {
    method[variable] <- ""
  } else if (is.factor(values) && nlevels(values) == 2) {
    method[variable] <- "logreg.boot"
  } else if (is.factor(values)) {
    method[variable] <- "cart"      # ordered or unordered, more than two categories
  } else {
    method[variable] <- "pmm"
  }
}

# Intended HE subject route at 16, from its three parts (the same rules as
# File 2): leaving full-time education, or no university route -> no
# university or FTE route; unsure or undecided -> itself; subject area given
# -> the subject domain.
method["w3_intention_route16"] <- paste0(
  "~ I(factor(",
  "ifelse(w3_continue_fte16 %in% 'Leave or not immediately return' | ",
  "w3_route_stage16 %in% 'No university route', ",
  "'No university or full-time education route', ",
  "ifelse(w3_route_stage16 %in% 'Unsure about university or subject undecided', ",
  "'Unsure about university or subject undecided', ",
  "ifelse(w3_route_stage16 %in% 'Subject area given', ",
  "paste('Subject given:', w3_intended_domain16), NA))), ",
  "levels = ", paste(deparse(levels(imputation_data$w3_intention_route16)), collapse = ""),
  "))"
)

# HE pathway at 20: not in HE, or in HE with its subject domain.
method["he_pathway20"] <- paste0(
  "~ I(factor(ifelse(he_participation20 == 'Not in HE by age 20', ",
  "'Not in HE by age 20', paste('In HE:', he_domain20)), ",
  "levels = ", paste(deparse(levels(imputation_data$he_pathway20)), collapse = ""),
  "))"
)

# Degree pathway at 25: no degree, or a degree with its subject domain.
method["degree_pathway25"] <- paste0(
  "~ I(factor(ifelse(degree_attained25 == 'No degree', ",
  "'No degree by age 25', paste('Degree:', degree_domain25)), ",
  "levels = ", paste(deparse(levels(imputation_data$degree_pathway25)), collapse = ""),
  "))"
)

# Outcome at 32: SRET among the employed if in paid work, otherwise 0.
method["sret_professional_32"] <- paste0(
  "~ I(ifelse(employed32 == 'In paid work', ",
  "as.numeric(sret_among_employed32 == 'SRET professional'), 0))"
)

is_passive       <- startsWith(method, "~")
imputed_by_model <- method != "" & !is_passive

# The order in which mice visits the variables (Block 3's order).
visit_sequence <- names(method)[method != ""]


###############################################################################
# Block 7. Predictor matrix
#
# Row = the variable being imputed; column = a possible predictor (1 = used).
# Start with every variable predicting every other, then apply these rules.
###############################################################################

variables <- names(imputation_data)
predictor_matrix <- matrix(1L, nrow = length(variables), ncol = length(variables),
                           dimnames = list(variables, variables))
diag(predictor_matrix) <- 0L

# 7a. Only variables imputed by a model need predictors. Complete variables
# and calculated (passive) variables have none.
predictor_matrix[!imputed_by_model, ] <- 0L

# 7b. The routed parts of the pathways don't predict other variables.
# - Parts with structural (not applicable) values keep them as missing, so
#   they can't be predictors: the Wave-3 route and subject domain, the HE
#   and degree subjects, and SRET among the employed.
# - The first parts of the HE and degree pathways are exact functions of the
#   calculated pathways, so using both would duplicate information.
# - The calculated pathways and the overall outcome carry the information of
#   their parts.
# Continuation in full-time education at 16, and being in paid work at 32,
# are first parts that aren't functions of their calculated variables, so
# they remain predictors of the other variables.
predictor_matrix[, c("w3_route_stage16", "w3_intended_domain16",
                     "he_participation20", "he_domain20",
                     "degree_attained25", "degree_domain25",
                     "sret_among_employed32")] <- 0L

# 7c. A calculated variable doesn't predict its own parts: it is computed
# from them, so it would feed their own values back into them.
predictor_matrix[c("w3_continue_fte16", "w3_route_stage16", "w3_intended_domain16"),
                 "w3_intention_route16"] <- 0L
predictor_matrix[c("he_participation20", "he_domain20"), "he_pathway20"] <- 0L
predictor_matrix[c("degree_attained25", "degree_domain25"), "degree_pathway25"] <- 0L
predictor_matrix[c("employed32", "sret_among_employed32"), "sret_professional_32"] <- 0L

# 7d. A routed part is fitted only to the people it applies to, for whom
# the earlier answer is (nearly) constant, so that answer isn't a predictor
# there:
#   SRET among the employed: fitted to people in paid work;
#   the Wave-3 route and subject domain: fitted to people staying on in
#   full-time education.
predictor_matrix["sret_among_employed32", "employed32"] <- 0L
predictor_matrix[c("w3_route_stage16", "w3_intended_domain16"), "w3_continue_fte16"] <- 0L

# 7e. Survey design: the log design weight and the stratum predict every
# imputed variable (the departure from Silverwood et al. described in the
# header).
predictor_matrix[imputed_by_model, c("log_design_weight", "stratum_group")] <- 1L


###############################################################################
# Block 8. Checks of the specification
#
# Each check stops the script if it fails. No data values are printed.
###############################################################################

# 8a. Each calculated variable's formula reproduces its observed values.
for (variable in names(method)[is_passive]) {
  rebuilt  <- model.frame(as.formula(method[[variable]]), imputation_data,
                          na.action = na.pass)[[1]]
  observed <- !is.na(imputation_data[[variable]])
  if (anyNA(rebuilt[observed]) ||
      any(as.character(rebuilt[observed]) !=
          as.character(imputation_data[[variable]][observed]))) {
    stop("The formula for ", variable, " doesn't reproduce its observed values.")
  }
}

# 8b. Every variable imputed by a model has the design predictors.
if (any(predictor_matrix[imputed_by_model, c("log_design_weight", "stratum_group")] != 1)) {
  stop("A variable is imputed without the design predictors.")
}

# 8c. No variable that keeps structural missing values is a predictor.
if (any(predictor_matrix[, c("w3_route_stage16", "w3_intended_domain16", "he_domain20",
                             "degree_domain25", "sret_among_employed32")] != 0)) {
  stop("A routed part with structural missing values is used as a predictor.")
}

# 8d. Every variable to be imputed has a method, and every method has
# something to impute.
if (any(colSums(where) > 0 & method == "") || any(colSums(where) == 0 & method != "")) {
  stop("The methods don't match the values to be imputed.")
}

# 8e. Bootstrap logistic models: every category of every categorical (or
# 0/1) predictor has at least 20 people among the rows the model is fitted
# to (the rows where the imputed variable is observed). Only observed
# predictor values are counted; imputation can only add to these counts, so
# this is a safe lower bound. Every model and predictor is checked, and
# each category below the minimum is listed by name in
# 07_category_check.csv (no counts). The script then stops if any was
# found, so that all sparse categories can be dealt with at once.
# First the (model, predictor) pairs to check, then one loop over them.
check_pairs <- data.frame()
for (variable in names(method)[method == "logreg.boot"]) {
  predictors_v <- names(which(predictor_matrix[variable, ] == 1))
  check_pairs <- rbind(check_pairs, data.frame(model_for = rep(variable, length(predictors_v)),
                                               predictor = predictors_v))
}
category_check <- data.frame(model_for = character(), predictor = character(),
                             category = character(), fewer_than = numeric())
for (j in seq_len(nrow(check_pairs))) {
  fitted_rows <- !is.na(imputation_data[[check_pairs$model_for[j]]])
  values <- imputation_data[[check_pairs$predictor[j]]][fitted_rows]
  is_categorical <- is.factor(values) || all(values %in% c(0, 1, NA))
  counts <- table(values)
  if (is_categorical && any(counts < minimum_category_size)) {
    category_check <- rbind(category_check, data.frame(
      model_for = check_pairs$model_for[j], predictor = check_pairs$predictor[j],
      category = names(counts)[counts < minimum_category_size],
      fewer_than = minimum_category_size
    ))
  }
}

# Disclosure check: names only. A listed category has fewer than 20 people
# among the model's rows; the exact count is not shown.
category_check$disclosure_checked <- rep(TRUE, nrow(category_check))
write.csv(category_check, file.path(review_dir, "07_category_check.csv"),
          row.names = FALSE)

if (nrow(category_check) > 0) {
  stop(nrow(category_check), " predictor categories have fewer than ",
       minimum_category_size, " people among the rows a bootstrap logistic ",
       "model is fitted to; see 07_category_check.csv. Collapse them ",
       "(Block 4) before the test run.")
}


###############################################################################
# Block 9. Record the specification (metadata)
###############################################################################

n_imputed <- colSums(where)   # values imputed per variable (per imputation)

variable_table$type <- NA_character_
for (i in seq_len(nrow(variable_table))) {
  values <- imputation_data[[variable_table$variable[i]]]
  variable_table$type[i] <- if (is.ordered(values)) {
    paste("ordered,", nlevels(values), "categories")
  } else if (is.factor(values)) {
    paste(nlevels(values), "categories")
  } else if (all(values %in% c(0, 1, NA))) {
    "0/1"
  } else {
    "numeric"
  }
}
variable_table$method        <- unname(method[variable_table$variable])
variable_table$visit_order   <- match(variable_table$variable, visit_sequence)
variable_table$n_predictors  <- unname(rowSums(predictor_matrix)[variable_table$variable])
variable_table$disclosure_checked <- TRUE

write.csv(variable_table, file.path(review_dir, "07_imputation_model.csv"),
          row.names = FALSE, na = "")
write.csv(data.frame(imputed_variable = rownames(predictor_matrix),
                     predictor_matrix, check.names = FALSE),
          file.path(review_dir, "07_predictor_matrix.csv"), row.names = FALSE)


###############################################################################
# Block 10. Test run
#
# 5 imputations of 30 iterations, with a fixed seed. mice prints its
# progress (iteration, imputation and variable names only). Warnings are
# collected rather than printed, and written to 07_warnings.csv.
###############################################################################

# mice finds two of its internal objects, "state" and "loggedEvents", by
# name. If the R workspace holds an object with either name (e.g. from other
# code run in the same session), mice picks it up instead and fails with
# "$ operator is invalid for atomic vectors" the first time it logs an
# event. Stop now with a clear message instead.
for (name in c("state", "loggedEvents")) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) {
    stop("An object called '", name, "' is in the R workspace, and mice would ",
         "mistake it for its own. Remove it with rm(", name, ") or restart R, ",
         "then run this block again.")
  }
}

cat("\nTest run:", test_m, "imputations of", iterations, "iterations.\n")
start_time <- Sys.time()

# The imputation, keeping its warnings: withCallingHandlers() runs mice();
# each time a warning occurs, the small function after "warning =" adds its
# message to test_warnings (the "<<-" writes to that list outside the
# function) and lets the imputation carry on.
test_warnings <- character()
test_imputation <- withCallingHandlers(
  mice(imputation_data,
       m               = test_m,
       maxit           = iterations,
       method          = method,
       predictorMatrix = predictor_matrix,
       where           = where,
       visitSequence   = visit_sequence,
       seed            = test_seed,
       printFlag       = TRUE),
  warning = function(w) {
    # Keep the warning's message and carry on.
    test_warnings <<- c(test_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

test_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

# The completed datasets, used by the checks and diagnostics below.
completed <- list()
for (i in seq_len(test_m)) {
  completed[[i]] <- complete(test_imputation, i)
}


###############################################################################
# Block 11. Routing and consistency checks on the completed data
#
# For each completed dataset:
#   a. observed values are unchanged;
#   b. only structural values are still missing;
#   c. each calculated variable equals its formula applied to its imputed
#      parts, for everyone;
#   d. the outcome is 0 or 1 for everyone.
# The pass/fail table is written first; the script then stops if any check
# failed.
###############################################################################

routing_checks <- data.frame()

for (i in seq_len(test_m)) {
  data_i <- completed[[i]]

  # a. Observed values unchanged.
  unchanged <- TRUE
  for (variable in variables) {
    observed <- !is.na(imputation_data[[variable]])
    if (any(as.character(data_i[[variable]][observed]) !=
            as.character(imputation_data[[variable]][observed]))) {
      unchanged <- FALSE
    }
  }

  # b. Only structural values still missing.
  only_structural_missing <- TRUE
  for (variable in variables) {
    still_missing  <- is.na(data_i[[variable]])
    expect_missing <- is.na(imputation_data[[variable]]) & !where[, variable]
    if (any(still_missing != expect_missing)) only_structural_missing <- FALSE
  }

  # c. Calculated variables equal their formulas.
  calculated_consistent <- TRUE
  for (variable in names(method)[is_passive]) {
    rebuilt <- model.frame(as.formula(method[[variable]]), data_i,
                           na.action = na.pass)[[1]]
    if (anyNA(rebuilt) ||
        any(as.character(rebuilt) != as.character(data_i[[variable]]))) {
      calculated_consistent <- FALSE
    }
  }

  # d. Outcome is 0 or 1.
  outcome_binary <- all(data_i$sret_professional_32 %in% c(0, 1))

  routing_checks <- rbind(routing_checks, data.frame(
    imputation = i,
    check = c("Observed values unchanged",
              "Only structural (not applicable) values still missing",
              "Pathways and outcome equal their parts",
              "Outcome is 0 or 1 for everyone"),
    passed = c(unchanged, only_structural_missing, calculated_consistent,
               outcome_binary)
  ))
}

# Disclosure check: pass/fail results only; no counts or data values.
routing_checks$disclosure_checked <- TRUE
write.csv(routing_checks, file.path(review_dir, "07_routing_checks.csv"),
          row.names = FALSE)

if (!all(routing_checks$passed)) {
  stop("A routing or consistency check failed; see 07_routing_checks.csv.")
}


###############################################################################
# Block 12. Trace plots
#
# For every imputed variable: the mean and SD of its imputed values at each
# iteration, one line per imputation. The lines should mix and show no
# trend by the last iterations. For categorical variables, mice records the
# mean and SD of the category numbers (1, 2, ...), which is enough to see a
# trend.
# All variables go to a private file. The review copy includes only
# variables with at least 10 imputed values, so every plotted mean and SD
# rests on at least 10 values.
###############################################################################

traced_variables <- visit_sequence
review_traced    <- traced_variables[n_imputed[traced_variables] >= disclosure_threshold]

for (copy in c("private", "review")) {
  if (copy == "private") {
    pdf(file.path(output_dir, "private_07_trace_plots_all.pdf"), width = 8.5, height = 11)
    plot_variables <- traced_variables
  } else {
    pdf(file.path(review_dir, "07_trace_plots.pdf"), width = 8.5, height = 11)
    plot_variables <- review_traced
  }
  par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
  for (variable in plot_variables) {
    chain_means <- test_imputation$chainMean[variable, , ]        # iterations x imputations
    chain_sds   <- sqrt(test_imputation$chainVar[variable, , ])
    matplot(chain_means, type = "l", lty = 1, xlab = "Iteration",
            ylab = "Mean of imputed values", main = variable, cex.main = 0.9)
    matplot(chain_sds, type = "l", lty = 1, xlab = "Iteration",
            ylab = "SD of imputed values", main = variable, cex.main = 0.9)
  }
  dev.off()
}


###############################################################################
# Block 12b. Convergence table
#
# A numeric check alongside the trace plots, for the mean and the SD of each
# variable's imputed values, over the second half of the iterations:
#   split R-hat (Vehtari et al. 2021): each imputation's sequence is cut
#     into two halves, and the variation between all the half-sequences is
#     compared with the variation within them. Values close to 1 mean they
#     agree. A trend makes the first and second halves differ, so split
#     R-hat also rises when all imputations drift together, which the
#     classic R-hat can miss.
#   trend: the average change per iteration, from a regression of the
#     statistic on iteration with a separate level for each imputation, and
#     its p-value. Successive iterations are correlated, so the p-value
#     overstates the evidence. It is therefore paired with an effect size:
#   drift_vs_spread: the total change over the second half (trend x number
#     of steps), divided by the typical spread between imputations at an
#     iteration. Above 1, the drift is larger than the differences between
#     the imputations themselves.
# A variable is flagged for a closer look if split R-hat is above 1.1, or if
# the trend has p < 0.01 and drift_vs_spread is above 1. Flags are a guide,
# not a rule: Oberman, van Buuren and
# Vink (2021) found valid estimates after fewer iterations than such
# diagnostics suggest. As for the trace plots, only variables with at least
# 10 imputed values are included.
###############################################################################

second_half <- (floor(iterations / 2) + 1):iterations        # e.g. 16-30
half_length <- floor(length(second_half) / 2)
early_part  <- head(second_half, half_length)                 # e.g. 11-15
late_part   <- tail(second_half, half_length)                 # e.g. 16-20

convergence <- data.frame(variable = character(), statistic = character(),
                          split_rhat = numeric(), trend_per_iteration = numeric(),
                          trend_p = numeric(), drift_vs_spread = numeric(),
                          flag = logical())
# One row per variable x statistic (none if the halves are too short).
trace_table <- expand.grid(statistic = c("mean", "sd"), variable = review_traced,
                           stringsAsFactors = FALSE)
if (half_length < 2) trace_table <- trace_table[0, ]

for (j in seq_len(nrow(trace_table))) {
  variable  <- trace_table$variable[j]
  statistic <- trace_table$statistic[j]
  values <- test_imputation$chainMean[variable, , ]                # iterations x imputations
  if (statistic == "sd") values <- sqrt(test_imputation$chainVar[variable, , ])

  # Split R-hat: 2 x test_m half-sequences of equal length.
  halves  <- cbind(values[early_part, ], values[late_part, ])
  n       <- nrow(halves)
  within  <- mean(apply(halves, 2, var))                           # variance within each half-sequence
  between <- n * var(colMeans(halves))                             # variance between half-sequence means
  split_rhat <- NA
  if (within > 0) split_rhat <- sqrt(((n - 1) / n * within + between / n) / within)

  # Trend over the second half of the iterations.
  trend_data <- data.frame(
    value      = as.vector(values[second_half, ]),
    iteration  = rep(second_half, times = test_m),
    imputation = factor(rep(seq_len(test_m), each = length(second_half)))
  )
  trend <- summary(lm(value ~ iteration + imputation, data = trend_data))$coefficients
  trend_slope <- trend["iteration", "Estimate"]
  trend_p     <- trend["iteration", "Pr(>|t|)"]

  # Drift over the second half, relative to the spread between imputations
  # (the SD across imputations, averaged over iterations).
  spread <- mean(apply(values[second_half, ], 1, sd))
  drift_vs_spread <- NA
  if (spread > 0) drift_vs_spread <- abs(trend_slope) * (length(second_half) - 1) / spread

  convergence <- rbind(convergence, data.frame(
    variable            = variable,
    statistic           = statistic,
    split_rhat          = round(split_rhat, 3),
    trend_per_iteration = signif(trend_slope, 3),
    trend_p             = signif(trend_p, 3),
    drift_vs_spread     = round(drift_vs_spread, 2),
    flag                = isTRUE(split_rhat > 1.1) ||
                          (isTRUE(trend_p < 0.01) && isTRUE(drift_vs_spread > 1))
  ))
}

# Disclosure check: summaries of the means and SDs of imputed values, only
# for variables with at least 10 imputed values; no counts.
convergence$disclosure_checked <- rep(TRUE, nrow(convergence))
write.csv(convergence, file.path(review_dir, "07_convergence.csv"), row.names = FALSE)


###############################################################################
# Block 13. Observed versus imputed distributions
#
# For every imputed variable, the distribution of its observed values and of
# its imputed values (all test imputations together):
#   categorical (including 0/1): percentage in each category;
#   numeric: mean, SD and quartiles.
# Large differences aren't errors in themselves (under missing at random the
# imputed values can differ from the observed ones), but they should be
# explainable. For a routed second part, the imputed values include people
# whose imputed earlier answer doesn't lead to it; the calculated pathway
# shows the values that are used.
###############################################################################

comparison <- data.frame()

for (variable in traced_variables) {
  observed_values <- imputation_data[[variable]][!is.na(imputation_data[[variable]])]
  imputed_values  <- c()
  for (i in seq_len(test_m)) {
    imputed_values <- c(imputed_values,
                        as.character(completed[[i]][[variable]][where[, variable]]))
  }
  is_categorical <- is.factor(imputation_data[[variable]]) ||
    all(imputation_data[[variable]] %in% c(0, 1, NA))

  if (is_categorical) {
    categories <- if (is.factor(observed_values)) levels(observed_values) else c("0", "1")
    for (category in categories) {
      observed_n     <- sum(as.character(observed_values) == category)
      imputed_mean_n <- sum(imputed_values == category) / test_m
      comparison <- rbind(comparison, data.frame(
        variable = variable, category = category, n_imputed = n_imputed[[variable]],
        observed_n = observed_n,
        observed_percent = round(100 * observed_n / length(observed_values), 1),
        imputed_mean_n = imputed_mean_n,
        imputed_percent = round(100 * imputed_mean_n / n_imputed[[variable]], 1),
        observed_mean = NA, observed_sd = NA, observed_p25 = NA,
        observed_median = NA, observed_p75 = NA,
        imputed_mean = NA, imputed_sd = NA, imputed_p25 = NA,
        imputed_median = NA, imputed_p75 = NA
      ))
    }
  } else {
    imputed_numbers <- as.numeric(imputed_values)
    observed_q <- quantile(observed_values, c(0.25, 0.5, 0.75), names = FALSE)
    imputed_q  <- quantile(imputed_numbers, c(0.25, 0.5, 0.75), names = FALSE)
    comparison <- rbind(comparison, data.frame(
      variable = variable, category = "(numeric)", n_imputed = n_imputed[[variable]],
      observed_n = length(observed_values), observed_percent = NA,
      imputed_mean_n = NA, imputed_percent = NA,
      observed_mean = round(mean(observed_values), 2),
      observed_sd = round(sd(observed_values), 2),
      observed_p25 = observed_q[1], observed_median = observed_q[2],
      observed_p75 = observed_q[3],
      imputed_mean = round(mean(imputed_numbers), 2),
      imputed_sd = round(sd(imputed_numbers), 2),
      imputed_p25 = imputed_q[1], imputed_median = imputed_q[2],
      imputed_p75 = imputed_q[3]
    ))
  }
}

write.csv(comparison, file.path(output_dir, "private_07_observed_vs_imputed.csv"),
          row.names = FALSE)

# --- Disclosure check: observed versus imputed --------------------------------
# Categories: a row's percentages are hidden if fewer than 10 people are
# observed in the category, or fewer than 10 values (average per
# imputation) are imputed into it. Secondary: within each variable, if
# exactly one row is hidden, the row with the next smallest observed count
# is hidden too (percentages add to 100). Counts are removed from the shared
# table, so it can't be combined with other tables' counts.
# Numeric rows: observed statistics only if at least 10 are observed;
# imputed statistics only if at least 10 values are imputed.
# The number imputed per variable is hidden if it is 1-9.
checked <- comparison
is_numeric_row <- checked$category == "(numeric)"

primary <- !is_numeric_row &
  ((checked$observed_n > 0 & checked$observed_n < disclosure_threshold) |
   (checked$imputed_mean_n > 0 & checked$imputed_mean_n < disclosure_threshold))
secondary <- rep(FALSE, nrow(checked))
for (v in unique(checked$variable[!is_numeric_row])) {
  rows <- which(checked$variable == v)
  if (sum(primary[rows]) == 1) {
    candidates <- rows[!primary[rows] & checked$observed_n[rows] > 0]
    if (length(candidates) > 0) {
      secondary[candidates[which.min(checked$observed_n[candidates])]] <- TRUE
    }
  }
}
checked[primary | secondary, c("observed_percent", "imputed_percent")] <- NA
checked$fewer_than_10 <- primary

numeric_observed_small <- is_numeric_row & checked$observed_n < disclosure_threshold
numeric_imputed_small  <- is_numeric_row & checked$n_imputed < disclosure_threshold
checked[numeric_observed_small, c("observed_mean", "observed_sd", "observed_p25",
                                  "observed_median", "observed_p75")] <- NA
checked[numeric_imputed_small, c("imputed_mean", "imputed_sd", "imputed_p25",
                                 "imputed_median", "imputed_p75")] <- NA

small_n_imputed <- checked$n_imputed > 0 & checked$n_imputed < disclosure_threshold
checked$n_imputed[small_n_imputed] <- NA

checked <- checked[, setdiff(names(checked), c("observed_n", "imputed_mean_n"))]
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "07_observed_vs_imputed.csv"),
          row.names = FALSE, na = "")

disclosure_log <- data.frame(
  table = c("07_observed_vs_imputed", "07_trace_plots"),
  primary_suppressed = c(sum(primary) + sum(numeric_observed_small) +
                           sum(numeric_imputed_small) +
                           length(unique(checked$variable[small_n_imputed])),
                         length(traced_variables) - length(review_traced)),
  secondary_suppressed = c(sum(secondary), 0L)
)
# -----------------------------------------------------------------------------


###############################################################################
# Block 14. Logged events and warnings
#
# mice logs every change it makes to the model while running, e.g. a
# predictor category removed because it is constant (or collinear) among the
# rows an equation is fitted to. Both tables hold variable names, method
# names and messages only.
###############################################################################

logged_events <- test_imputation$loggedEvents
if (is.null(logged_events)) {
  logged_events <- data.frame(it = integer(), im = integer(), dep = character(),
                              meth = character(), out = character())
}
logged_events$disclosure_checked <- rep(TRUE, nrow(logged_events))
write.csv(logged_events, file.path(review_dir, "07_logged_events.csv"),
          row.names = FALSE)

warning_table <- data.frame(message = unique(test_warnings))
warning_table$times <- as.integer(table(test_warnings)[warning_table$message])
warning_table$disclosure_checked <- rep(TRUE, nrow(warning_table))
write.csv(warning_table, file.path(review_dir, "07_warnings.csv"), row.names = FALSE)


###############################################################################
# Block 15. Preliminary fraction of missing information (FMI)
#
# Each analysis model from File 5 (M0-M4 and the employed-only M0) is fitted
# to each completed test dataset, in the primary analysis sample: the full
# cohort with the initial design weight, PSU and stratum. For each
# coefficient, Rubin's rules give:
#   within-imputation variance   W = mean of the squared standard errors
#   between-imputation variance  B = variance of the estimates
#   total variance               T = W + (1 + 1/m) B
#   FMI (approximate)            (1 + 1/m) B / T
# With 5 imputations this is a rough estimate. The main run uses at least
# 50 imputations, rounded up to the next 10 above the largest FMI as a
# percentage (CLS user guide). File 12 checks the rule again.
###############################################################################

fmi_models <- analysis_models[analysis_models$model != "M0 baseline (own complete cases)", ]
fmi_table  <- data.frame()

for (m in seq_len(nrow(fmi_models))) {
  model_outcome <- fmi_models$outcome[m]
  model_formula <- as.formula(paste(model_outcome, "~", fmi_models$predictors[m]))
  estimates <- NULL
  variances <- NULL
  smallest_outcome_group <- Inf

  for (i in seq_len(test_m)) {
    analysis_data <- cbind(design_data, completed[[i]])
    # As in File 5: HE likelihood enters as unordered categories, and SRET
    # among the employed is a 0/1 outcome.
    analysis_data$he_likelihood16 <- factor(analysis_data$he_likelihood16, ordered = FALSE)
    analysis_data$sret_among_employed32 <-
      as.numeric(analysis_data$sret_among_employed32 == "SRET professional")
    design <- svydesign(ids = ~psu, strata = ~stratum, weights = ~design_weight,
                        nest = TRUE, data = analysis_data)
    if (model_outcome == "sret_among_employed32") {
      design <- subset(design, employed32 == "In paid work")
    }
    fit <- svyglm(model_formula, design = design, family = quasibinomial(),
                  control = glm.control(maxit = 100))
    estimates <- cbind(estimates, coef(fit))
    variances <- cbind(variances, diag(vcov(fit)))

    y <- design$variables[[model_outcome]]
    smallest_outcome_group <- min(smallest_outcome_group, sum(y == 1), sum(y == 0))
  }

  within  <- rowMeans(variances)
  between <- apply(estimates, 1, var)
  total   <- within + (1 + 1 / test_m) * between
  fmi     <- (1 + 1 / test_m) * between / total

  fmi_table <- rbind(fmi_table, data.frame(
    model                 = fmi_models$model[m],
    coefficients          = length(fmi),
    median_fmi            = round(median(fmi, na.rm = TRUE), 3),
    largest_fmi           = round(max(fmi, na.rm = TRUE), 3),
    term_with_largest_fmi = names(fmi)[which.max(fmi)],
    smallest_outcome_group = smallest_outcome_group
  ))
}

largest_fmi  <- max(fmi_table$largest_fmi)
production_m <- max(minimum_m, 10 * ceiling(100 * largest_fmi / 10))

write.csv(fmi_table, file.path(output_dir, "private_07_fmi.csv"), row.names = FALSE)

# --- Disclosure check: FMI ---------------------------------------------------
# A model's row is shown only if it has at least 10 events and 10 non-events
# in every completed dataset. The event counts are removed from the shared
# table.
checked <- fmi_table
model_too_small <- checked$smallest_outcome_group < disclosure_threshold
checked[model_too_small, c("coefficients", "median_fmi", "largest_fmi",
                           "term_with_largest_fmi")] <- NA
checked <- checked[, setdiff(names(checked), "smallest_outcome_group")]
checked$production_m <- production_m
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "07_fmi.csv"), row.names = FALSE, na = "")
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "07_fmi", primary_suppressed = sum(model_too_small),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------

disclosure_log <- rbind(disclosure_log, data.frame(
  table = c("07_imputation_model", "07_predictor_matrix", "07_category_check",
            "07_routing_checks", "07_convergence", "07_logged_events",
            "07_warnings"),
  primary_suppressed = 0L, secondary_suppressed = 0L))
disclosure_log$threshold <- disclosure_threshold
write.csv(disclosure_log, file.path(review_dir, "07_disclosure_log.csv"),
          row.names = FALSE)


###############################################################################
# Block 16. Save
#
# The specification holds everything the imputation runs need (Files 8a,
# 9a, 11a and 13a): the data, the design table in the same row order, the
# methods, predictor matrix, where matrix, visit sequence, iterations and
# number of imputations.
###############################################################################

imputation_spec <- list(
  data             = imputation_data,
  design_data      = design_data,
  variables        = variable_table,
  method           = method,
  predictor_matrix = predictor_matrix,
  where            = where,
  visit_sequence   = visit_sequence,
  iterations       = iterations,
  production_m     = production_m
)
saveRDS(imputation_spec, file.path(data_dir, "07_imputation_spec.rds"))
saveRDS(test_imputation, file.path(data_dir, "07_test_imputation.rds"))


###############################################################################
# Block 17. Summary
#
# Counts of variables, events and checks only. No data values.
###############################################################################

cat("\nFile 7 finished.\n")
cat("  Variables in the imputation model:    ", ncol(imputation_data), "\n")
cat("    imputed by a model:                 ", sum(imputed_by_model), "\n")
cat("    calculated from their parts:        ", sum(is_passive), "\n")
cat("    complete (predictors only):         ", sum(method == ""), "\n")
cat("  Auxiliaries (File 6):                 ", length(auxiliary_list), "\n")
cat("  Stratum predictor categories:         ", nlevels(imputation_data$stratum_group),
    "(", length(small_strata), "small strata pooled )\n")
cat("  Methods: pmm", sum(method == "pmm"), "| logreg.boot", sum(method == "logreg.boot"),
    "| cart", sum(method == "cart"), "\n")
cat("  Bootstrap logistic category check:     passed (at least",
    minimum_category_size, "per category)\n")
cat("  Test run minutes (all imputations):   ", round(test_minutes, 1), "\n")
cat("    about minutes per imputation:       ", round(test_minutes / test_m, 1), "\n")
cat("  Logged events:                        ", nrow(logged_events), "\n")
cat("  Distinct warnings:                    ", nrow(warning_table), "\n")
cat("  Routing and consistency checks:        all passed\n")
cat("  Convergence flags (07_convergence.csv):", sum(convergence$flag), "of",
    nrow(convergence), "variable x statistic rows\n")
cat("  Largest preliminary FMI:              ", largest_fmi, "\n")
cat("  Imputations for the main run:         ", production_m, "\n")
cat("  Review before File 8a: 07_trace_plots.pdf, 07_observed_vs_imputed.csv,\n")
cat("    07_logged_events.csv and 07_warnings.csv\n")
cat("  Saved:", file.path(data_dir, "07_imputation_spec.rds"), "\n")
