###############################################################################
# 16. Pipeline technical check
#
# Purpose
#   Run the per-fold prediction pipeline of File 17a on one outer split
#   (repeat 1, outer fold 1) at every landmark, and check that it works
#   before the full cross-validation. It checks, and never looks at how well
#   the models predict:
#     - the imputation-model category check, for every outer training set and
#       landmark (a bootstrap logistic model needs every predictor category
#       among its rows; where one doesn't, that variable is imputed with a
#       classification tree in that training set);
#     - leakage: test outcomes are masked, the imputation models are fitted
#       to training rows only, and imputations of the training rows don't
#       depend on the test rows at all;
#     - the imputation completes, respects the routing, and converges (trace
#       plots and a convergence table);
#     - warnings and logged events are captured;
#     - the fitted models are valid and their predictions are well formed;
#     - run time, to plan File 17a.
#   No outer-test performance measure is computed. Pilots check technical
#   behaviour only, so that no result informs the design.
#
# The per-fold pipeline (as in File 17a; one outer split x landmark)
#   1. Training rows are the outer training set; test rows are the outer
#      test fold.
#   2. Auxiliary screen (CLS user guide): each CLS predictor of non-response
#      known by the landmark is tested for association with the landmark's
#      most-missing feature, among training rows only (F test, chi-square test
#      or correlation test, by type), and kept if p < 0.05.
#   3. Multiple imputation of the landmark's features (Wahl et al. 2016;
#      Mertens et al. 2020, section 3.1):
#        - variables: the features, the routed parts of the combined features
#          (which are recalculated from their parts at every iteration, as in
#          File 7), the kept auxiliaries and, in the primary analysis, the
#          outcome;
#        - the outcome is used for training rows and set to missing for test
#          rows (CLS user guide: the outcome belongs in the imputation
#          model); the imputed test outcomes are discarded;
#        - the imputation models are fitted to training rows only (mice's
#          "ignore"); test rows are imputed from those models;
#        - methods, routing and predictor rules as in File 7, except that
#          there are no design variables (see "Survey weights" below):
#          predictive mean matching (numeric), bootstrap logistic
#          regression (two categories), classification trees (more
#          categories); structural
#          (not applicable) values are never imputed. Where a two-category
#          variable's predictors have a category with fewer than 20 people
#          among that training set's rows, a classification tree is used for
#          it instead (checked on training rows only);
#        - in the outcome-omitted sensitivity analysis, the outcome is left
#          out of every imputation equation.
#   4. For each imputed dataset: rare categories (fewer than 20 training
#      people) are pooled, learned from training rows only; the elastic net
#      is tuned by grouped inner cross-validation on the training rows (whole
#      PSUs; File 14) over alpha and lambda (binomial deviance); the elastic
#      net, an unpenalized logistic regression and a null (prevalence) model
#      are fitted to the training rows; and test-row probabilities are
#      predicted.
#   5. Predictions are averaged over the imputed datasets (Mertens et al.
#      2020, Approach 1).
#   Models are fitted unweighted; performance (File 18) is weighted by the
#   Wave-9 final weight. Main activity at 25 is grouped into the same six
#   groups as in File 7 before imputation.
#
# Survey weights and the imputation model (a deliberate departure)
#   Silverwood et al. (2024) weighted their imputation models by the initial
#   design weights, and the imputation analysis (File 7) includes the log
#   design weight and the sampling stratum as predictors in every equation.
#   Neither is done here: these imputation models are unweighted and contain
#   no design variables. An imputation model should match the analysis that
#   uses the completed data. The imputation analysis estimates
#   survey-weighted population quantities, so its imputation carries the
#   survey design. The prediction models here are fitted unweighted, and an
#   imputation model without design variables matches them. The prediction
#   analysis refers to the Wave-9 target population through the Wave-9 final
#   weight, applied when performance is estimated (File 18), not through the
#   imputation. The design oversampled deprived schools and minority ethnic
#   pupils; ethnic group and family background are features at every
#   landmark.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   The review tables hold names, pass/fail checks, p-values, counts of
#   model coefficients (not of people), run times and convergence summaries.
#   Trace plots and the convergence table show only variables with at least
#   10 imputed values. Category-check failures are listed by name, without
#   counts. No prediction or performance measure is written.
#
# Reads
#   data_dir/14_prediction_design.rds
#
# Writes
#   review_dir/16_checks.csv             pass/fail technical checks
#   review_dir/16_category_check.csv     sparse predictor categories, every
#                                        split x landmark (names only)
#   review_dir/16_auxiliary_screen.csv   screen results for the check split
#   review_dir/16_convergence.csv        convergence table
#   review_dir/16_trace_plots.pdf        trace plots
#   review_dir/16_model_checks.csv       tuning and fitting summaries
#   review_dir/16_timing.csv             run times and projected File 17a time
#   review_dir/16_warnings.csv           warnings and logged events
#
# References
#   Mertens BJA, Banzato E, de Wreede LC (2020). Biometrical Journal 62(3):
#   724-741.
#   Wahl S, et al. (2016). BMC Medical Research Methodology 16: 144.
#   Silverwood RJ, et al. (2024). Longitudinal and Life Course Studies 15(2):
#   227-250.
#   Silverwood R, et al. (2026). Handling missing data in the CLS cohort
#   studies: User guide (Version 2). UCL Centre for Longitudinal Studies.
###############################################################################

suppressPackageStartupMessages(library(mice))
suppressPackageStartupMessages(library(glmnet))


###############################################################################
# Block 1. Paths and settings
###############################################################################

paths      <- readRDS("project_paths.rds")
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs

disclosure_threshold  <- 10
check_repeat          <- 1       # the outer split used for the full check
check_fold            <- 1
imputations           <- 5       # imputations per outer training set
iterations            <- 30      # iterations per imputation (as in File 17a)
alpha_grid            <- c(0, 0.25, 0.5, 0.75, 1)   # elastic net mixing
minimum_category_size <- 20      # bootstrap-logistic check and rare-level pooling
screen_p_cutoff       <- 0.05    # auxiliary screen
base_seed             <- 1600    # the check's imputation seeds (File 17a has its own)


###############################################################################
# Block 2. Read the prediction design (File 14)
###############################################################################

design             <- readRDS(file.path(data_dir, "14_prediction_design.rds"))
cohort             <- design$cohort
landmark_variables <- design$landmark_variables
settings           <- design$settings
outcome            <- settings$outcome
landmark_order     <- settings$landmark_order


###############################################################################
# Block 3. Prepare the variables (as File 7)
#
#   - HE likelihood enters as unordered categories;
#   - main activity at 25 is grouped into six groups by meaning (as File 7);
#   - two-category 0/1 variables are stored as factors for the bootstrap
#     logistic imputation: SRET at 25 and, for imputation, the outcome.
###############################################################################

cohort$he_likelihood16 <- factor(cohort$he_likelihood16, ordered = FALSE)

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
not_grouped <- setdiff(levels(cohort$activity25), names(activity25_groups))
if (length(not_grouped) > 0) {
  stop("Main activity at 25 has categories without a group: ",
       paste(not_grouped, collapse = ", "))
}
cohort$activity25 <- factor(unname(activity25_groups[as.character(cohort$activity25)]),
                            levels = unique(activity25_groups))

cohort$sret_professional_25 <- factor(cohort$sret_professional_25, levels = 0:1,
                                      labels = c("Not SRET professional", "SRET professional"))
cohort$outcome_for_imputation <- factor(cohort[[outcome]], levels = 0:1,
                                        labels = c("Not SRET at 32", "SRET at 32"))

routed_parts <- c("w3_route_stage16", "w3_intended_domain16", "he_participation20",
                  "he_domain20", "degree_attained25", "degree_domain25")


###############################################################################
# Block 4. Imputation specification for every outer split x landmark
#
# For each repeat, outer fold and landmark (primary analysis), cheaply,
# without imputing:
#   a. training and test rows;
#   b. the auxiliary screen on training rows;
#   c. the variables, where matrix, methods and predictor matrix;
#   d. the category check of every bootstrap logistic model, among its
#      training rows (observed predictor values only; a safe lower bound).
# The specifications of the check split are kept for Block 5.
###############################################################################

category_check  <- data.frame()
screen_results  <- data.frame()
check_specs     <- list()

# One row per repeat x outer fold x landmark (landmark changing fastest).
split_table <- expand.grid(landmark = landmark_order, outer_fold = seq_len(settings$n_outer_folds),
                           repeat_number = seq_len(settings$n_repeats), stringsAsFactors = FALSE)

for (k in seq_len(nrow(split_table))) {
  r        <- split_table$repeat_number[k]
  f        <- split_table$outer_fold[k]
  landmark <- split_table$landmark[k]
  is_check_split <- r == check_repeat && f == check_fold

  # a. Training and test rows, and the landmark's variables.
  training <- cohort[[paste0("outer_fold_r", r)]] != f
  test     <- !training
  rows_l   <- landmark_variables[landmark_variables$landmark == landmark, ]
  features <- rows_l$variable[rows_l$role == "Feature"]
  part_of  <- rows_l[startsWith(rows_l$role, "Part of"), ]
  combined <- unique(sub("^Part of ", "", part_of$role))
  auxiliary_candidates <- rows_l$variable[rows_l$role == "Auxiliary"]

  # b. Auxiliary screen: the feature with the most missing values among
  # training rows, tested against each candidate on training rows. The test
  # depends on the types: correlation (both numeric), F test (one numeric),
  # chi-square test (both categorical).
  missing_counts <- c()
  for (v in features) missing_counts[v] <- sum(is.na(cohort[[v]][training]))
  screen_target  <- names(which.max(missing_counts))
  y_values <- cohort[[screen_target]][training]
  kept_auxiliaries <- c()
  for (candidate in auxiliary_candidates) {
    x_values <- cohort[[candidate]][training]
    both <- !is.na(y_values) & !is.na(x_values)
    y <- y_values[both]
    x <- x_values[both]
    if (is.factor(y)) y <- droplevels(y)
    if (is.factor(x)) x <- droplevels(x)
    can_test <- length(unique(y)) > 1 && length(unique(x)) > 1
    p_value <- NA
    test_used <- NA
    if (can_test && is.numeric(y) && is.numeric(x)) {
      p_value <- cor.test(y, x)$p.value
      test_used <- "correlation"
    }
    if (can_test && is.numeric(y) && !is.numeric(x)) {
      p_value <- anova(lm(y ~ x))[["Pr(>F)"]][1]
      test_used <- "F test"
    }
    if (can_test && !is.numeric(y) && is.numeric(x)) {
      p_value <- anova(lm(x ~ y))[["Pr(>F)"]][1]
      test_used <- "F test"
    }
    if (can_test && !is.numeric(y) && !is.numeric(x)) {
      p_value <- suppressWarnings(chisq.test(table(y, x))$p.value)
      test_used <- "chi-square"
    }
    keep <- isTRUE(p_value < screen_p_cutoff)
    if (keep) kept_auxiliaries <- c(kept_auxiliaries, candidate)
    if (is_check_split) {
      screen_results <- rbind(screen_results, data.frame(
        landmark = landmark, screened_against = screen_target,
        auxiliary = candidate, test = test_used, p_value = signif(p_value, 3),
        kept = keep))
    }
  }

  # c. The variables, in visiting order: auxiliaries, then the features
  # (each combined feature after its parts), then the outcome.
  imputation_variables <- kept_auxiliaries
  for (feature in features) {
    if (feature %in% combined) {
      parts <- part_of$variable[part_of$role == paste("Part of", feature)]
      imputation_variables <- c(imputation_variables, setdiff(parts, imputation_variables))
    }
    imputation_variables <- c(imputation_variables, setdiff(feature, imputation_variables))
  }
  imputation_variables <- c(imputation_variables, "outcome_for_imputation")

  imputation_data <- cohort[, imputation_variables]
  imputation_data$outcome_for_imputation[test] <- NA     # masked for test rows

  # Where: item and wave non-response (File 3), and the masked outcome.
  where <- matrix(FALSE, nrow = nrow(imputation_data), ncol = ncol(imputation_data),
                  dimnames = list(NULL, imputation_variables))
  for (v in imputation_variables) {
    state_name <- paste0("state_", v)
    if (state_name %in% names(cohort)) {
      where[, v] <- cohort[[state_name]] %in% c("Item non-response", "Wave non-response")
    } else {
      where[, v] <- is.na(imputation_data[[v]])
    }
  }

  # Values marked for imputation are set to missing first. This matters for
  # the activity counts at 16-19: a history that stops early gives a count
  # over the observed months only, which File 3 classes as item non-response,
  # so the count over all months is imputed.
  for (v in imputation_variables) {
    imputation_data[[v]][where[, v]] <- NA
  }

  # Checks: every missing value is either imputed or structural (not
  # applicable), as in File 7; and every variable to be imputed has observed
  # values among the training rows (an imputation model needs something to
  # learn from).
  for (v in imputation_variables) {
    state_name <- paste0("state_", v)
    structural <- rep(FALSE, nrow(cohort))
    if (state_name %in% names(cohort)) {
      structural <- cohort[[state_name]] == "Structural (not applicable)"
    }
    if (any(is.na(imputation_data[[v]]) & !where[, v] & !structural)) {
      stop(v, " has missing values that are neither imputed nor structural (",
           landmark, ", repeat ", r, ", fold ", f, ").")
    }
    if (any(where[, v]) && sum(training & !is.na(imputation_data[[v]])) == 0) {
      stop(v, " has no observed values among the training rows (", landmark,
           ", repeat ", r, ", fold ", f, "), so it can't be imputed.")
    }
  }

  # Methods by type; combined features are recalculated (as File 7).
  method <- setNames(rep("", length(imputation_variables)), imputation_variables)
  for (v in imputation_variables) {
    values <- imputation_data[[v]]
    if (any(where[, v]) && is.factor(values) && nlevels(values) == 2) method[v] <- "logreg.boot"
    if (any(where[, v]) && is.factor(values) && nlevels(values) != 2) method[v] <- "cart"
    if (any(where[, v]) && !is.factor(values))                        method[v] <- "pmm"
  }
  if ("w3_intention_route16" %in% combined) {
    method["w3_intention_route16"] <- paste0(
      "~ I(factor(",
      "ifelse(w3_continue_fte16 %in% 'Leave or not immediately return' | ",
      "w3_route_stage16 %in% 'No university route', ",
      "'No university or full-time education route', ",
      "ifelse(w3_route_stage16 %in% 'Unsure about university or subject undecided', ",
      "'Unsure about university or subject undecided', ",
      "ifelse(w3_route_stage16 %in% 'Subject area given', ",
      "paste('Subject given:', w3_intended_domain16), NA))), ",
      "levels = ", paste(deparse(levels(cohort$w3_intention_route16)), collapse = ""), "))")
  }
  if ("he_pathway20" %in% combined) {
    method["he_pathway20"] <- paste0(
      "~ I(factor(ifelse(he_participation20 == 'Not in HE by age 20', ",
      "'Not in HE by age 20', paste('In HE:', he_domain20)), ",
      "levels = ", paste(deparse(levels(cohort$he_pathway20)), collapse = ""), "))")
  }
  if ("degree_pathway25" %in% combined) {
    method["degree_pathway25"] <- paste0(
      "~ I(factor(ifelse(degree_attained25 == 'No degree', ",
      "'No degree by age 25', paste('Degree:', degree_domain25)), ",
      "levels = ", paste(deparse(levels(cohort$degree_pathway25)), collapse = ""), "))")
  }
  is_passive       <- startsWith(method, "~")
  imputed_by_model <- method != "" & !is_passive

  # Predictor matrix (File 7's rules 7a-7d, for the variables present).
  predictor_matrix <- matrix(1L, length(imputation_variables), length(imputation_variables),
                             dimnames = list(imputation_variables, imputation_variables))
  diag(predictor_matrix) <- 0L
  predictor_matrix[!imputed_by_model, ] <- 0L                              # 7a
  predictor_matrix[, intersect(routed_parts, imputation_variables)] <- 0L  # 7b
  if ("w3_intention_route16" %in% combined) {                              # 7c
    predictor_matrix[intersect(c("w3_continue_fte16", "w3_route_stage16",
                                 "w3_intended_domain16"), imputation_variables),
                     "w3_intention_route16"] <- 0L
    predictor_matrix[intersect(c("w3_route_stage16", "w3_intended_domain16"),
                               imputation_variables), "w3_continue_fte16"] <- 0L   # 7d
  }
  if ("he_pathway20" %in% combined) {
    predictor_matrix[c("he_participation20", "he_domain20"), "he_pathway20"] <- 0L
  }
  if ("degree_pathway25" %in% combined) {
    predictor_matrix[c("degree_attained25", "degree_domain25"), "degree_pathway25"] <- 0L
  }

  # d. Category check: training rows where the imputed variable is observed,
  # observed predictor values only. A bootstrap logistic model that fails it
  # can meet a bootstrap sample without a category, and its coefficients
  # can't then be estimated. For such a variable, in this training set only,
  # a classification tree is used instead (the CLS user guide's alternative
  # for sparse categories); every switch is listed in 16_category_check.csv.
  # First the (model, predictor) pairs to check, then one loop over them.
  check_pairs <- data.frame()
  for (v in names(method)[method == "logreg.boot"]) {
    predictors_v <- names(which(predictor_matrix[v, ] == 1))
    check_pairs <- rbind(check_pairs, data.frame(model_for = rep(v, length(predictors_v)),
                                                 predictor = predictors_v))
  }
  failing_models <- c()
  for (j in seq_len(nrow(check_pairs))) {
    fitted_rows <- training & !is.na(imputation_data[[check_pairs$model_for[j]]])
    values <- imputation_data[[check_pairs$predictor[j]]][fitted_rows]
    if (is.factor(values) && any(table(values) < minimum_category_size)) {
      counts <- table(values)
      category_check <- rbind(category_check, data.frame(
        repeat_number = r, outer_fold = f, landmark = landmark,
        model_for = check_pairs$model_for[j], predictor = check_pairs$predictor[j],
        category = names(counts)[counts < minimum_category_size],
        action = "classification tree used for this model"))
      failing_models <- c(failing_models, check_pairs$model_for[j])
    }
  }
  method[unique(failing_models)] <- "cart"

  # The check split's specification is kept for Blocks 6 and 7.
  if (is_check_split) {
    check_specs[[landmark]] <- list(
      training = training, test = test, features = features, combined = combined,
      imputation_data = imputation_data, where = where, method = method,
      predictor_matrix = predictor_matrix,
      visit_sequence = imputation_variables[method[imputation_variables] != ""])
  }
}

category_check$disclosure_checked <- rep(TRUE, nrow(category_check))
write.csv(category_check, file.path(review_dir, "16_category_check.csv"), row.names = FALSE)
screen_results$disclosure_checked <- TRUE
write.csv(screen_results, file.path(review_dir, "16_auxiliary_screen.csv"), row.names = FALSE)


###############################################################################
# Block 5. Workspace check (as Files 7 and 8a)
###############################################################################

for (name in c("state", "loggedEvents")) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) {
    stop("An object called '", name, "' is in the R workspace, and mice would ",
         "mistake it for its own. Remove it with rm(", name, ") or restart R.")
  }
}


###############################################################################
# Block 6. Leakage test of the imputation
#
# If the imputation models are fitted to training rows only, the imputed
# values of the training rows can't depend on the test rows' values. The
# L16 imputation is run twice with the same seed (1 imputation, 3
# iterations): once as it is, and once with every variable's observed values
# shuffled among the test rows (same missing pattern). The training rows'
# completed data must be identical.
###############################################################################

spec <- check_specs[["L16"]]
scrambled <- spec$imputation_data
set.seed(base_seed)
for (v in names(scrambled)) {
  observed_test <- which(spec$test & !is.na(scrambled[[v]]))
  if (length(observed_test) > 1) {
    scrambled[[v]][observed_test] <- scrambled[[v]][sample(observed_test)]
  }
}
leakage_runs <- list()
for (version in c("original", "scrambled")) {
  data_v <- spec$imputation_data
  if (version == "scrambled") data_v <- scrambled
  leakage_runs[[version]] <- suppressWarnings(mice(
    data_v, m = 1, maxit = 3, method = spec$method,
    predictorMatrix = spec$predictor_matrix, where = spec$where,
    ignore = spec$test, visitSequence = spec$visit_sequence,
    seed = base_seed, printFlag = FALSE))
}
training_original <- complete(leakage_runs$original, 1)[spec$training, ]
training_scrambled <- complete(leakage_runs$scrambled, 1)[spec$training, ]
leakage_free <- identical(training_original, training_scrambled)

checks <- data.frame(
  landmark = "L16", strategy = "primary",
  check = "Training-row imputations don't depend on test rows",
  passed = leakage_free)


###############################################################################
# Block 7. The pipeline on the check split
#
# The primary analysis at every landmark, and the outcome-omitted
# sensitivity analysis at L16. For each: imputation (with warnings kept),
# checks, convergence summaries, model fitting and prediction, and timing.
###############################################################################

runs <- rbind(
  data.frame(landmark = landmark_order, strategy = "primary"),
  data.frame(landmark = "L16", strategy = "outcome omitted")
)

convergence  <- data.frame()
model_checks <- data.frame()
timing       <- data.frame()
all_warnings <- data.frame()
trace_data   <- list()

for (k in seq_len(nrow(runs))) {
  landmark <- runs$landmark[k]
  strategy <- runs$strategy[k]
  spec     <- check_specs[[landmark]]
  data_k   <- spec$imputation_data
  method_k <- spec$method
  pm_k     <- spec$predictor_matrix
  where_k  <- spec$where
  visit_k  <- spec$visit_sequence

  # The outcome-omitted analysis leaves the outcome out of every equation.
  if (strategy == "outcome omitted") {
    keep <- setdiff(names(data_k), "outcome_for_imputation")
    data_k   <- data_k[, keep]
    method_k <- method_k[keep]
    pm_k     <- pm_k[keep, keep]
    where_k  <- where_k[, keep]
    visit_k  <- setdiff(visit_k, "outcome_for_imputation")
  }

  cat("\n", landmark, "(", strategy, "): imputing ...\n")
  start_time <- Sys.time()

  # Imputation, keeping its warnings. withCallingHandlers() runs mice(); each
  # time a warning occurs, the small function after "warning =" adds its
  # message to run_warnings (the "<<-" writes to that list outside the
  # function) and lets the imputation carry on.
  run_warnings <- character()
  imputation <- withCallingHandlers(
    mice(data_k, m = imputations, maxit = iterations, method = method_k,
         predictorMatrix = pm_k, where = where_k, ignore = spec$test,
         visitSequence = visit_k, seed = base_seed + k, printFlag = FALSE),
    warning = function(w) {
      run_warnings <<- c(run_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  imputation_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  # Checks on the completed data.
  # Primary: the test outcomes are missing before imputation. Outcome
  # omitted: the outcome isn't in the imputation at all.
  if (strategy == "primary") {
    masked <- all(is.na(data_k$outcome_for_imputation[spec$test]))
  } else {
    masked <- !("outcome_for_imputation" %in% names(data_k))
  }
  # The completed datasets, and one loop per check: observed values
  # unchanged (every imputation x variable); only structural values missing;
  # combined features equal their parts (every imputation x combined
  # feature).
  completed <- list()
  for (i in seq_len(imputations)) completed[[i]] <- complete(imputation, i)

  unchanged <- TRUE
  observed_table <- expand.grid(variable = names(data_k), imputation = seq_len(imputations),
                                stringsAsFactors = FALSE)
  for (j in seq_len(nrow(observed_table))) {
    v <- observed_table$variable[j]
    data_i <- completed[[observed_table$imputation[j]]]
    observed <- !is.na(data_k[[v]])
    if (any(as.character(data_i[[v]][observed]) != as.character(data_k[[v]][observed]))) {
      unchanged <- FALSE
    }
  }

  only_structural <- TRUE
  for (i in seq_len(imputations)) {
    if (any(is.na(completed[[i]]) != (is.na(data_k) & !where_k))) only_structural <- FALSE
  }

  passive_ok <- TRUE
  passive_table <- expand.grid(variable = names(method_k)[startsWith(method_k, "~")],
                               imputation = seq_len(imputations), stringsAsFactors = FALSE)
  for (j in seq_len(nrow(passive_table))) {
    v <- passive_table$variable[j]
    data_i <- completed[[passive_table$imputation[j]]]
    rebuilt <- model.frame(as.formula(method_k[[v]]), data_i, na.action = na.pass)[[1]]
    if (anyNA(rebuilt) || any(as.character(rebuilt) != as.character(data_i[[v]]))) {
      passive_ok <- FALSE
    }
  }
  checks <- rbind(checks, data.frame(
    landmark = landmark, strategy = strategy,
    check = c("Test outcomes unavailable to the imputation",
              "Observed values unchanged",
              "Only structural (not applicable) values still missing",
              "Combined features equal their parts"),
    passed = c(masked, unchanged, only_structural, passive_ok)))

  # Convergence over the second half of the iterations (as File 7), for
  # variables imputed by a model with at least 10 imputed values.
  n_imputed <- colSums(where_k)
  traced <- visit_k[!startsWith(method_k[visit_k], "~") &
                      n_imputed[visit_k] >= disclosure_threshold]
  # For each variable and statistic (mean, SD of the imputed values):
  # split R-hat over the second half of the iterations, and the trend across
  # iterations relative to the spread between chains. Flagged if split R-hat
  # is above 1.1, or if there is a clear trend (p < 0.01) larger than the
  # spread.
  second_half <- (floor(iterations / 2) + 1):iterations
  half_length <- floor(length(second_half) / 2)
  trace_table <- expand.grid(statistic = c("mean", "sd"), variable = traced, stringsAsFactors = FALSE)
  for (j in seq_len(nrow(trace_table))) {
    v <- trace_table$variable[j]
    statistic <- trace_table$statistic[j]
    values <- imputation$chainMean[v, , ]
    if (statistic == "sd") values <- sqrt(imputation$chainVar[v, , ])
    halves  <- cbind(values[head(second_half, half_length), , drop = FALSE],
                     values[tail(second_half, half_length), , drop = FALSE])
    n       <- nrow(halves)
    within  <- mean(apply(halves, 2, var))           # variance within each half-chain
    between <- n * var(colMeans(halves))             # variance between half-chain means
    split_rhat <- NA
    if (isTRUE(within > 0)) split_rhat <- sqrt(((n - 1) / n * within + between / n) / within)
    trend_data <- data.frame(value = as.vector(values[second_half, , drop = FALSE]),
                             iteration = rep(second_half, times = imputations),
                             chain = factor(rep(seq_len(imputations), each = length(second_half))))
    trend <- summary(lm(value ~ iteration + chain, data = trend_data))$coefficients
    spread <- mean(apply(values[second_half, , drop = FALSE], 1, sd))   # spread between chains
    drift <- NA
    if (isTRUE(spread > 0)) drift <- abs(trend["iteration", "Estimate"]) * (length(second_half) - 1) / spread
    convergence <- rbind(convergence, data.frame(
      landmark = landmark, strategy = strategy, variable = v, statistic = statistic,
      split_rhat = round(split_rhat, 3), trend_p = signif(trend["iteration", "Pr(>|t|)"], 3),
      drift_vs_spread = round(drift, 2),
      flag = isTRUE(split_rhat > 1.1) ||
        (isTRUE(trend["iteration", "Pr(>|t|)"] < 0.01) && isTRUE(drift > 1))))
  }
  trace_data[[paste(landmark, strategy)]] <- list(imputation = imputation, traced = traced)

  # Warnings and logged events.
  if (length(run_warnings) > 0) {
    all_warnings <- rbind(all_warnings, data.frame(
      landmark = landmark, strategy = strategy, source = "imputation warning",
      message = unique(run_warnings),
      times = as.integer(table(run_warnings)[unique(run_warnings)])))
  }
  if (!is.null(imputation$loggedEvents)) {
    all_warnings <- rbind(all_warnings, data.frame(
      landmark = landmark, strategy = strategy, source = "logged event",
      message = paste(imputation$loggedEvents$dep, imputation$loggedEvents$meth,
                      imputation$loggedEvents$out),
      times = 1L))
  }

  # Model fitting and prediction, on each imputed dataset.
  start_time <- Sys.time()
  y_training <- cohort[[outcome]][spec$training]
  inner_fold <- cohort[[paste0("inner_fold_r", check_repeat, "_f", check_fold)]][spec$training]
  if (anyNA(inner_fold)) stop("Inner folds don't cover the training rows.")
  predictions <- matrix(NA, nrow = sum(spec$test), ncol = 3,
                        dimnames = list(NULL, c("elastic_net", "logistic", "null")))
  prediction_sum <- matrix(0, nrow = sum(spec$test), ncol = 3)

  # Rare categories to pool, learned from the training rows' observed
  # values, so the same in every imputed dataset: categories with fewer than
  # 20 training people go into "Other (rare categories)", or into the most
  # common category if together they still have fewer than 20.
  pooling <- data.frame()
  for (v in spec$features) {
    if (!is.factor(cohort[[v]])) next
    counts <- table(cohort[[v]][spec$training])
    rare <- names(counts)[counts < minimum_category_size]
    new_label <- names(which.max(counts))
    if (sum(counts[rare]) >= minimum_category_size) new_label <- "Other (rare categories)"
    if (length(rare) > 0) {
      pooling <- rbind(pooling, data.frame(variable = v, category = rare, new_label = new_label))
    }
  }

  for (i in seq_len(imputations)) {
    features_i <- completed[[i]][, spec$features]

    # Rare categories pooled.
    for (j in seq_len(nrow(pooling))) {
      values <- as.character(features_i[[pooling$variable[j]]])
      values[values == pooling$category[j]] <- pooling$new_label[j]
      features_i[[pooling$variable[j]]] <- factor(values)
    }

    # A categorical feature left with a single category (after pooling)
    # carries no information in this training set and is left out.
    for (v in names(features_i)) {
      if (is.factor(features_i[[v]]) && nlevels(droplevels(features_i[[v]])) < 2) features_i[[v]] <- NULL
    }
    X <- model.matrix(~ ., data = features_i)[, -1, drop = FALSE]
    X_training <- X[spec$training, , drop = FALSE]
    X_test     <- X[spec$test, , drop = FALSE]

    # Elastic net tuned by grouped inner cross-validation: one fit per
    # alpha, then the alpha with the lowest cross-validated deviance (the
    # first, if tied).
    cv_fits <- list()
    deviances <- c()
    for (a in seq_along(alpha_grid)) {
      cv_fits[[a]] <- cv.glmnet(X_training, y_training, family = "binomial", alpha = alpha_grid[a],
                                foldid = inner_fold, type.measure = "deviance")
      deviances[a] <- min(cv_fits[[a]]$cvm)
    }
    best <- list(alpha = alpha_grid[which.min(deviances)], fit = cv_fits[[which.min(deviances)]])
    lambda_path <- best$fit$lambda
    p_elastic <- as.vector(predict(best$fit, newx = X_test, s = "lambda.min", type = "response"))

    # Unpenalized logistic regression, keeping its warnings (as for the
    # imputation above).
    glm_warnings <- character()
    glm_fit <- withCallingHandlers(
      glm.fit(cbind(1, X_training), y_training, family = binomial()),
      warning = function(w) {
        glm_warnings <<- c(glm_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    beta <- glm_fit$coefficients
    beta[is.na(beta)] <- 0
    p_logistic <- as.vector(plogis(cbind(1, X_test) %*% beta))

    # Null model: the training prevalence.
    p_null <- rep(mean(y_training), sum(spec$test))

    prediction_sum <- prediction_sum + cbind(p_elastic, p_logistic, p_null)

    model_checks <- rbind(model_checks, data.frame(
      landmark = landmark, strategy = strategy, imputation = i,
      columns = ncol(X_training),
      chosen_alpha = best$alpha,
      lambda_at_largest = best$fit$lambda.min == max(lambda_path),
      lambda_at_smallest = best$fit$lambda.min == min(lambda_path),
      nonzero_coefficients = sum(as.vector(as.matrix(coef(best$fit, s = "lambda.min")))[-1] != 0),
      logistic_converged = glm_fit$converged,
      logistic_aliased_columns = sum(is.na(glm_fit$coefficients)),
      logistic_warnings = length(glm_warnings)))
    if (length(glm_warnings) > 0) {
      all_warnings <- rbind(all_warnings, data.frame(
        landmark = landmark, strategy = strategy, source = "logistic benchmark warning",
        message = unique(glm_warnings),
        times = as.integer(table(glm_warnings)[unique(glm_warnings)])))
    }
  }
  predictions <- prediction_sum / imputations      # averaged (Approach 1)
  modelling_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  # Checks on the predictions (their form only; they are not scored).
  checks <- rbind(checks, data.frame(
    landmark = landmark, strategy = strategy,
    check = c("Predictions are complete and between 0 and 1",
              "Models fitted to training rows only"),
    passed = c(!anyNA(predictions) && all(predictions >= 0 & predictions <= 1),
               nrow(X_training) == sum(spec$training) &&
                 length(y_training) == sum(spec$training))))

  timing <- rbind(timing, data.frame(
    landmark = landmark, strategy = strategy,
    imputation_minutes = round(imputation_minutes, 1),
    modelling_minutes = round(modelling_minutes, 1)))
  cat("  done: imputation", round(imputation_minutes, 1), "min; models",
      round(modelling_minutes, 1), "min\n")
}


###############################################################################
# Block 8. Trace plots
#
# Mean and SD of the imputed values by iteration, one line per imputation,
# for variables with at least 10 imputed values.
###############################################################################

pdf(file.path(review_dir, "16_trace_plots.pdf"), width = 8.5, height = 11)
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
for (name in names(trace_data)) {
  imputation <- trace_data[[name]]$imputation
  for (v in trace_data[[name]]$traced) {
    matplot(imputation$chainMean[v, , ], type = "l", lty = 1, xlab = "Iteration",
            ylab = "Mean of imputed values", main = paste(name, v), cex.main = 0.8)
    matplot(sqrt(imputation$chainVar[v, , ]), type = "l", lty = 1, xlab = "Iteration",
            ylab = "SD of imputed values", main = paste(name, v), cex.main = 0.8)
  }
}
dev.off()


###############################################################################
# Block 9. Review tables and projected File 17a time
#
# Projected File 17a time: each primary landmark's time here, times the 25
# outer splits, plus the outcome-omitted analysis (repeat 1: 5 splits per
# landmark, taken as the primary time). Permutation importance in File 17a
# adds some modelling time.
###############################################################################

checks$disclosure_checked <- TRUE
write.csv(checks, file.path(review_dir, "16_checks.csv"), row.names = FALSE)
convergence$disclosure_checked <- rep(TRUE, nrow(convergence))
write.csv(convergence, file.path(review_dir, "16_convergence.csv"), row.names = FALSE)
model_checks$disclosure_checked <- TRUE
write.csv(model_checks, file.path(review_dir, "16_model_checks.csv"), row.names = FALSE)
if (nrow(all_warnings) == 0) {
  all_warnings <- data.frame(landmark = character(), strategy = character(),
                             source = character(), message = character(), times = integer())
}
all_warnings$disclosure_checked <- rep(TRUE, nrow(all_warnings))
write.csv(all_warnings, file.path(review_dir, "16_warnings.csv"), row.names = FALSE)

primary_timing <- timing[timing$strategy == "primary", ]
task_hours <- (primary_timing$imputation_minutes + primary_timing$modelling_minutes) / 60
projected_17a_hours <- sum(task_hours) * settings$n_repeats * settings$n_outer_folds +
  sum(task_hours) * settings$n_outer_folds
timing$disclosure_checked <- TRUE
write.csv(rbind(timing,
                data.frame(landmark = "All", strategy = "projected File 17a (hours, one core)",
                           imputation_minutes = NA, modelling_minutes = round(projected_17a_hours, 1),
                           disclosure_checked = TRUE)),
          file.path(review_dir, "16_timing.csv"), row.names = FALSE)


###############################################################################
# Block 10. Summary
#
# Checks and counts only. No predictions or performance.
###############################################################################

cat("\nFile 16 finished.\n")
cat("  Technical checks passed:              ", sum(checks$passed), "of", nrow(checks), "\n")
cat("  Training imputations independent of test rows:", leakage_free, "\n")
cat("  Sparse categories (all splits x landmarks):", nrow(category_check), "\n")
cat("  Convergence flags:                    ", sum(convergence$flag), "of",
    nrow(convergence), "rows\n")
cat("  Projected File 17a time (hours, one core):", round(projected_17a_hours, 1), "\n")
cat("  Review: 16_checks.csv, 16_category_check.csv, 16_convergence.csv,\n")
cat("    16_trace_plots.pdf, 16_model_checks.csv, 16_timing.csv, 16_warnings.csv,\n")
cat("    16_auxiliary_screen.csv\n")
