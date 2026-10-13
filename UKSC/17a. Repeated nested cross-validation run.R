###############################################################################
# 17a. Repeated nested cross-validation run
#
# Purpose
#   Fit and apply the prediction models in repeated, grouped, nested
#   cross-validation (File 14's folds). Each task is one analysis strategy x
#   repeat x outer fold x landmark:
#     tasks   1-125  primary analysis: 5 repeats x 5 folds x 5 landmarks;
#     tasks 126-150  outcome-omitted sensitivity analysis: repeat 1 only,
#                    identical except that the outcome is left out of every
#                    imputation equation (Sisk et al. 2023).
#   Each task saves out-of-fold predictions for its test fold, predictions
#   with each feature and each block of features permuted (for importance,
#   File 19), and technical diagnostics. Performance is computed in Files 18
#   and 19, not here.
#
# The per-fold pipeline (tested in File 16; the same steps)
#   1. Training rows are the outer training set; test rows are the outer
#      test fold.
#   2. Auxiliary screen (CLS user guide), on training rows only: each CLS
#      predictor of non-response known by the landmark is tested for
#      association with the landmark's most-missing feature (F test,
#      chi-square test or correlation test, by type), and kept if p < 0.05.
#   3. Multiple imputation of the landmark's features (Wahl et al. 2016;
#      Mertens et al. 2020, section 3.1): 5 imputations of 30 iterations.
#        - Variables: the features, the routed parts of the combined
#          features (recalculated from their parts at every iteration, as in
#          File 7), the kept auxiliaries and, in the primary analysis, the
#          outcome.
#        - The outcome is used for training rows and set to missing for test
#          rows (CLS user guide: the outcome belongs in the imputation
#          model); the imputed test outcomes are discarded. In the
#          outcome-omitted analysis the outcome isn't in the imputation.
#        - The imputation models are fitted to training rows only (mice's
#          "ignore"); test rows are imputed from those models.
#        - Methods, routing and predictor rules as in File 7, except that
#          there are no design variables (see "Survey weights" below):
#          predictive mean matching (numeric), bootstrap logistic
#          regression (two categories), classification trees (more
#          categories); structural
#          (not applicable) values are never imputed; values marked for
#          imputation are set to missing first (partial activity counts at
#          16-19 are imputed over the whole period). Where a two-category
#          variable's predictors have a category with fewer than 20 people
#          among that training set's rows, a classification tree is used for
#          it instead (checked on training rows only). Main activity at 25 is
#          grouped into File 7's six groups.
#   4. For each imputed dataset: rare categories (fewer than 20 training
#      people) are pooled, learned from training rows only; the elastic net
#      is tuned by grouped inner cross-validation (File 14's PSU-whole inner
#      folds) over alpha in {0, 0.25, 0.5, 0.75, 1} and the lambda path, on
#      binomial deviance, at the minimum; the elastic net, an unpenalized
#      logistic regression and a null (prevalence) model are fitted to the
#      training rows; test-row probabilities are predicted.
#   5. Predictions are averaged over the imputed datasets (Mertens et al.
#      2020, Approach 1). Models are fitted unweighted; performance (File 18)
#      is weighted by the Wave-9 final weight.
#   6. Permutation importance (predictions only): for the elastic net, each
#      feature, and each block of features, is shuffled among the test rows
#      (one random order per feature or block, the same in every imputed
#      dataset), and predictions are made again and averaged. File 19 scores
#      the loss of performance.
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
# Running on several computers or R sessions
#   Give each computer or R session its own tasks_to_run (Block 1), e.g.
#   1:50, 51:100 and 101:150; each block of 50 has a similar mix of
#   landmarks, so the loads are similar. A task whose file exists is
#   skipped. Each task's seeds are fixed (base seed + task number), so
#   results don't depend on which computer ran which task. Copy every task
#   file into tasks_dir/17a on one computer and run File 17b there. File 0
#   records package versions; File 17b combines tasks only if they match.
#   Only copy data between computers covered by the UKDS End User Licence
#   and your institution's data security rules.
#
# Reads
#   data_dir/14_prediction_design.rds
#
# Writes
#   tasks_dir/17a/task_NNN.rds   one per task: out-of-fold predictions,
#       permuted predictions and diagnostics (participant level; never shared)
#
# References
#   Mertens BJA, Banzato E, de Wreede LC (2020). Biometrical Journal 62(3):
#   724-741.
#   Sisk R, et al. (2023). Statistical Methods in Medical Research 32(8):
#   1461-1477.
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

paths     <- readRDS("project_paths.rds")
data_dir  <- paths$data_dir     # datasets passed between scripts
tasks_dir <- paths$tasks_dir    # one file per task in split runs

# Which tasks this computer (or R session) runs: "all", or task numbers,
# e.g. 1:50.
tasks_to_run <- "all"

imputations           <- 5       # imputations per outer training set
iterations            <- 30      # iterations per imputation
alpha_grid            <- c(0, 0.25, 0.5, 0.75, 1)   # elastic net mixing
minimum_category_size <- 20      # bootstrap-logistic check and rare-level pooling
screen_p_cutoff       <- 0.05    # auxiliary screen
base_seed             <- 1700    # task t uses seed base_seed + t

task_folder <- file.path(tasks_dir, "17a")
dir.create(task_folder, recursive = TRUE, showWarnings = FALSE)


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
# Block 3. Prepare the variables (as File 16)
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
# Block 4. Task table
#
# Primary: repeat x fold x landmark (125 tasks); then the outcome-omitted
# sensitivity analysis for repeat 1 (25 tasks).
###############################################################################

task_table <- rbind(
  expand.grid(landmark = landmark_order, outer_fold = seq_len(settings$n_outer_folds),
              repeat_number = seq_len(settings$n_repeats), strategy = "primary",
              stringsAsFactors = FALSE),
  expand.grid(landmark = landmark_order, outer_fold = seq_len(settings$n_outer_folds),
              repeat_number = 1, strategy = "outcome omitted",
              stringsAsFactors = FALSE)
)
task_table <- task_table[, c("strategy", "repeat_number", "outer_fold", "landmark")]
task_table$task <- seq_len(nrow(task_table))
task_table$file <- file.path(task_folder, sprintf("task_%03d.rds", task_table$task))

if (identical(tasks_to_run, "all")) tasks_to_run <- task_table$task
if (!all(tasks_to_run %in% task_table$task)) {
  stop("tasks_to_run includes task numbers outside 1-", nrow(task_table), ".")
}


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
# Block 6. Run the tasks
#
# Each task is one commented block of steps a-j. A failed check stops the
# script and names the task.
###############################################################################

cat("\nCross-validation run:", length(tasks_to_run), "of", nrow(task_table),
    "tasks on this computer.\n")

for (t in tasks_to_run) {
  task_file <- task_table$file[t]
  if (file.exists(task_file)) {
    cat("  Task", t, "already done; skipped.\n")
    next
  }
  strategy <- task_table$strategy[t]
  r        <- task_table$repeat_number[t]
  f        <- task_table$outer_fold[t]
  landmark <- task_table$landmark[t]
  cat("  Task", t, "(", strategy, ", repeat", r, ", fold", f, ",", landmark, ") started",
      format(Sys.time(), "%H:%M"), "\n")
  start_time <- Sys.time()

  # a. Training and test rows, and the landmark's variables.
  training <- cohort[[paste0("outer_fold_r", r)]] != f
  test     <- !training
  rows_l   <- landmark_variables[landmark_variables$landmark == landmark, ]
  features <- rows_l$variable[rows_l$role == "Feature"]
  part_of  <- rows_l[startsWith(rows_l$role, "Part of"), ]
  combined <- unique(sub("^Part of ", "", part_of$role))
  auxiliary_candidates <- rows_l$variable[rows_l$role == "Auxiliary"]

  # b. Auxiliary screen on training rows, against the most-missing feature.
  # The test depends on the types: correlation (both numeric), F test (one
  # numeric), chi-square test (both categorical).
  missing_counts <- c()
  for (v in features) missing_counts[v] <- sum(is.na(cohort[[v]][training]))
  screen_target  <- names(which.max(missing_counts))
  y_values <- cohort[[screen_target]][training]
  kept_auxiliaries <- c()
  screen <- data.frame()
  for (candidate in auxiliary_candidates) {
    x_values <- cohort[[candidate]][training]
    both <- !is.na(y_values) & !is.na(x_values)
    y <- y_values[both]
    x <- x_values[both]
    if (is.factor(y)) y <- droplevels(y)
    if (is.factor(x)) x <- droplevels(x)
    can_test <- length(unique(y)) > 1 && length(unique(x)) > 1
    p_value <- NA
    if (can_test && is.numeric(y) && is.numeric(x)) {
      p_value <- cor.test(y, x)$p.value
    }
    if (can_test && is.numeric(y) && !is.numeric(x)) {
      p_value <- anova(lm(y ~ x))[["Pr(>F)"]][1]
    }
    if (can_test && !is.numeric(y) && is.numeric(x)) {
      p_value <- anova(lm(x ~ y))[["Pr(>F)"]][1]
    }
    if (can_test && !is.numeric(y) && !is.numeric(x)) {
      p_value <- suppressWarnings(chisq.test(table(y, x))$p.value)
    }
    keep <- isTRUE(p_value < screen_p_cutoff)
    if (keep) kept_auxiliaries <- c(kept_auxiliaries, candidate)
    screen <- rbind(screen, data.frame(screened_against = screen_target,
                                       auxiliary = candidate, p_value = p_value, kept = keep))
  }

  # c. The imputation variables, in visiting order: auxiliaries, then the
  # features (each combined feature after its parts), then (primary only)
  # the outcome.
  imputation_variables <- kept_auxiliaries
  for (feature in features) {
    if (feature %in% combined) {
      parts <- part_of$variable[part_of$role == paste("Part of", feature)]
      imputation_variables <- c(imputation_variables, setdiff(parts, imputation_variables))
    }
    imputation_variables <- c(imputation_variables, setdiff(feature, imputation_variables))
  }
  if (strategy == "primary") {
    imputation_variables <- c(imputation_variables, "outcome_for_imputation")
  }

  imputation_data <- cohort[, imputation_variables]
  if (strategy == "primary") {
    imputation_data$outcome_for_imputation[test] <- NA     # masked for test rows
  }

  # d. Where (item and wave non-response, and the masked outcome); values
  # marked for imputation set to missing; checks.
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
  for (v in imputation_variables) {
    imputation_data[[v]][where[, v]] <- NA
  }
  for (v in imputation_variables) {
    state_name <- paste0("state_", v)
    structural <- rep(FALSE, nrow(cohort))
    if (state_name %in% names(cohort)) {
      structural <- cohort[[state_name]] == "Structural (not applicable)"
    }
    if (any(is.na(imputation_data[[v]]) & !where[, v] & !structural)) {
      stop("Task ", t, ": ", v, " has missing values that are neither imputed nor structural.")
    }
    if (any(where[, v]) && sum(training & !is.na(imputation_data[[v]])) == 0) {
      stop("Task ", t, ": ", v, " has no observed values among the training rows.")
    }
  }

  # e. Methods and predictor matrix (as File 16), with the category check
  # and the classification-tree fallback.
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

  # Category check: for each variable imputed by bootstrap logistic
  # regression, the categories of each categorical predictor among the
  # training rows its model is fitted to. A category with fewer than 20
  # people switches that variable to a classification tree. First the
  # (model, predictor) pairs to check, then one loop over them.
  check_pairs <- data.frame()
  for (v in names(method)[method == "logreg.boot"]) {
    predictors_v <- names(which(predictor_matrix[v, ] == 1))
    check_pairs <- rbind(check_pairs, data.frame(model_for = rep(v, length(predictors_v)),
                                                 predictor = predictors_v))
  }
  switches <- data.frame()
  for (j in seq_len(nrow(check_pairs))) {
    fitted_rows <- training & !is.na(imputation_data[[check_pairs$model_for[j]]])
    values <- imputation_data[[check_pairs$predictor[j]]][fitted_rows]
    if (is.factor(values) && any(table(values) < minimum_category_size)) {
      counts <- table(values)
      switches <- rbind(switches, data.frame(model_for = check_pairs$model_for[j],
                                             predictor = check_pairs$predictor[j],
                                             category = names(counts)[counts < minimum_category_size]))
    }
  }
  if (nrow(switches) > 0) method[unique(switches$model_for)] <- "cart"
  visit_sequence <- imputation_variables[method[imputation_variables] != ""]

  # f. Imputation (training rows fit the models; test rows are imputed),
  # keeping its warnings. withCallingHandlers() runs mice(); each time a
  # warning occurs, the small function after "warning =" adds its message to
  # task_warnings (the "<<-" writes to that list outside the function) and
  # lets the imputation carry on.
  task_warnings <- character()
  imputation <- withCallingHandlers(
    mice(imputation_data, m = imputations, maxit = iterations, method = method,
         predictorMatrix = predictor_matrix, where = where, ignore = test,
         visitSequence = visit_sequence, seed = base_seed + t, printFlag = FALSE),
    warning = function(w) {
      task_warnings <<- c(task_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  imputation_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  # g. Checks on the completed data (as File 16): observed values unchanged
  # (every imputation x variable); only structural values missing; combined
  # features equal their parts (every imputation x combined feature).
  completed <- list()
  for (i in seq_len(imputations)) completed[[i]] <- complete(imputation, i)

  observed_table <- expand.grid(variable = names(imputation_data), imputation = seq_len(imputations),
                                stringsAsFactors = FALSE)
  for (j in seq_len(nrow(observed_table))) {
    v <- observed_table$variable[j]
    data_i <- completed[[observed_table$imputation[j]]]
    observed <- !is.na(imputation_data[[v]])
    if (any(as.character(data_i[[v]][observed]) != as.character(imputation_data[[v]][observed]))) {
      stop("Task ", t, ": observed values of ", v, " changed.")
    }
  }
  for (i in seq_len(imputations)) {
    if (any(is.na(completed[[i]]) != (is.na(imputation_data) & !where))) {
      stop("Task ", t, ": values other than structural ones are missing after imputation.")
    }
  }
  passive_table <- expand.grid(variable = names(method)[is_passive], imputation = seq_len(imputations),
                               stringsAsFactors = FALSE)
  for (j in seq_len(nrow(passive_table))) {
    v <- passive_table$variable[j]
    data_i <- completed[[passive_table$imputation[j]]]
    rebuilt <- model.frame(as.formula(method[[v]]), data_i, na.action = na.pass)[[1]]
    if (anyNA(rebuilt) || any(as.character(rebuilt) != as.character(data_i[[v]]))) {
      stop("Task ", t, ": ", v, " doesn't equal its parts.")
    }
  }

  # Convergence summary over the second half of the iterations (as File 7),
  # for each variable (at least 10 imputed values) and statistic (mean, SD
  # of the imputed values): split R-hat, and the trend across iterations
  # relative to the spread between chains. Flagged if split R-hat is above
  # 1.1, or if there is a clear trend (p < 0.01) larger than the spread.
  convergence <- data.frame()
  n_imputed <- colSums(where)
  traced <- visit_sequence[!startsWith(method[visit_sequence], "~") &
                             n_imputed[visit_sequence] >= 10]
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
      variable = v, statistic = statistic, split_rhat = split_rhat,
      trend_p = trend["iteration", "Pr(>|t|)"], drift_vs_spread = drift,
      flag = isTRUE(split_rhat > 1.1) ||
        (isTRUE(trend["iteration", "Pr(>|t|)"] < 0.01) && isTRUE(drift > 1))))
  }

  # h. Models on each imputed dataset, and test-row predictions.
  y_training <- cohort[[outcome]][training]
  inner_fold <- cohort[[paste0("inner_fold_r", r, "_f", f)]][training]
  if (anyNA(inner_fold)) stop("Task ", t, ": inner folds don't cover the training rows.")

  blocks <- unique(rows_l$block[rows_l$role == "Feature"])
  block_members <- c()
  for (b in blocks) {
    block_members[b] <- paste(rows_l$variable[rows_l$role == "Feature" & rows_l$block == b],
                              collapse = " | ")
  }
  permutation_units <- rbind(
    data.frame(unit = paste("Feature:", features), type = "feature", members = features),
    data.frame(unit = paste("Block:", blocks), type = "block", members = unname(block_members))
  )
  # One random order of the test rows per unit, the same in every imputation.
  set.seed(base_seed + 100000 + t)
  permutation_orders <- list()
  for (u in permutation_units$unit) permutation_orders[[u]] <- sample(sum(test))

  prediction_sum  <- matrix(0, nrow = sum(test), ncol = 3,
                            dimnames = list(NULL, c("elastic_net", "logistic", "null")))
  permutation_sum <- matrix(0, nrow = sum(test), ncol = nrow(permutation_units),
                            dimnames = list(NULL, permutation_units$unit))
  model_checks <- data.frame()

  # Rare categories to pool, learned from the training rows' observed
  # values, so the same in every imputed dataset: categories with fewer than
  # 20 training people go into "Other (rare categories)", or into the most
  # common category if together they still have fewer than 20.
  pooling <- data.frame()
  for (v in features) {
    if (!is.factor(cohort[[v]])) next
    counts <- table(cohort[[v]][training])
    rare <- names(counts)[counts < minimum_category_size]
    new_label <- names(which.max(counts))
    if (sum(counts[rare]) >= minimum_category_size) new_label <- "Other (rare categories)"
    if (length(rare) > 0) {
      pooling <- rbind(pooling, data.frame(variable = v, category = rare, new_label = new_label))
    }
  }

  for (i in seq_len(imputations)) {
    features_i <- completed[[i]][, features]

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

    # Model matrix. model.matrix records, for each column, the number of the
    # feature it came from ("assign"), which gives each column's feature.
    X <- model.matrix(~ ., data = features_i)
    column_feature <- names(features_i)[attr(X, "assign")[-1]]
    X <- X[, -1, drop = FALSE]
    X_training <- X[training, , drop = FALSE]
    X_test     <- X[test, , drop = FALSE]

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
    p_elastic <- as.vector(predict(best$fit, newx = X_test, s = "lambda.min", type = "response"))

    # Unpenalized logistic regression (keeping its warnings, as for the
    # imputation above) and the null model.
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
    p_null <- rep(mean(y_training), sum(test))

    prediction_sum <- prediction_sum + cbind(p_elastic, p_logistic, p_null)

    # Permutation predictions (elastic net).
    for (k in seq_len(nrow(permutation_units))) {
      members <- strsplit(permutation_units$members[k], " | ", fixed = TRUE)[[1]]
      columns <- which(column_feature %in% members)
      X_permuted <- X_test
      X_permuted[, columns] <- X_test[permutation_orders[[k]], columns, drop = FALSE]
      permutation_sum[, k] <- permutation_sum[, k] +
        as.vector(predict(best$fit, newx = X_permuted, s = "lambda.min", type = "response"))
    }

    lambda_path <- best$fit$lambda
    model_checks <- rbind(model_checks, data.frame(
      imputation = i, columns = ncol(X_training), chosen_alpha = best$alpha,
      lambda_at_largest = best$fit$lambda.min == max(lambda_path),
      lambda_at_smallest = best$fit$lambda.min == min(lambda_path),
      nonzero_coefficients = sum(as.vector(as.matrix(coef(best$fit, s = "lambda.min")))[-1] != 0),
      logistic_converged = glm_fit$converged,
      logistic_aliased_columns = sum(is.na(glm_fit$coefficients)),
      logistic_warnings = length(glm_warnings)))
    task_warnings <- c(task_warnings, glm_warnings)
  }
  predictions_average  <- prediction_sum / imputations      # averaged (Approach 1)
  permutation_average  <- permutation_sum / imputations
  if (anyNA(predictions_average) || any(predictions_average < 0 | predictions_average > 1)) {
    stop("Task ", t, ": predictions are missing or outside 0-1.")
  }

  # i. Out-of-fold predictions, with whether any feature (or routed part)
  # was imputed for the test row (for File 18's calibration split).
  feature_columns <- intersect(c(features, part_of$variable), colnames(where))
  any_imputed <- rowSums(where[test, feature_columns, drop = FALSE]) > 0
  predictions <- data.frame(
    NSID = cohort$NSID[test], outcome = cohort[[outcome]][test],
    w9_final_weight = cohort$w9_final_weight[test],
    psu = cohort$psu[test], stratum = cohort$stratum[test],
    strategy = strategy, repeat_number = r, outer_fold = f, landmark = landmark,
    p_elastic_net = predictions_average[, "elastic_net"],
    p_logistic = predictions_average[, "logistic"],
    p_null = predictions_average[, "null"],
    all_features_observed = !any_imputed
  )

  # j. Save the task.
  task_result <- list(
    task = t, strategy = strategy, repeat_number = r, outer_fold = f, landmark = landmark,
    predictions = predictions,
    permutation_units = permutation_units,
    permuted_predictions = permutation_average,
    screen = screen, switches = switches, convergence = convergence,
    model_checks = model_checks,
    logged_events = if (is.null(imputation$loggedEvents)) 0L else nrow(imputation$loggedEvents),
    warnings = table(task_warnings),
    imputation_minutes = imputation_minutes,
    total_minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins")),
    r_version = as.character(getRversion()),
    mice_version = as.character(packageVersion("mice")),
    glmnet_version = as.character(packageVersion("glmnet")),
    computer = gsub("[^A-Za-z0-9_-]", "_", Sys.info()[["nodename"]])
  )
  temporary_file <- paste0(task_file, ".partial")
  saveRDS(task_result, temporary_file)
  file.rename(temporary_file, task_file)
  cat("  Task", t, "finished in", round(task_result$total_minutes, 1), "minutes\n")
}


###############################################################################
# Block 7. Summary
#
# Task counts only. No predictions or performance.
###############################################################################

tasks_done <- sum(file.exists(task_table$file))
cat("\nFile 17a finished on this computer.\n")
cat("  Tasks in the table:        ", nrow(task_table), "\n")
cat("  Task files now present:    ", tasks_done, "\n")
if (tasks_done == nrow(task_table)) cat("  All tasks done. Next: File 17b.\n")
