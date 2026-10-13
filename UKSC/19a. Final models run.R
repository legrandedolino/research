###############################################################################
# 19a. Final models run
#
# Purpose
#   Fit each landmark's final elastic net to the whole prediction cohort, so
#   that its coefficients can be reported (TRIPOD+AI 22). The steps are
#   those of the cross-validation's primary analysis (File 17a), with every
#   person in the training set and no test rows:
#     1. Auxiliary screen (CLS user guide): each CLS predictor of
#        non-response known by the landmark is tested for association with
#        the landmark's most-missing feature (F test, chi-square test or
#        correlation test, by type), and kept if p < 0.05.
#     2. Multiple imputation of the landmark's features: 5 imputations of 30
#        iterations. The outcome is in the imputation model (CLS user
#        guide). Methods, routing and predictor rules as in File 7, except
#        that there are no design variables (see "Survey weights" below):
#        predictive mean matching (numeric), bootstrap logistic regression
#        (two categories), classification trees (more categories);
#        structural (not applicable) values are never imputed; values
#        marked for imputation are set to missing first. Where a
#        two-category variable's predictors have a category with fewer than
#        20 people, a classification tree is used for it instead.
#     3. For each imputed dataset: rare categories (fewer than 20 people)
#        pooled; a categorical feature left with one category dropped; the
#        elastic net tuned by grouped cross-validation over alpha in
#        {0, 0.25, 0.5, 0.75, 1} and the lambda path (binomial deviance,
#        lambda at the minimum), with repeat 1's outer folds from File 14
#        (whole schools) as the tuning folds.
#   File 19b combines the coefficients over the imputed datasets. The
#   models describe how the predictions are formed. They are not offered
#   for use with individuals, and their performance is the cross-validated
#   performance of File 18.
#
# Survey weights and the imputation model (a deliberate departure)
#   Silverwood et al. (2024) weighted their imputation models by the initial
#   design weights, and the imputation analysis (File 7) includes the log
#   design weight and the sampling stratum as predictors in every equation.
#   Neither is done here, as in the cross-validation (File 17a): these
#   imputation models are unweighted and contain no design variables. An
#   imputation model should match the analysis that uses the completed
#   data. The final models are fitted unweighted, and an imputation model
#   without design variables matches them. The design oversampled deprived
#   schools and minority ethnic pupils; ethnic group and family background
#   are features at every landmark.
#
# Running on several computers or R sessions
#   One task per landmark (5 tasks). Give each computer or R session its
#   own tasks_to_run (Block 1), e.g. c(1, 2, 5) and c(3, 4), which give
#   similar loads. A task whose file exists is skipped. Each task's seed is
#   fixed (base seed + task number), so results don't depend on which
#   computer ran it. Copy every task file into tasks_dir/19a on one computer
#   and run File 19b there. Only copy data between computers covered by the
#   UKDS End User Licence and your institution's data security rules.
#
# Warnings
#   Any warnings from the imputation or the models are printed at the end
#   of the run. They contain no data values and can be shared for review.
#
# Reads
#   data_dir/14_prediction_design.rds
#
# Writes
#   tasks_dir/19a/task_N.rds   one per landmark: coefficients, counts per
#       model term, tuning and diagnostics (private; never shared)
#
# References
#   Collins GS, et al. (2024). TRIPOD+AI statement. BMJ 385: e078378.
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
# e.g. c(1, 2, 5).
tasks_to_run <- "all"

imputations           <- 5       # imputations (as in File 17a)
iterations            <- 30      # iterations per imputation (as in File 17a)
alpha_grid            <- c(0, 0.25, 0.5, 0.75, 1)   # elastic net mixing
minimum_category_size <- 20      # bootstrap-logistic check and rare-level pooling
screen_p_cutoff       <- 0.05    # auxiliary screen
base_seed             <- 1900    # task t uses seed base_seed + t

task_folder <- file.path(tasks_dir, "19a")
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

n <- nrow(cohort)
outcome_values <- cohort[[outcome]]            # 1 = SRET at 32, 0 = not
tuning_folds   <- cohort$outer_fold_r1         # whole-school folds from File 14


###############################################################################
# Block 3. Prepare the variables (as File 17a)
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
# Block 4. Task table: one task per landmark
###############################################################################

task_table <- data.frame(task = seq_along(landmark_order), landmark = landmark_order)
task_table$file <- file.path(task_folder, sprintf("task_%d.rds", task_table$task))

if (identical(tasks_to_run, "all")) tasks_to_run <- task_table$task
if (!all(tasks_to_run %in% task_table$task)) {
  stop("tasks_to_run includes task numbers outside 1-", nrow(task_table), ".")
}


###############################################################################
# Block 5. Workspace check (as File 17a)
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
# Each task is one commented block of steps a-i. A failed check stops the
# script and names the task.
###############################################################################

cat("\nFinal models:", length(tasks_to_run), "of", nrow(task_table), "tasks on this computer.\n")

for (t in tasks_to_run) {
  task_file <- task_table$file[t]
  if (file.exists(task_file)) {
    cat("  Task", t, "already done; skipped.\n")
    next
  }
  landmark <- task_table$landmark[t]
  cat("  Task", t, "(", landmark, ") started", format(Sys.time(), "%H:%M"), "\n")
  start_time <- Sys.time()

  # a. The landmark's variables.
  rows_l   <- landmark_variables[landmark_variables$landmark == landmark, ]
  features <- rows_l$variable[rows_l$role == "Feature"]
  part_of  <- rows_l[startsWith(rows_l$role, "Part of"), ]
  combined <- unique(sub("^Part of ", "", part_of$role))
  auxiliary_candidates <- rows_l$variable[rows_l$role == "Auxiliary"]

  # b. Auxiliary screen against the most-missing feature. The test depends
  # on the types: correlation (both numeric), F test (one numeric),
  # chi-square test (both categorical).
  missing_counts <- c()
  for (v in features) missing_counts[v] <- sum(is.na(cohort[[v]]))
  screen_target <- names(which.max(missing_counts))
  kept_auxiliaries <- c()
  screen <- data.frame()
  for (candidate in auxiliary_candidates) {
    both <- !is.na(cohort[[screen_target]]) & !is.na(cohort[[candidate]])
    y_s <- cohort[[screen_target]][both]
    x_s <- cohort[[candidate]][both]
    if (is.factor(y_s)) y_s <- droplevels(y_s)
    if (is.factor(x_s)) x_s <- droplevels(x_s)
    can_test <- length(unique(y_s)) > 1 && length(unique(x_s)) > 1
    p_value <- NA
    if (can_test && is.numeric(y_s) && is.numeric(x_s)) {
      p_value <- cor.test(y_s, x_s)$p.value
    }
    if (can_test && is.numeric(y_s) && !is.numeric(x_s)) {
      p_value <- anova(lm(y_s ~ x_s))[["Pr(>F)"]][1]
    }
    if (can_test && !is.numeric(y_s) && is.numeric(x_s)) {
      p_value <- anova(lm(x_s ~ y_s))[["Pr(>F)"]][1]
    }
    if (can_test && !is.numeric(y_s) && !is.numeric(x_s)) {
      p_value <- suppressWarnings(chisq.test(table(y_s, x_s))$p.value)
    }
    keep <- isTRUE(p_value < screen_p_cutoff)
    if (keep) kept_auxiliaries <- c(kept_auxiliaries, candidate)
    screen <- rbind(screen, data.frame(screened_against = screen_target, auxiliary = candidate,
                                       p_value = p_value, kept = keep))
  }

  # c. The imputation variables, in visiting order: auxiliaries, then the
  # features (each combined feature after its parts), then the outcome.
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

  # d. Where to impute (item and wave non-response); values marked for
  # imputation set to missing first; checks.
  where <- matrix(FALSE, nrow = n, ncol = length(imputation_variables),
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
    structural <- rep(FALSE, n)
    if (state_name %in% names(cohort)) {
      structural <- cohort[[state_name]] == "Structural (not applicable)"
    }
    if (any(is.na(imputation_data[[v]]) & !where[, v] & !structural)) {
      stop("Task ", t, ": ", v, " has missing values that are neither imputed nor structural.")
    }
    if (any(where[, v]) && all(is.na(imputation_data[[v]]))) {
      stop("Task ", t, ": ", v, " has no observed values.")
    }
  }

  # e. Methods and predictor matrix (as File 17a).
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
  is_passive       <- setNames(startsWith(method, "~"), names(method))
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
  # regression, the categories of each categorical predictor among the rows
  # its model is fitted to. A category with fewer than 20 people switches
  # that variable to a classification tree.
  check_pairs <- data.frame()
  for (v in names(method)[method == "logreg.boot"]) {
    predictors_v <- names(which(predictor_matrix[v, ] == 1))
    check_pairs <- rbind(check_pairs, data.frame(model_for = rep(v, length(predictors_v)),
                                                 predictor = predictors_v))
  }
  switches <- data.frame()
  for (k in seq_len(nrow(check_pairs))) {
    fitted_rows <- !is.na(imputation_data[[check_pairs$model_for[k]]])
    values <- imputation_data[[check_pairs$predictor[k]]][fitted_rows]
    if (is.factor(values) && any(table(values) < minimum_category_size)) {
      counts <- table(values)
      switches <- rbind(switches, data.frame(model_for = check_pairs$model_for[k],
                                             predictor = check_pairs$predictor[k],
                                             category = names(counts)[counts < minimum_category_size]))
    }
  }
  if (nrow(switches) > 0) method[unique(switches$model_for)] <- "cart"
  visit_sequence <- imputation_variables[method[imputation_variables] != ""]

  # f. Imputation.
  imputation <- mice(imputation_data, m = imputations, maxit = iterations, method = method,
                     predictorMatrix = predictor_matrix, where = where,
                     visitSequence = visit_sequence, seed = base_seed + t, printFlag = FALSE)
  imputation_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  # g. Checks on the completed data (as File 17a): observed values
  # unchanged; only structural values missing; combined features equal
  # their parts.
  completed <- list()
  for (i in seq_len(imputations)) completed[[i]] <- complete(imputation, i)

  for (i in seq_len(imputations)) {
    if (any(is.na(completed[[i]]) != (is.na(imputation_data) & !where))) {
      stop("Task ", t, ": values other than structural ones are missing after imputation.")
    }
  }
  check_table <- expand.grid(imputation = seq_len(imputations), variable = imputation_variables,
                             stringsAsFactors = FALSE)
  for (k in seq_len(nrow(check_table))) {
    v <- check_table$variable[k]
    observed <- !is.na(imputation_data[[v]])
    values <- completed[[check_table$imputation[k]]][[v]]
    if (any(as.character(values[observed]) != as.character(imputation_data[[v]][observed]))) {
      stop("Task ", t, ": observed values of ", v, " changed.")
    }
  }
  passive_table <- expand.grid(imputation = seq_len(imputations),
                               variable = names(method)[is_passive], stringsAsFactors = FALSE)
  for (k in seq_len(nrow(passive_table))) {
    data_k <- completed[[passive_table$imputation[k]]]
    v <- passive_table$variable[k]
    rebuilt <- model.frame(as.formula(method[[v]]), data_k, na.action = na.pass)[[1]]
    if (anyNA(rebuilt) || any(as.character(rebuilt) != as.character(data_k[[v]]))) {
      stop("Task ", t, ": ", v, " doesn't equal its parts.")
    }
  }

  # Convergence summary over the second half of the iterations (as File
  # 17a): split R-hat, and the trend across iterations relative to the
  # spread between chains. Variables with at least 10 imputed values.
  n_imputed <- colSums(where)
  traced <- visit_sequence[!is_passive[visit_sequence] & n_imputed[visit_sequence] >= 10]
  second_half <- (floor(iterations / 2) + 1):iterations
  half_length <- floor(length(second_half) / 2)
  convergence <- expand.grid(variable = traced, statistic = c("mean", "sd"), stringsAsFactors = FALSE)
  convergence$split_rhat      <- rep(NA, nrow(convergence))
  convergence$trend_p         <- rep(NA, nrow(convergence))
  convergence$drift_vs_spread <- rep(NA, nrow(convergence))
  convergence$flag            <- rep(NA, nrow(convergence))
  for (k in seq_len(nrow(convergence))) {
    v <- convergence$variable[k]
    values <- imputation$chainMean[v, , ]
    if (convergence$statistic[k] == "sd") values <- sqrt(imputation$chainVar[v, , ])
    halves  <- cbind(values[head(second_half, half_length), , drop = FALSE],
                     values[tail(second_half, half_length), , drop = FALSE])
    within  <- mean(apply(halves, 2, var))
    between <- nrow(halves) * var(colMeans(halves))
    if (isTRUE(within > 0)) {
      convergence$split_rhat[k] <- sqrt(((nrow(halves) - 1) / nrow(halves) * within +
                                           between / nrow(halves)) / within)
    }
    trend_data <- data.frame(value = as.vector(values[second_half, , drop = FALSE]),
                             iteration = rep(second_half, times = imputations),
                             chain = factor(rep(seq_len(imputations), each = length(second_half))))
    trend <- summary(lm(value ~ iteration + chain, data = trend_data))$coefficients
    spread <- mean(apply(values[second_half, , drop = FALSE], 1, sd))
    convergence$trend_p[k] <- trend["iteration", "Pr(>|t|)"]
    if (isTRUE(spread > 0)) {
      convergence$drift_vs_spread[k] <- abs(trend["iteration", "Estimate"]) *
        (length(second_half) - 1) / spread
    }
    # Flagged if split R-hat is above 1.1, or if there is a clear trend
    # (p < 0.01) larger than the spread between chains.
    convergence$flag[k] <- isTRUE(convergence$split_rhat[k] > 1.1) ||
      (isTRUE(convergence$trend_p[k] < 0.01) && isTRUE(convergence$drift_vs_spread[k] > 1))
  }

  # h. The elastic net on each imputed dataset.

  # Rare categories to pool, from the observed values (as File 17a), so the
  # same in every imputed dataset: categories with fewer than 20 people go
  # into "Other (rare categories)", or into the largest category if together
  # they still have fewer than 20.
  pooling <- data.frame()
  for (v in features) {
    if (!is.factor(cohort[[v]])) next
    counts <- table(cohort[[v]])
    rare <- names(counts)[counts < minimum_category_size]
    new_label <- names(which.max(counts))
    if (sum(counts[rare]) >= minimum_category_size) new_label <- "Other (rare categories)"
    if (length(rare) > 0) {
      pooling <- rbind(pooling, data.frame(variable = v, category = rare, new_label = new_label))
    }
  }

  coefficients <- data.frame()
  tuning <- data.frame()
  for (i in seq_len(imputations)) {
    features_i <- completed[[i]][, features]

    # Rare categories pooled.
    for (k in seq_len(nrow(pooling))) {
      values <- as.character(features_i[[pooling$variable[k]]])
      values[values == pooling$category[k]] <- pooling$new_label[k]
      features_i[[pooling$variable[k]]] <- factor(values)
    }

    # A categorical feature left with a single category carries no
    # information and is left out.
    for (v in names(features_i)) {
      if (is.factor(features_i[[v]]) && nlevels(droplevels(features_i[[v]])) < 2) features_i[[v]] <- NULL
    }

    # Model matrix. model.matrix records, for each column, the number of the
    # feature it came from ("assign"), which gives each column's feature.
    X <- model.matrix(~ ., data = features_i)
    column_feature <- names(features_i)[attr(X, "assign")[-1]]
    X <- X[, -1, drop = FALSE]

    # Elastic net tuned by grouped cross-validation (whole-school folds):
    # one fit per alpha, then the alpha with the lowest cross-validated
    # deviance (the first, if tied).
    cv_fits <- list()
    deviances <- c()
    for (a in seq_along(alpha_grid)) {
      cv_fits[[a]] <- cv.glmnet(X, outcome_values, family = "binomial", alpha = alpha_grid[a],
                                foldid = tuning_folds, type.measure = "deviance")
      deviances[a] <- min(cv_fits[[a]]$cvm)
    }
    best <- list(alpha = alpha_grid[which.min(deviances)], fit = cv_fits[[which.min(deviances)]])
    beta <- as.matrix(coef(best$fit, s = "lambda.min"))[, 1]   # intercept first

    # Each term's coefficient, with the column's SD and, for 0/1 columns,
    # the people with the term and the events and non-events among them
    # (for the disclosure check in File 19b).
    terms_i <- data.frame(imputation = i, term = names(beta),
                          feature = c("(Intercept)", column_feature), coefficient = unname(beta),
                          column_sd = NA, binary = FALSE, people_with_term = NA,
                          events = NA, non_events = NA)
    for (j in seq_len(ncol(X))) {
      has_term <- X[, j] == 1
      terms_i$column_sd[j + 1] <- sd(X[, j])
      terms_i$binary[j + 1] <- all(X[, j] %in% c(0, 1))
      terms_i$people_with_term[j + 1] <- sum(has_term)
      terms_i$events[j + 1] <- sum(has_term & outcome_values == 1)
      terms_i$non_events[j + 1] <- sum(has_term & outcome_values == 0)
    }
    coefficients <- rbind(coefficients, terms_i)

    lambda_path <- best$fit$lambda
    tuning <- rbind(tuning, data.frame(
      imputation = i, columns = ncol(X), chosen_alpha = best$alpha,
      lambda = best$fit$lambda.min,
      lambda_at_largest = best$fit$lambda.min == max(lambda_path),
      lambda_at_smallest = best$fit$lambda.min == min(lambda_path),
      nonzero_coefficients = sum(beta[-1] != 0)))
  }

  # i. Save the task.
  task_result <- list(
    task = t, landmark = landmark,
    coefficients = coefficients, tuning = tuning,
    screen = screen, switches = switches, convergence = convergence,
    logged_events = if (is.null(imputation$loggedEvents)) 0L else nrow(imputation$loggedEvents),
    events = sum(outcome_values == 1), non_events = sum(outcome_values == 0),
    imputation_minutes = imputation_minutes,
    total_minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins")),
    r_version = as.character(getRversion()),
    mice_version = as.character(packageVersion("mice")),
    glmnet_version = as.character(packageVersion("glmnet")),
    computer = gsub("[^A-Za-z0-9_-]", "_", Sys.info()[["nodename"]])
  )
  temporary_file <- paste0(task_file, ".partial")
  saveRDS(task_result, temporary_file)
  invisible(file.rename(temporary_file, task_file))
  cat("  Task", t, "finished in", round(task_result$total_minutes, 1), "minutes\n")
}


###############################################################################
# Block 7. Summary
#
# Task counts only. No coefficients or predictions.
###############################################################################

tasks_done <- sum(file.exists(task_table$file))
cat("\nFile 19a finished on this computer.\n")
cat("  Tasks in the table:        ", nrow(task_table), "\n")
cat("  Task files now present:    ", tasks_done, "\n")
if (tasks_done == nrow(task_table)) cat("  All tasks done. Next: File 19b.\n")
