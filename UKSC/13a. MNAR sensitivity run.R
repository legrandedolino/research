###############################################################################
# 13a. MNAR sensitivity run
#
# Purpose
#   The main imputations (File 8) assume the age-32 outcome is missing at
#   random (MAR) given the imputation model. That can't be checked from the
#   data. Following the CLS user guide, this script asks how the key results
#   would change if the outcome were missing not at random (MNAR): the
#   imputed values are modified according to MNAR scenarios, and the analysis
#   model is fitted again (a pattern-mixture approach; Carpenter and Kenward
#   2013). File 13b combines the results and finds the tipping points: how
#   far from MAR the data would have to be before the conclusions change
#   (CLS user guide: "At what value ... is there no longer evidence of the
#   association of interest? How plausible is such a value in practice?").
#
# The scenarios
#   The outcome is imputed in two parts: whether a person is in paid work at
#   32, and, among those in paid work, whether the job is SRET. The scenarios
#   concern the second part. Scenario delta: among people in paid work, those
#   whose occupation is missing have log-odds of an SRET job that differ by
#   delta from otherwise similar people whose occupation was observed (an
#   odds ratio of exp(delta)). delta = 0 is MAR.
#   Five families, each over delta = -2 to +2 in steps of 0.5 (odds ratios
#   0.14 to 7.4):
#     uniform       delta applies to everyone with a missing outcome;
#     liking        only to those who liked science "a lot" at 14;
#     competence    only to those who rated themselves "very good" at science;
#     favourite     only to those who named science as a favourite subject;
#     maths         only to those who liked mathematics "a lot".
#   A uniform shift mainly moves the prevalence. An association can only be
#   challenged by a shift that differs with the exposure, which is what the
#   other four families do (the CLS user guide notes that the shift may
#   differ between subgroups). delta = 0 is shared by all families.
#
# How each scenario is applied, for each of the main imputations
#   1. The equation for SRET among the employed is fitted again to the people
#      whose occupation was observed, with the same predictors as the main
#      imputation model and the same method: logistic regression on a
#      bootstrap sample (no pseudo-observations are added), which carries the
#      uncertainty in its coefficients into the imputations.
#   2. For everyone whose outcome was missing, delta (or 0 outside the
#      scenario's group) is added to the linear predictor, and a new value is
#      drawn. The overall outcome is then recalculated from employment and SRET
#      among the employed, as in the main model. Everything else in the
#      completed dataset is left as imputed.
#   3. The analysis is fitted as in File 12's primary analysis: the whole
#      target population, the initial design weight, PSU and the original
#      strata.
#   Common random numbers: in every scenario, imputation i uses the same seed
#   (base seed + i), so the bootstrap sample and the random draws are the
#   same and only delta changes. The results then change smoothly with
#   delta, rather than jumping with Monte Carlo noise.
#
# Key estimands (the primary analysis only)
#   Weighted SRET prevalence, overall and among people in paid work; the
#   standardized risk differences of the three science contrasts in the
#   baseline model (File 5); the odds ratio for liking mathematics (per
#   point) in the baseline model.
#
# Survey weights and the imputation model (a deliberate departure)
#   Silverwood et al. (2024) weighted their imputation models by the initial
#   design weights. As in the main imputation (File 7), that is not done
#   here: the refitted equation is unweighted, with the log design weight and
#   the sampling stratum among its predictors, and the survey design is
#   applied in the analysis. Weight interactions are not included, which may
#   overstate standard errors (Silverwood et al. 2024).
#
# Running on several computers or R sessions
#   Each task is one scenario (all imputations). Give each computer or R
#   session its own tasks_to_run (Block 1); files that already exist are
#   skipped. Copy every task file into tasks_dir/13a on one computer and run
#   File 13b there. Results don't depend on which computer ran which task.
#
# Reads
#   data_dir/07_imputation_spec.rds
#   data_dir/08_imputations_main.rds
#   data_dir/05_complete_case_results.rds   (the baseline model and contrasts)
#
# Writes
#   tasks_dir/13a/scenario_NNN.rds   one per scenario: the estimates, their
#       variances and, for the disclosure checks, event counts, for every
#       imputation (no participant-level data), plus run time and versions
#
# References
#   Carpenter JR, Kenward MG (2013). Multiple Imputation and its Application.
#   Wiley.
#   Silverwood R, et al. (2026). Handling missing data in the CLS cohort
#   studies: User guide (Version 2). UCL Centre for Longitudinal Studies.
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

paths     <- readRDS("project_paths.rds")
data_dir  <- paths$data_dir     # datasets passed between scripts
tasks_dir <- paths$tasks_dir    # one file per task in split runs

# Which tasks (scenarios) this computer or R session runs: "all", or task
# numbers, e.g. 1:14.
tasks_to_run <- "all"

base_seed   <- 13000                   # imputation i uses seed base_seed + i
delta_grid  <- seq(-2, 2, by = 0.5)    # log odds ratio scale

task_folder <- file.path(tasks_dir, "13a")
dir.create(task_folder, recursive = TRUE, showWarnings = FALSE)

# As in Files 5 and 12: a stratum with a single PSU in a subset is measured
# against the overall mean.
options(survey.lonely.psu = "adjust")


###############################################################################
# Block 2. Read the specification, the main imputations and the model
###############################################################################

imputation_spec <- readRDS(file.path(data_dir, "07_imputation_spec.rds"))
imputation      <- readRDS(file.path(data_dir, "08_imputations_main.rds"))
file5           <- readRDS(file.path(data_dir, "05_complete_case_results.rds"))

baseline_model <- file5$models[file5$models$model == "M0 baseline", ]
contrasts      <- file5$contrasts
design_data    <- imputation_spec$design_data


###############################################################################
# Block 3. Scenario table
#
# One row per scenario: its family, delta, and which people with a missing
# outcome receive delta. The uniform family includes delta = 0 (MAR); the
# other families use the non-zero values only.
###############################################################################

scenarios <- rbind(
  data.frame(family = "uniform", delta = delta_grid),
  data.frame(family = "liking",     delta = delta_grid[delta_grid != 0]),
  data.frame(family = "competence", delta = delta_grid[delta_grid != 0]),
  data.frame(family = "favourite",  delta = delta_grid[delta_grid != 0]),
  data.frame(family = "maths",      delta = delta_grid[delta_grid != 0])
)
scenarios$task <- seq_len(nrow(scenarios))
scenarios$file <- file.path(task_folder, sprintf("scenario_%03d.rds", scenarios$task))

if (identical(tasks_to_run, "all")) tasks_to_run <- scenarios$task
if (!all(tasks_to_run %in% scenarios$task)) {
  stop("tasks_to_run includes task numbers outside 1-", nrow(scenarios), ".")
}


###############################################################################
# Block 4. The equation for SRET among the employed
#
# Its predictors are those of the main imputation model (the predictor
# matrix), and its method must be the bootstrap logistic regression.
# "Missing" rows are those the main model imputed (the where matrix);
# "observed" rows are those with an observed occupation.
###############################################################################

if (imputation_spec$method[["sret_among_employed32"]] != "logreg.boot") {
  stop("The main model doesn't impute SRET among the employed with ",
       "logreg.boot; this script assumes it does.")
}
sret_predictors <- names(which(imputation_spec$predictor_matrix["sret_among_employed32", ] == 1))
sret_right_side <- as.formula(paste("~", paste(sret_predictors, collapse = " + ")))

observed_rows <- which(!is.na(imputation_spec$data$sret_among_employed32))
missing_rows  <- which(imputation_spec$where[, "sret_among_employed32"])
observed_y    <- as.numeric(imputation_spec$data$sret_among_employed32[observed_rows] ==
                              "SRET professional")
sret_levels   <- levels(imputation_spec$data$sret_among_employed32)

baseline_formula <- as.formula(paste("sret_professional_32 ~", baseline_model$predictors))


###############################################################################
# Block 5. Run the scenarios
#
# For each scenario, and each main imputation i (seed base_seed + i):
#   a. the completed dataset;
#   b. the SRET-among-employed equation, refitted on a bootstrap sample of
#      the observed rows (a new bootstrap sample is drawn if a category is
#      missing from it, so that every coefficient can be estimated);
#   c. new values for the missing rows, with delta added for the scenario's
#      group, and the overall outcome recalculated;
#   d. the analysis (File 12's primary analysis) and the key estimands, each
#      with its variance, the design's degrees of freedom and event counts.
###############################################################################

cat("\nMNAR sensitivity run:", length(tasks_to_run), "of", nrow(scenarios),
    "scenarios on this computer;", imputation$m, "imputations each.\n")

for (t in tasks_to_run) {
  task_file <- scenarios$file[t]
  if (file.exists(task_file)) {
    cat("  Scenario", t, "already done; skipped.\n")
    next
  }
  family <- scenarios$family[t]
  delta  <- scenarios$delta[t]
  cat("  Scenario", t, "(", family, ", delta", delta, ") started",
      format(Sys.time(), "%H:%M"), "\n")
  start_time <- Sys.time()

  results <- data.frame()
  for (i in seq_len(imputation$m)) {
    set.seed(base_seed + i)                     # the same in every scenario

    # a. The completed dataset.
    data_i <- complete(imputation, i)

    # b. Refit the equation on a bootstrap sample of the observed rows.
    X <- model.matrix(sret_right_side, model.frame(sret_right_side, data_i, na.action = na.pass))
    if (anyNA(X)) stop("A predictor of SRET among the employed is missing in imputation ", i, ".")
    beta <- NA
    attempts <- 0
    while (anyNA(beta)) {
      attempts <- attempts + 1
      if (attempts > 20) stop("The SRET equation couldn't be estimated in imputation ", i, ".")
      boot <- sample(length(observed_rows), replace = TRUE)
      fit  <- glm.fit(X[observed_rows[boot], , drop = FALSE], observed_y[boot],
                      family = binomial())
      beta <- fit$coefficients
    }

    # c. New values for the missing rows, with delta added to the log-odds
    # for the scenario's group: everyone with a missing value (uniform), or
    # only those at the high level of one exposure.
    linear_predictor <- as.vector(X[missing_rows, , drop = FALSE] %*% beta)
    random_draw      <- runif(length(missing_rows))
    in_group <- rep(TRUE, length(missing_rows))                     # uniform
    if (family == "liking")     in_group <- data_i$science_liking14[missing_rows] == 4
    if (family == "competence") in_group <- data_i$science_competence14[missing_rows] == 4
    if (family == "favourite")  in_group <- data_i$science_favourite14[missing_rows] == "Science is a favourite subject"
    if (family == "maths")      in_group <- data_i$math_liking14[missing_rows] == 4
    shift <- rep(0, length(missing_rows))
    shift[in_group] <- delta
    new_sret <- random_draw < plogis(linear_predictor + shift)
    new_values <- rep("Other occupation", length(missing_rows))
    new_values[new_sret] <- "SRET professional"
    data_i$sret_among_employed32[missing_rows] <- factor(new_values, levels = sret_levels)

    # The overall outcome recalculated: SRET among the employed if in paid
    # work, otherwise 0.
    in_paid_work <- data_i$employed32 == "In paid work"
    data_i$sret_professional_32 <- as.numeric(data_i$sret_among_employed32 == "SRET professional")
    data_i$sret_professional_32[!in_paid_work] <- 0
    if (anyNA(data_i$sret_professional_32)) {
      stop("The recalculated outcome has missing values in imputation ", i, ".")
    }

    # d. The primary analysis (as File 12).
    analysis_data <- cbind(design_data, data_i)
    analysis_data$sret_among_employed32 <-
      as.numeric(analysis_data$sret_among_employed32 == "SRET professional")
    design <- svydesign(ids = ~psu, strata = ~stratum, weights = ~design_weight,
                        nest = TRUE, data = analysis_data)
    design_employed <- subset(design, employed32 == "In paid work")

    y_all      <- design$variables$sret_professional_32
    y_employed <- design_employed$variables$sret_among_employed32
    prevalence          <- svymean(~sret_professional_32, design)
    prevalence_employed <- svymean(~sret_among_employed32, design_employed)

    fit <- svyglm(baseline_formula, design = design, family = quasibinomial(),
                  control = glm.control(maxit = 100))
    w   <- weights(design)

    estimates <- data.frame(
      estimand = c("Prevalence: SRET overall", "Prevalence: SRET among the employed",
                   "Odds ratio: liking mathematics (per point)"),
      scale    = c("proportion", "proportion", "log odds ratio"),
      exposure = c(NA, NA, "math_liking14"),
      estimate = c(coef(prevalence)[[1]], coef(prevalence_employed)[[1]],
                   coef(fit)[["math_liking14"]]),
      variance = c(SE(prevalence)[[1]]^2, SE(prevalence_employed)[[1]]^2,
                   vcov(fit)["math_liking14", "math_liking14"]),
      df_complete = c(degf(design), degf(design_employed), degf(design)),
      events      = c(sum(y_all == 1), sum(y_employed == 1), sum(y_all == 1)),
      non_events  = c(sum(y_all == 0), sum(y_employed == 0), sum(y_all == 0))
    )

    # Standardized risk differences of the science contrasts (as Files 5, 12).
    for (k in seq_len(nrow(contrasts))) {
      variable  <- contrasts$variable[k]
      data_high <- design$variables
      data_low  <- design$variables
      numeric_item <- is.numeric(data_high[[variable]])
      if (numeric_item)  data_high[[variable]] <- as.numeric(contrasts$high[k])
      if (numeric_item)  data_low[[variable]]  <- as.numeric(contrasts$low[k])
      if (!numeric_item) data_high[[variable]] <- factor(contrasts$high[k], levels = levels(data_high[[variable]]))
      if (!numeric_item) data_low[[variable]]  <- factor(contrasts$low[k],  levels = levels(data_low[[variable]]))
      predictor_terms <- delete.response(terms(fit))
      X_high <- model.matrix(predictor_terms, data = data_high, xlev = fit$xlevels)
      X_low  <- model.matrix(predictor_terms, data = data_low,  xlev = fit$xlevels)
      p_high <- as.vector(plogis(X_high %*% coef(fit)))
      p_low  <- as.vector(plogis(X_low  %*% coef(fit)))
      rd <- sum(w * (p_high - p_low)) / sum(w)
      gradient <- colSums(w * (p_high * (1 - p_high) * X_high -
                               p_low  * (1 - p_low)  * X_low)) / sum(w)
      estimates <- rbind(estimates, data.frame(
        estimand = paste("Risk difference:", contrasts$contrast[k]),
        scale = "risk difference", exposure = variable, estimate = rd,
        variance = as.numeric(t(gradient) %*% vcov(fit) %*% gradient),
        df_complete = degf(design),
        events = sum(y_all == 1), non_events = sum(y_all == 0)
      ))
    }

    estimates$imputation <- i
    results <- rbind(results, estimates)
  }

  task_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  task_result <- list(
    task = t, family = family, delta = delta,
    results      = results,
    minutes      = task_minutes,
    r_version    = as.character(getRversion()),
    mice_version = as.character(packageVersion("mice")),
    computer     = gsub("[^A-Za-z0-9_-]", "_", Sys.info()[["nodename"]])
  )
  temporary_file <- paste0(task_file, ".partial")
  saveRDS(task_result, temporary_file)
  file.rename(temporary_file, task_file)
  cat("  Scenario", t, "finished in", round(task_minutes, 1), "minutes\n")
}


###############################################################################
# Block 6. Summary
#
# Task counts only. No data values.
###############################################################################

tasks_done <- sum(file.exists(scenarios$file))
cat("\nFile 13a finished on this computer.\n")
cat("  Scenarios in the table:     ", nrow(scenarios), "\n")
cat("  Scenario files now present: ", tasks_done, "\n")
if (tasks_done == nrow(scenarios)) {
  cat("  All scenarios done. Next: File 13b.\n")
}
