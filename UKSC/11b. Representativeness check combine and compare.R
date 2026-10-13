###############################################################################
# 11b. Representativeness check combine and compare
#
# Purpose
#   Combine the imputations of File 11a (the age-14 analysis variables
#   masked among Wave-9 non-respondents and imputed again), check them, and
#   compare the distributions of the masked variables estimated four ways,
#   following Silverwood et al. (2024):
#     truth                    all known values, design weight;
#     complete cases           Wave-9 respondents' values only, design
#                              weight: the bias from attrition;
#     complete cases, W9FINWT  Wave-9 respondents' values with the Wave-9
#                              non-response weight: how far weighting alone
#                              restores representativeness (the approach of
#                              the prediction analysis);
#     imputation               respondents' values plus the imputed values
#                              of the masked non-respondents, design weight,
#                              combined with Rubin's rules.
#   All four are estimated among the people whose value is actually known.
#   If the imputation estimate is close to the truth while the complete-case
#   estimate is not, the imputation model restores representativeness for
#   these variables. This does not prove that the age-32 data are missing at
#   random. Silverwood et al. could not restore the distribution of sex;
#   whether that recurs is reported.
#
# Survey weights and the imputation model (a deliberate departure)
#   Silverwood et al. (2024) weighted their imputation models by the initial
#   design weights. That is not done here: R's mice package has no weighted
#   imputation methods for categorical variables. The log of the initial
#   design weight and the sampling stratum are predictors in every
#   imputation equation instead (File 7), the imputation models are
#   unweighted, and survey weights, PSUs and strata are applied when each
#   completed dataset is analysed (here). Weight interactions are not
#   included, which may overstate standard errors (Silverwood et al. 2024).
#
# Checks before combining (the script stops if any fails; as File 8b)
#   1. Every task in the table has a file.
#   2. Each file is the task it claims to be (task number, seed, one
#      imputation, the specification's iterations).
#   3. Each imputation used the masked specification saved by File 11a,
#      unchanged.
#   4. Only structural (not applicable) values are still missing.
#   5. Every task was run with the same versions of R and mice.
#   6. The combined object gives back every task's completed dataset
#      exactly (after repairing a quirk of mice::ibind; Block 6).
#
# Diagnostics (CLS user guide; as File 10, for the masked variables)
#   Trace plots and the convergence table for the masked variables.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   Percentages and means only; no counts. A category's row is hidden
#   unless at least 10 people have that value in each group it describes:
#   everyone with a known value, the Wave-9 respondents (with and without a
#   Wave-9 weight), and the non-respondents whose value was masked.
#   Secondary suppression within each variable (percentages add to 100).
#   Trace plots and convergence rows only for variables with at least 10
#   imputed values.
#
# Reads
#   data_dir/07_imputation_spec.rds   (the true values)
#   data_dir/11_masked_spec.rds
#   tasks_dir/11a/imputation_NNN.rds
#
# Writes
#   data_dir/11_masked_imputations.rds      combined imputations (participant
#                                           level)
#   data_dir/11_representativeness.rds      exact comparison results
#   output_dir/private_11_*                 exact tables and all-variable plots
#   review_dir/11_task_log.csv, 11_logged_events.csv, 11_warnings.csv
#   review_dir/11_trace_plots.pdf, 11_convergence.csv
#   review_dir/11_representativeness.csv, 11_imputation_number_check.csv,
#   review_dir/11_disclosure_log.csv
#
# References
#   Barnard J, Rubin DB (1999). Biometrika 86(4): 948-955.
#   Rubin DB (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
#   Silverwood RJ, et al. (2024). Longitudinal and Life Course Studies 15(2):
#   227-250.
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
tasks_dir  <- paths$tasks_dir    # one file per task in split runs

disclosure_threshold <- 10       # smallest count that may be shown
minimum_m            <- 50       # CLS user guide: at least 50 imputations
base_seed            <- 11000    # as in File 11a: task t used seed base_seed + t
task_folder          <- file.path(tasks_dir, "11a")

# As in Files 5 and 12: a stratum with a single PSU in a subset is measured
# against the overall mean.
options(survey.lonely.psu = "adjust")


###############################################################################
# Block 2. Read the specifications and build the task table
###############################################################################

spec_file <- file.path(data_dir, "11_masked_spec.rds")
if (!file.exists(spec_file)) {
  stop("11_masked_spec.rds is missing. Run File 11a on this computer first ",
       "(it saves the specification before running any task).")
}
masked_spec     <- readRDS(spec_file)
imputation_spec <- readRDS(file.path(data_dir, "07_imputation_spec.rds"))   # true values
masked_variables <- masked_spec$masked_variables

task_table <- data.frame(task = seq_len(masked_spec$production_m))
task_table$seed <- base_seed + task_table$task
task_table$file <- file.path(task_folder,
                             sprintf("imputation_%03d.rds", task_table$task))

# Which values should still be missing after imputation: the structural
# (not applicable) values.
expected_missing <- is.na(masked_spec$data) & !masked_spec$where


###############################################################################
# Block 3. Check 1: every task file is present
###############################################################################

missing_tasks <- task_table$task[!file.exists(task_table$file)]
if (length(missing_tasks) > 0) {
  stop(length(missing_tasks), " task files are missing (tasks ",
       paste(missing_tasks, collapse = ", "), "). Run them with File 11a, ",
       "or copy their files into ", task_folder, ".")
}


###############################################################################
# Block 4. Checks 2-4, and combining, one task file at a time
#
# Each file is read and checked, then added to the combined imputations with
# mice::ibind, in task order. A failed check stops the script and names the
# task, so that task's file can be deleted and rerun in File 11a.
###############################################################################

combined      <- NULL
task_log      <- data.frame()
logged_events <- data.frame()
all_warnings  <- character()

for (t in task_table$task) {
  result <- readRDS(task_table$file[t])
  imputation <- result$imputation

  # 2. The right task, seed, number of imputations and iterations.
  if (!identical(as.integer(result$task), as.integer(t)) ||
      !identical(as.numeric(result$seed), as.numeric(task_table$seed[t])) ||
      !is.mids(imputation) || imputation$m != 1 ||
      imputation$iteration != masked_spec$iterations) {
    stop("Task ", t, "'s file is not the expected task (task number, seed, ",
         "one imputation or iterations differ). Delete it and rerun the task.")
  }

  # 3. The masked specification, unchanged. mice stores the predictor
  # matrix with a different number type and drops the row names of the
  # where matrix, so those two are compared by value.
  if (!identical(imputation$data, masked_spec$data) ||
      !identical(imputation$method, masked_spec$method) ||
      !isTRUE(all.equal(imputation$predictorMatrix, masked_spec$predictor_matrix)) ||
      !identical(colnames(imputation$where), colnames(masked_spec$where)) ||
      any(imputation$where != masked_spec$where) ||
      !identical(imputation$visitSequence, masked_spec$visit_sequence)) {
    stop("Task ", t, " was not run with the current masked specification. ",
         "Delete its file and rerun the task with the current specification.")
  }

  # 4. Only structural values still missing.
  completed <- complete(imputation, 1)
  if (any(is.na(completed) != expected_missing)) {
    stop("Task ", t, "'s completed data has missing values other than the ",
         "structural (not applicable) ones. Delete its file and rerun the task.")
  }

  # Add the imputation to the combined set.
  if (is.null(combined)) {
    combined <- imputation
  } else {
    combined <- ibind(combined, imputation)
  }

  # Record the task (no data values).
  events <- imputation$loggedEvents
  n_events <- if (is.null(events)) 0L else nrow(events)
  if (n_events > 0) {
    events$task <- t
    logged_events <- rbind(logged_events, events)
  }
  all_warnings <- c(all_warnings, result$warnings)

  task_log <- rbind(task_log, data.frame(
    task          = t,
    seed          = result$seed,
    computer      = result$computer,
    minutes       = round(result$minutes, 1),
    warnings      = length(result$warnings),
    logged_events = n_events,
    r_version     = result$r_version,
    mice_version  = result$mice_version
  ))
}


###############################################################################
# Block 5. Check 5: the same R and mice versions for every task
###############################################################################

if (length(unique(task_log$r_version)) > 1 ||
    length(unique(task_log$mice_version)) > 1) {
  stop("Tasks were run with different versions of R (",
       paste(unique(task_log$r_version), collapse = ", "), ") or mice (",
       paste(unique(task_log$mice_version), collapse = ", "), "). Align the ",
       "versions (File 0 records them on each computer) and rerun the tasks ",
       "that used the other versions.")
}


###############################################################################
# Block 6. Repair the combined object, check it, and save it
#
# In mice 3.19, ibind drops the empty entries of variables with nothing
# imputed, and mice::complete() then matches the remaining entries to the
# wrong variables by position. The full list is restored in the data's
# column order, and check 6 confirms that the combined object gives back
# every task's completed dataset exactly.
###############################################################################

if (combined$m != nrow(task_table)) {
  stop("The combined imputations number ", combined$m, ", not ",
       nrow(task_table), ".")
}

full_imp <- vector("list", ncol(combined$data))
names(full_imp) <- names(combined$data)
for (variable in names(combined$imp)) {
  full_imp[[variable]] <- combined$imp[[variable]]
}
combined$imp <- full_imp

# 6. Each completed dataset from the combined object equals that task's own.
for (t in task_table$task) {
  task_imputation <- readRDS(task_table$file[t])$imputation
  if (!identical(complete(combined, t), complete(task_imputation, 1))) {
    stop("Completed dataset ", t, " from the combined imputations differs ",
         "from task ", t, "'s own. The combined object is not usable.")
  }
}

saveRDS(combined, file.path(data_dir, "11_masked_imputations.rds"))


###############################################################################
# Block 7. Task log, logged events and warnings
#
# Disclosure check: task numbers, seeds, computer names, run times, version
# numbers, counts of warnings and logged events, variable and method names,
# and warning messages. None holds a count of people or a data value.
###############################################################################

task_log$disclosure_checked <- TRUE
write.csv(task_log, file.path(review_dir, "11_task_log.csv"), row.names = FALSE)

if (nrow(logged_events) == 0) {
  logged_events <- data.frame(it = integer(), im = integer(), dep = character(),
                              meth = character(), out = character(),
                              task = integer())
}
logged_events$disclosure_checked <- rep(TRUE, nrow(logged_events))
write.csv(logged_events, file.path(review_dir, "11_logged_events.csv"),
          row.names = FALSE)

warning_table <- data.frame(message = unique(all_warnings))
warning_table$times <- as.integer(table(all_warnings)[warning_table$message])
warning_table$disclosure_checked <- rep(TRUE, nrow(warning_table))
write.csv(warning_table, file.path(review_dir, "11_warnings.csv"), row.names = FALSE)


###############################################################################
# Block 8. Diagnostics for the masked variables
#
# Trace plots (each imputation a thin grey line, their average a black
# line) and the convergence table (split R-hat, trend and drift_vs_spread
# over the second half of the iterations), as in File 10. Only the masked
# variables are shown; the rest of the model is the main model, diagnosed in
# File 10. Review copies include only variables with at least 10 imputed
# values.
###############################################################################

n_imputed <- colSums(masked_spec$where)
traced    <- masked_variables
reviewed  <- traced[n_imputed[traced] >= disclosure_threshold]

for (copy in c("private", "review")) {
  if (copy == "private") {
    pdf(file.path(output_dir, "private_11_trace_plots_all.pdf"), width = 8.5, height = 11)
    plot_variables <- traced
  } else {
    pdf(file.path(review_dir, "11_trace_plots.pdf"), width = 8.5, height = 11)
    plot_variables <- reviewed
  }
  par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
  for (variable in plot_variables) {
    chain_means <- combined$chainMean[variable, , ]        # iterations x imputations
    chain_sds   <- sqrt(combined$chainVar[variable, , ])
    matplot(chain_means, type = "l", lty = 1, col = grey(0.6, alpha = 0.4),
            xlab = "Iteration", ylab = "Mean of imputed values",
            main = variable, cex.main = 0.9)
    lines(rowMeans(chain_means), lwd = 2)
    matplot(chain_sds, type = "l", lty = 1, col = grey(0.6, alpha = 0.4),
            xlab = "Iteration", ylab = "SD of imputed values",
            main = variable, cex.main = 0.9)
    lines(rowMeans(chain_sds), lwd = 2)
  }
  dev.off()
}

iterations  <- combined$iteration
m           <- combined$m
second_half <- (floor(iterations / 2) + 1):iterations        # e.g. 16-30
half_length <- floor(length(second_half) / 2)
early_part  <- head(second_half, half_length)
late_part   <- tail(second_half, half_length)

convergence <- data.frame()
for (variable in reviewed) {
  for (statistic in c("mean", "sd")) {
    if (statistic == "mean") {
      values <- combined$chainMean[variable, , ]            # iterations x imputations
    } else {
      values <- sqrt(combined$chainVar[variable, , ])
    }

    # Split R-hat (needs at least 2 iterations in each half).
    split_rhat <- NA
    if (half_length >= 2) {
      halves  <- cbind(values[early_part, , drop = FALSE],
                       values[late_part, , drop = FALSE])
      n       <- nrow(halves)
      within  <- mean(apply(halves, 2, var))
      between <- n * var(colMeans(halves))
      if (within > 0) split_rhat <- sqrt(((n - 1) / n * within + between / n) / within)
    }

    # Trend over the second half, and drift relative to the spread between
    # imputations.
    trend_data <- data.frame(
      value      = as.vector(values[second_half, , drop = FALSE]),
      iteration  = rep(second_half, times = m),
      imputation = factor(rep(seq_len(m), each = length(second_half)))
    )
    trend <- summary(lm(value ~ iteration + imputation, data = trend_data))$coefficients
    trend_slope <- trend["iteration", "Estimate"]
    trend_p     <- trend["iteration", "Pr(>|t|)"]
    spread <- mean(apply(values[second_half, , drop = FALSE], 1, sd))
    drift_vs_spread <- if (spread > 0) abs(trend_slope) * (length(second_half) - 1) / spread else NA

    convergence <- rbind(convergence, data.frame(
      variable = variable, statistic = statistic, imputations = m,
      split_rhat = round(split_rhat, 3),
      trend_per_iteration = signif(trend_slope, 3),
      trend_p = signif(trend_p, 3),
      drift_vs_spread = round(drift_vs_spread, 2),
      flag = isTRUE(split_rhat > 1.1) ||
             (isTRUE(trend_p < 0.01) && isTRUE(drift_vs_spread > 1))
    ))
  }
}

# Disclosure check: summaries of the means and SDs of imputed values, only
# for variables with at least 10 imputed values; no counts of people.
convergence$disclosure_checked <- rep(TRUE, nrow(convergence))
write.csv(convergence, file.path(review_dir, "11_convergence.csv"), row.names = FALSE)


###############################################################################
# Block 9. The four estimates of each masked variable's distribution
#
# For each masked variable, among the people whose value is known ("known"):
#   truth, and the two complete-case estimates, from the true values;
#   imputation: each completed dataset (true values for respondents,
#     imputed values for masked non-respondents) gives an estimate and its
#     variance, combined with Rubin's rules (Barnard-Rubin degrees of
#     freedom, as in File 12).
# Categorical variables: the weighted percentage in each category. Numeric
# variables: the weighted mean. (For a single mean, survey returns the
# standard error as a 1 x 1 matrix, so standard errors are converted to
# plain numbers.)
# The counts behind each row are kept for the disclosure check.
###############################################################################

true_data <- cbind(imputation_spec$design_data, imputation_spec$data[masked_variables])
responded <- true_data$responded_w9 %in% TRUE
for (variable in masked_variables) {
  true_data[[paste0("known_", variable)]] <- !is.na(true_data[[variable]])
}

design_all <- svydesign(ids = ~psu, strata = ~stratum, weights = ~design_weight,
                        nest = TRUE, data = true_data)
respondent_rows <- responded & !is.na(true_data$w9_final_weight)
design_w9 <- svydesign(ids = ~psu, strata = ~stratum, weights = ~w9_final_weight,
                       nest = TRUE, data = true_data[respondent_rows, ])

# The completed datasets' values of the masked variables, one table each.
completed_values <- list()
for (i in seq_len(m)) {
  completed_values[[i]] <- complete(combined, i)[masked_variables]
}

comparison <- data.frame()
imputation_number <- data.frame()

for (variable in masked_variables) {
  known <- paste0("known_", variable)
  estimate_formula <- as.formula(paste("~", variable))
  is_categorical <- is.factor(true_data[[variable]])

  # Truth and complete cases (true values).
  truth    <- svymean(estimate_formula, subset(design_all, get(known)))
  cc       <- svymean(estimate_formula, subset(design_all, get(known) & responded_w9 %in% TRUE))
  cc_w9    <- svymean(estimate_formula, subset(design_w9, get(known)))
  truth_df <- degf(subset(design_all, get(known)))
  cc_df    <- degf(subset(design_all, get(known) & responded_w9 %in% TRUE))
  cc_w9_df <- degf(subset(design_w9, get(known)))

  # Imputation: one estimate per completed dataset.
  estimates <- NULL
  variances <- NULL
  for (i in seq_len(m)) {
    data_i <- true_data
    data_i[[variable]] <- completed_values[[i]][[variable]]
    design_i <- svydesign(ids = ~psu, strata = ~stratum, weights = ~design_weight,
                          nest = TRUE, data = data_i)
    estimate_i <- svymean(estimate_formula, subset(design_i, get(known)))
    estimates  <- cbind(estimates, coef(estimate_i))
    variances  <- cbind(variances, as.numeric(SE(estimate_i))^2)
  }

  # Rubin's rules, per category (or for the mean).
  q_bar   <- rowMeans(estimates)
  within  <- rowMeans(variances)
  between <- apply(estimates, 1, var)
  total   <- within + (1 + 1 / m) * between
  rvi     <- (1 + 1 / m) * between / within
  lambda  <- (1 + 1 / m) * between / total
  df_observed <- (truth_df + 1) / (truth_df + 3) * truth_df * (1 - lambda)
  # Barnard-Rubin degrees of freedom, for each category at once (as File
  # 12): with no missing information (lambda = 0) they are df_observed.
  df_old  <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  df      <- ifelse(is.finite(df_old), df_old * df_observed / (df_old + df_observed), df_observed)
  fmi     <- (rvi + 2 / (df + 3)) / (rvi + 1)
  mi_se   <- sqrt(total)
  mi_critical <- qt(0.975, df = df)

  # Categories (or "(mean)"), and the counts behind each row.
  if (is_categorical) {
    categories <- levels(true_data[[variable]])
    scale      <- 100                                      # percentages
    values_known <- true_data[[variable]]
    n_known      <- as.vector(table(factor(values_known, levels = categories)))
    n_respondent <- as.vector(table(factor(values_known[responded], levels = categories)))
    n_respondent_weighted <- as.vector(table(factor(values_known[respondent_rows], levels = categories)))
  } else {
    categories <- "(mean)"
    scale      <- 1
    n_known      <- sum(!is.na(true_data[[variable]]))
    n_respondent <- sum(!is.na(true_data[[variable]][responded]))
    n_respondent_weighted <- sum(!is.na(true_data[[variable]][respondent_rows]))
  }

  comparison <- rbind(comparison, data.frame(
    variable = variable, category = categories,
    scale = if (is_categorical) "percent" else "mean",
    truth = scale * coef(truth),
    truth_low  = scale * (coef(truth) - qt(0.975, truth_df) * as.numeric(SE(truth))),
    truth_high = scale * (coef(truth) + qt(0.975, truth_df) * as.numeric(SE(truth))),
    complete_cases = scale * coef(cc),
    complete_cases_low  = scale * (coef(cc) - qt(0.975, cc_df) * as.numeric(SE(cc))),
    complete_cases_high = scale * (coef(cc) + qt(0.975, cc_df) * as.numeric(SE(cc))),
    complete_cases_w9finwt = scale * coef(cc_w9),
    complete_cases_w9finwt_low  = scale * (coef(cc_w9) - qt(0.975, cc_w9_df) * as.numeric(SE(cc_w9))),
    complete_cases_w9finwt_high = scale * (coef(cc_w9) + qt(0.975, cc_w9_df) * as.numeric(SE(cc_w9))),
    imputation = scale * q_bar,
    imputation_low  = scale * (q_bar - mi_critical * mi_se),
    imputation_high = scale * (q_bar + mi_critical * mi_se),
    imputation_fmi = fmi,
    # Differences from the truth (percentage points, or units of the mean).
    difference_complete_cases = scale * (coef(cc) - coef(truth)),
    difference_complete_cases_w9finwt = scale * (coef(cc_w9) - coef(truth)),
    difference_imputation = scale * (q_bar - coef(truth)),
    # For the disclosure check (not shown in the checked table).
    n_known = n_known, n_respondent = n_respondent,
    n_respondent_weighted = n_respondent_weighted,
    n_masked = n_known - n_respondent,
    row.names = NULL
  ))

  # (An empty category has no variance, so its FMI is undefined and skipped.)
  imputation_number <- rbind(imputation_number, data.frame(
    variable = variable, largest_fmi = max(fmi, na.rm = TRUE), imputations = m,
    imputations_needed = max(minimum_m, ceiling(100 * max(fmi, na.rm = TRUE))),
    smallest_n = min(n_known, n_respondent, n_known - n_respondent)
  ))
}
imputation_number$rule_met <- imputation_number$imputations >= imputation_number$imputations_needed

write.csv(comparison, file.path(output_dir, "private_11_representativeness.csv"),
          row.names = FALSE)
saveRDS(list(comparison = comparison, imputation_number = imputation_number),
        file.path(data_dir, "11_representativeness.rds"))


###############################################################################
# Block 10. Disclosure checks and review tables
###############################################################################

disclosure_log <- data.frame(table = character(), primary_suppressed = integer(),
                             secondary_suppressed = integer())

# --- Disclosure check: the comparison -------------------------------------------
# A row is hidden unless at least 10 people have that value (category) in
# each group it describes: everyone with a known value, the Wave-9
# respondents (with and without a Wave-9 weight), and the masked
# non-respondents. Rows with a zero count describe an empty category and are
# not hidden for that reason. Secondary: within each variable, if exactly
# one row is hidden, the row with the next smallest known count is hidden
# too (percentages add to 100). The counts are removed from the shared
# table.
checked <- comparison
counts  <- checked[, c("n_known", "n_respondent", "n_respondent_weighted", "n_masked")]
primary <- apply(counts > 0 & counts < disclosure_threshold, 1, any)
secondary <- rep(FALSE, nrow(checked))
for (v in unique(checked$variable[checked$scale == "percent"])) {
  rows <- which(checked$variable == v)
  if (sum(primary[rows]) == 1) {
    candidates <- rows[!primary[rows] & checked$n_known[rows] > 0]
    if (length(candidates) > 0) {
      secondary[candidates[which.min(checked$n_known[candidates])]] <- TRUE
    }
  }
}
value_columns <- setdiff(names(checked), c("variable", "category", "scale",
                                           "n_known", "n_respondent",
                                           "n_respondent_weighted", "n_masked"))
checked[primary | secondary, value_columns] <- NA
for (column in value_columns) {
  checked[[column]] <- round(checked[[column]], if (column == "imputation_fmi") 3 else 2)
}
checked <- checked[, c("variable", "category", "scale", value_columns)]
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "11_representativeness.csv"),
          row.names = FALSE, na = "")
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "11_representativeness", primary_suppressed = sum(primary),
  secondary_suppressed = sum(secondary)))
# -----------------------------------------------------------------------------

# --- Disclosure check: number of imputations -------------------------------------
# FMI summaries per variable; hidden if the variable has a group under 10.
checked <- imputation_number
hide <- checked$smallest_n < disclosure_threshold
checked[hide, c("largest_fmi", "imputations_needed", "rule_met")] <- NA
checked$largest_fmi <- round(checked$largest_fmi, 3)
checked <- checked[, c("variable", "largest_fmi", "imputations",
                       "imputations_needed", "rule_met")]
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "11_imputation_number_check.csv"),
          row.names = FALSE, na = "")
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "11_imputation_number_check", primary_suppressed = sum(hide),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------

disclosure_log <- rbind(disclosure_log, data.frame(
  table = c("11_task_log", "11_logged_events", "11_warnings",
            "11_trace_plots and 11_convergence"),
  primary_suppressed = c(0L, 0L, 0L, length(traced) - length(reviewed)),
  secondary_suppressed = 0L))
disclosure_log$threshold <- disclosure_threshold
write.csv(disclosure_log, file.path(review_dir, "11_disclosure_log.csv"), row.names = FALSE)


###############################################################################
# Block 11. Summary
#
# Counts and checks only. No estimates, which are read from the checked
# table.
###############################################################################

cat("\nFile 11b finished.\n")
cat("  Imputations combined:                 ", combined$m, "\n")
cat("  Checks (files, specification, missing values, versions,\n")
cat("    combined equals each task): all passed\n")
cat("  Logged events (all tasks):            ", nrow(logged_events), "\n")
cat("  Distinct warnings (all tasks):        ", nrow(warning_table), "\n")
cat("  Masked variables compared:            ", length(masked_variables), "\n")
cat("  Convergence flags (masked variables): ", sum(convergence$flag), "of",
    nrow(convergence), "rows\n")
cat("  Imputation-number rule (m >= 50 and >= largest FMI %):",
    if (all(imputation_number$rule_met)) "met" else "NOT met", "\n")
cat("  Review: 11_representativeness.csv, 11_imputation_number_check.csv,\n")
cat("    11_trace_plots.pdf, 11_convergence.csv, 11_task_log.csv\n")
