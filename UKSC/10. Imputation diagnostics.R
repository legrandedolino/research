###############################################################################
# 10. Imputation diagnostics
#
# Purpose
#   Check the main imputations (File 8) and the supplementary imputations
#   without the predictors of non-response (File 9) before any pooled
#   result is interpreted (CLS user guide):
#     a. routing and consistency: every completed dataset keeps the
#        observed values, leaves only structural (not applicable) values
#        missing, and has pathways and outcome equal to their parts;
#     b. convergence: trace plots of the mean and SD of the imputed values
#        at each iteration, one line per imputation (van Buuren 2018,
#        section 6.5.2), and a numeric convergence table;
#     c. observed versus imputed distributions, for both sets of
#        imputations side by side.
#   The test run in File 7 had 5 imputations. With every imputation of the
#   main run, the convergence statistics are much more precise.
#
# Reading the diagnostics
#   Trace plots: each imputation is a thin grey line and their average a
#   black line. The grey lines should overlap in one band without a trend
#   in the later iterations; the black line should be flat there. A black
#   line that is still moving means all the imputations drift together.
#   Convergence table: split R-hat (Vehtari et al. 2021) compares the
#   variation between the half-sequences of all imputations with the
#   variation within them; values close to 1 mean they agree. The trend and
#   drift_vs_spread measure how far the imputations move over the second
#   half of the iterations, relative to the spread between imputations.
#   Flags are a guide: closely related questions that are missing together
#   (e.g. the age-16 questions) mix slowly, which raises split R-hat
#   without a trend, and valid estimates can be reached before such
#   diagnostics settle (Oberman, van Buuren and Vink 2021).
#   Observed versus imputed: differences are expected under missing at
#   random when the people with missing values differ from the others, but
#   they should be explainable.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   - Trace plots and the convergence table in review_dir include only
#     variables with at least 10 imputed values; the private copies include
#     every variable.
#   - Observed-versus-imputed rows show percentages only. A category's
#     percentages are hidden if fewer than 10 people are observed in it or
#     fewer than 10 values (average per imputation) are imputed into it,
#     with secondary suppression within each variable and imputation set.
#   - Routing checks report counts of imputations passing, not people.
#
# Reads
#   data_dir/07_imputation_spec.rds
#   data_dir/08_imputations_main.rds
#   data_dir/09_supplementary_spec.rds
#   data_dir/09_imputations_supplementary.rds
#
# Writes
#   review_dir/10_routing_checks.csv
#   output_dir/private_10_trace_plots_main_all.pdf and
#     private_10_trace_plots_supplementary_all.pdf   every variable
#   review_dir/10_trace_plots_main.pdf and
#     10_trace_plots_supplementary.pdf              at least 10 imputed values
#   review_dir/10_convergence.csv
#   output_dir/private_10_observed_vs_imputed.csv
#   review_dir/10_observed_vs_imputed.csv
#   review_dir/10_disclosure_log.csv
###############################################################################

suppressPackageStartupMessages(library(mice))


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

disclosure_threshold <- 10    # smallest count that may be shown
rhat_flag            <- 1.1   # split R-hat above this is flagged
trend_p_flag         <- 0.01  # a trend is flagged if p is below this ...
drift_flag           <- 1     # ... and the drift exceeds the spread


###############################################################################
# Block 2. Read the two sets of imputations and their specifications
#
# Each set is kept with its specification (data before imputation, where
# matrix and methods) and a label. The loops below run over the two sets.
###############################################################################

imputation_sets <- list(
  main = list(
    label      = "Main (File 8)",
    spec       = readRDS(file.path(data_dir, "07_imputation_spec.rds")),
    imputation = readRDS(file.path(data_dir, "08_imputations_main.rds"))
  ),
  supplementary = list(
    label      = "Supplementary, without predictors of non-response (File 9)",
    spec       = readRDS(file.path(data_dir, "09_supplementary_spec.rds")),
    imputation = readRDS(file.path(data_dir, "09_imputations_supplementary.rds"))
  )
)

for (set in names(imputation_sets)) {
  spec <- imputation_sets[[set]]$spec
  imputation_sets[[set]]$n_imputed <- colSums(spec$where)      # values imputed per variable
  imputation_sets[[set]]$traced    <- spec$visit_sequence[
    spec$method[spec$visit_sequence] != "" &
    !startsWith(spec$method[spec$visit_sequence], "~")
  ]                                                             # variables imputed by a model
}


###############################################################################
# Block 3. Routing and consistency checks on every completed dataset
#
# For each completed dataset of each set:
#   a. observed values are unchanged;
#   b. only structural values are still missing;
#   c. each calculated variable equals its formula applied to its parts;
#   d. the outcome is 0 or 1 for everyone.
# The table counts the imputations passing each check. It is written first;
# the script then stops if any check failed.
###############################################################################

routing_checks <- data.frame()

for (set in names(imputation_sets)) {
  spec       <- imputation_sets[[set]]$spec
  imputation <- imputation_sets[[set]]$imputation
  passive    <- names(spec$method)[startsWith(spec$method, "~")]
  passed     <- matrix(FALSE, nrow = imputation$m, ncol = 4)

  for (i in seq_len(imputation$m)) {
    data_i <- complete(imputation, i)

    # a. Observed values unchanged.
    unchanged <- TRUE
    for (variable in names(spec$data)) {
      observed <- !is.na(spec$data[[variable]])
      changed  <- any(as.character(data_i[[variable]][observed]) !=
                        as.character(spec$data[[variable]][observed]))
      if (changed) unchanged <- FALSE
    }

    # b. Only structural values still missing.
    expected_missing <- is.na(spec$data) & !spec$where
    only_structural_missing <- !any(is.na(data_i) != expected_missing)

    # c. Calculated variables equal their formulas.
    calculated_consistent <- TRUE
    for (variable in passive) {
      rebuilt <- model.frame(as.formula(spec$method[[variable]]), data_i,
                             na.action = na.pass)[[1]]
      differs <- anyNA(rebuilt) || any(as.character(rebuilt) != as.character(data_i[[variable]]))
      if (differs) calculated_consistent <- FALSE
    }

    # d. Outcome is 0 or 1.
    outcome_binary <- all(data_i$sret_professional_32 %in% c(0, 1))

    passed[i, ] <- c(unchanged, only_structural_missing, calculated_consistent,
                     outcome_binary)
  }

  routing_checks <- rbind(routing_checks, data.frame(
    imputation_set = imputation_sets[[set]]$label,
    check = c("Observed values unchanged",
              "Only structural (not applicable) values still missing",
              "Pathways and outcome equal their parts",
              "Outcome is 0 or 1 for everyone"),
    imputations         = imputation$m,
    imputations_passing = colSums(passed)
  ))
}

# Disclosure check: counts of imputations, not of people; no data values.
routing_checks$disclosure_checked <- TRUE
write.csv(routing_checks, file.path(review_dir, "10_routing_checks.csv"),
          row.names = FALSE)

if (any(routing_checks$imputations_passing < routing_checks$imputations)) {
  stop("A routing or consistency check failed; see 10_routing_checks.csv.")
}


###############################################################################
# Block 4. Trace plots
#
# For every variable imputed by a model: the mean and SD of its imputed
# values at each iteration. Each imputation is a thin grey line; the black
# line is their average. For categorical variables, mice records the mean
# and SD of the category numbers (1, 2, ...), which is enough to see a
# trend. All variables go to a private file; the review copy includes only
# variables with at least 10 imputed values.
###############################################################################

for (set in names(imputation_sets)) {
  imputation <- imputation_sets[[set]]$imputation
  traced     <- imputation_sets[[set]]$traced
  n_imputed  <- imputation_sets[[set]]$n_imputed

  for (copy in c("private", "review")) {
    if (copy == "private") {
      pdf(file.path(output_dir, paste0("private_10_trace_plots_", set, "_all.pdf")),
          width = 8.5, height = 11)
      plot_variables <- traced
    } else {
      pdf(file.path(review_dir, paste0("10_trace_plots_", set, ".pdf")),
          width = 8.5, height = 11)
      plot_variables <- traced[n_imputed[traced] >= disclosure_threshold]
    }
    par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
    for (variable in plot_variables) {
      chain_means <- imputation$chainMean[variable, , ]        # iterations x imputations
      chain_sds   <- sqrt(imputation$chainVar[variable, , ])
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
}


###############################################################################
# Block 5. Convergence table
#
# For the mean and the SD of each variable's imputed values, over the
# second half of the iterations, using every imputation:
#   split R-hat: each imputation's second half is cut into two, and the
#     variation between all the half-sequences is compared with the
#     variation within them;
#   trend: the average change per iteration, from a regression of the
#     statistic on iteration with a separate level for each imputation, and
#     its p-value (successive iterations are correlated, so the p-value
#     overstates the evidence);
#   drift_vs_spread: the total change over the second half divided by the
#     typical spread between imputations at an iteration.
# Flagged if split R-hat is above 1.1, or if the trend has p < 0.01 and
# drift_vs_spread is above 1. Only variables with at least 10 imputed
# values are included.
###############################################################################

convergence <- data.frame()

for (set in names(imputation_sets)) {
  imputation <- imputation_sets[[set]]$imputation
  traced     <- imputation_sets[[set]]$traced
  n_imputed  <- imputation_sets[[set]]$n_imputed
  iterations <- imputation$iteration
  m          <- imputation$m

  second_half <- (floor(iterations / 2) + 1):iterations        # e.g. 16-30
  half_length <- floor(length(second_half) / 2)
  early_part  <- head(second_half, half_length)
  late_part   <- tail(second_half, half_length)

  # One row per variable (at least 10 imputed values) x statistic.
  trace_table <- expand.grid(statistic = c("mean", "sd"),
                             variable = traced[n_imputed[traced] >= disclosure_threshold],
                             stringsAsFactors = FALSE)
  for (j in seq_len(nrow(trace_table))) {
    variable  <- trace_table$variable[j]
    statistic <- trace_table$statistic[j]
    values <- imputation$chainMean[variable, , ]                 # iterations x imputations
    if (statistic == "sd") values <- sqrt(imputation$chainVar[variable, , ])

    # Split R-hat: 2 x m half-sequences of equal length (needs at least 2
    # iterations in each half).
    split_rhat <- NA
    if (half_length >= 2) {
      halves  <- cbind(values[early_part, , drop = FALSE],
                       values[late_part, , drop = FALSE])
      n       <- nrow(halves)
      within  <- mean(apply(halves, 2, var))                     # variance within each half-sequence
      between <- n * var(colMeans(halves))                       # variance between half-sequence means
      if (within > 0) split_rhat <- sqrt(((n - 1) / n * within + between / n) / within)
    }

    # Trend over the second half of the iterations.
    trend_data <- data.frame(
      value      = as.vector(values[second_half, , drop = FALSE]),
      iteration  = rep(second_half, times = m),
      imputation = factor(rep(seq_len(m), each = length(second_half)))
    )
    trend <- summary(lm(value ~ iteration + imputation, data = trend_data))$coefficients
    trend_slope <- trend["iteration", "Estimate"]
    trend_p     <- trend["iteration", "Pr(>|t|)"]

    # Drift over the second half, relative to the spread between imputations
    # (the SD across imputations, averaged over iterations).
    spread <- mean(apply(values[second_half, , drop = FALSE], 1, sd))
    drift_vs_spread <- NA
    if (spread > 0) drift_vs_spread <- abs(trend_slope) * (length(second_half) - 1) / spread

    convergence <- rbind(convergence, data.frame(
      imputation_set      = imputation_sets[[set]]$label,
      variable            = variable,
      statistic           = statistic,
      imputations         = m,
      split_rhat          = round(split_rhat, 3),
      trend_per_iteration = signif(trend_slope, 3),
      trend_p             = signif(trend_p, 3),
      drift_vs_spread     = round(drift_vs_spread, 2),
      flag                = isTRUE(split_rhat > rhat_flag) ||
                            (isTRUE(trend_p < trend_p_flag) &&
                             isTRUE(drift_vs_spread > drift_flag))
    ))
  }
}

# Disclosure check: summaries of the means and SDs of imputed values, only
# for variables with at least 10 imputed values; no counts of people.
convergence$disclosure_checked <- rep(TRUE, nrow(convergence))
write.csv(convergence, file.path(review_dir, "10_convergence.csv"), row.names = FALSE)


###############################################################################
# Block 6. Observed versus imputed distributions
#
# For every variable imputed by a model, in each set: the distribution of
# its observed values and of its imputed values (all imputations together):
#   categorical (including 0/1): percentage in each category;
#   numeric: mean, SD and quartiles.
# For a routed later part, the imputed values include people whose imputed
# earlier answer doesn't lead to it; the calculated pathway shows the values
# that are used.
###############################################################################

comparison <- data.frame()

for (set in names(imputation_sets)) {
  spec       <- imputation_sets[[set]]$spec
  imputation <- imputation_sets[[set]]$imputation
  n_imputed  <- imputation_sets[[set]]$n_imputed
  m          <- imputation$m

  # The rows to write, in order: one per category of a categorical variable
  # (0/1 variables have categories "0" and "1"), and one "(numeric)" row per
  # numeric variable.
  row_table <- data.frame()
  for (variable in imputation_sets[[set]]$traced) {
    is_categorical <- is.factor(spec$data[[variable]]) || all(spec$data[[variable]] %in% c(0, 1, NA))
    categories <- "(numeric)"
    if (is_categorical) categories <- c("0", "1")
    if (is.factor(spec$data[[variable]])) categories <- levels(spec$data[[variable]])
    row_table <- rbind(row_table, data.frame(variable = variable, category = categories))
  }

  for (j in seq_len(nrow(row_table))) {
    variable <- row_table$variable[j]
    category <- row_table$category[j]
    observed_values <- spec$data[[variable]][!is.na(spec$data[[variable]])]
    # mice keeps each variable's imputed values as a table: one row per
    # imputed value, one column per imputation.
    imputed_values  <- as.character(unlist(imputation$imp[[variable]]))

    if (category != "(numeric)") {
      observed_n     <- sum(as.character(observed_values) == category)
      imputed_mean_n <- sum(imputed_values == category) / m
      comparison <- rbind(comparison, data.frame(
        imputation_set = imputation_sets[[set]]$label,
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
    } else {
      imputed_numbers <- as.numeric(imputed_values)
      observed_q <- quantile(observed_values, c(0.25, 0.5, 0.75), names = FALSE)
      imputed_q  <- quantile(imputed_numbers, c(0.25, 0.5, 0.75), names = FALSE)
      comparison <- rbind(comparison, data.frame(
        imputation_set = imputation_sets[[set]]$label,
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
}

write.csv(comparison, file.path(output_dir, "private_10_observed_vs_imputed.csv"),
          row.names = FALSE)

# --- Disclosure check: observed versus imputed --------------------------------
# Categories: a row's percentages are hidden if fewer than 10 people are
# observed in the category, or fewer than 10 values (average per
# imputation) are imputed into it. Secondary: within each variable and
# imputation set, if exactly one row is hidden, the row with the next
# smallest observed count is hidden too (percentages add to 100). Counts
# are removed from the shared table, so it can't be combined with other
# tables' counts.
# Numeric rows: observed statistics only if at least 10 are observed;
# imputed statistics only if at least 10 values are imputed.
# The number imputed per variable is hidden if it is 1-9.
checked <- comparison
is_numeric_row <- checked$category == "(numeric)"

primary <- !is_numeric_row &
  ((checked$observed_n > 0 & checked$observed_n < disclosure_threshold) |
   (checked$imputed_mean_n > 0 & checked$imputed_mean_n < disclosure_threshold))
secondary <- rep(FALSE, nrow(checked))
group <- paste(checked$imputation_set, checked$variable)
for (g in unique(group[!is_numeric_row])) {
  rows <- which(group == g)
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
write.csv(checked, file.path(review_dir, "10_observed_vs_imputed.csv"),
          row.names = FALSE, na = "")

hidden_from_plots <- 0
for (set in names(imputation_sets)) {
  traced <- imputation_sets[[set]]$traced
  hidden_from_plots <- hidden_from_plots +
    sum(imputation_sets[[set]]$n_imputed[traced] < disclosure_threshold)
}
disclosure_log <- data.frame(
  table = c("10_observed_vs_imputed", "10_trace_plots and 10_convergence"),
  primary_suppressed = c(sum(primary) + sum(numeric_observed_small) +
                           sum(numeric_imputed_small) + sum(small_n_imputed),
                         hidden_from_plots),
  secondary_suppressed = c(sum(secondary), 0L),
  threshold = disclosure_threshold
)
write.csv(disclosure_log, file.path(review_dir, "10_disclosure_log.csv"), row.names = FALSE)
# -----------------------------------------------------------------------------


###############################################################################
# Block 7. Summary
#
# Counts of imputations, variables and flags only. No data values.
###############################################################################

cat("\nFile 10 finished.\n")
for (set in names(imputation_sets)) {
  label <- imputation_sets[[set]]$label
  cat(" ", label, "\n")
  cat("    imputations:                        ", imputation_sets[[set]]$imputation$m, "\n")
  cat("    iterations:                         ", imputation_sets[[set]]$imputation$iteration, "\n")
  cat("    routing and consistency checks:      all passed\n")
  rows <- convergence$imputation_set == label
  cat("    convergence flags:                  ", sum(convergence$flag[rows]), "of",
      sum(rows), "variable x statistic rows\n")
  cat("      of which with a trend (p < 0.01, drift > spread):",
      sum(rows & convergence$trend_p < trend_p_flag &
            convergence$drift_vs_spread > drift_flag, na.rm = TRUE), "\n")
}
cat("  Review: 10_routing_checks.csv, 10_convergence.csv,\n")
cat("    10_trace_plots_main.pdf, 10_trace_plots_supplementary.pdf,\n")
cat("    10_observed_vs_imputed.csv\n")
