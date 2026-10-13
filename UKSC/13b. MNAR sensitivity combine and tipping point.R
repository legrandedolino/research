###############################################################################
# 13b. MNAR sensitivity combine and tipping point
#
# Purpose
#   Combine the scenarios of File 13a with Rubin's rules and report how the
#   key estimates change as the data move away from missing at random (MAR),
#   and the tipping points: the smallest departure, in each direction, at
#   which a conclusion changes (CLS user guide).
#
# Conclusions and tipping points (fixed before the results)
#   Associations (the three science risk differences, and the odds ratio for
#   liking mathematics): the conclusion is whether the 95% confidence
#   interval excludes no association (0 for a risk difference, 1 for an odds
#   ratio). The tipping point in each direction is the delta closest to 0 at
#   which that conclusion differs from the MAR scenario (delta = 0).
#   Prevalences: the tipping point is the delta closest to 0 at which the
#   estimate falls outside the MAR scenario's 95% confidence interval.
#   Each family is read for the estimand it targets: uniform for every
#   estimand; liking, competence and favourite for their science contrast;
#   maths for the mathematics odds ratio. Every family x estimand is also
#   shown in full.
#   If a conclusion doesn't change within delta = -2 to +2 (odds ratios 0.14
#   to 7.4), it is reported as robust within that range. Larger departures
#   are implausible for otherwise similar people.
#   How plausible a tipping point is must be judged in context: delta of
#   +/-0.5 (odds ratios of about 0.6 and 1.6) is a moderate departure from
#   MAR; +/-1 or more (about 0.37 and 2.7) is a large one.
#
# Rubin's rules and checks
#   As File 12: Barnard-Rubin degrees of freedom with the survey design's as
#   the complete-data value. Checks: every scenario file is present, with
#   every imputation and the same R and mice versions. The MAR scenario
#   (uniform, delta = 0) is compared with File 12's primary estimates; they
#   differ only by Monte Carlo noise, because the outcome is drawn again.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   Estimates are shown only if every completed dataset in every scenario has
#   at least 10 events and 10 non-events for that estimand. No counts are
#   shown. The plots show only pooled estimates and their intervals.
#
# Reads
#   tasks_dir/13a/scenario_NNN.rds
#   data_dir/12_pooled_estimates.rds   (for the MAR comparison)
#
# Writes
#   data_dir/13_mnar_results.rds       exact pooled results
#   output_dir/private_13_*.csv        exact tables (never shared)
#   review_dir/13_mnar_estimates.csv, 13_tipping_points.csv,
#     13_mnar_plots.pdf, 13_task_log.csv, 13_disclosure_log.csv
###############################################################################


###############################################################################
# Block 1. Paths and settings
###############################################################################

paths      <- readRDS("project_paths.rds")
output_dir <- paths$output_dir   # exact (private) tables are written here
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs
tasks_dir  <- paths$tasks_dir    # one file per task in split runs

disclosure_threshold <- 10
delta_grid  <- seq(-2, 2, by = 0.5)     # as in File 13a
families    <- c("uniform", "liking", "competence", "favourite", "maths")
task_folder <- file.path(tasks_dir, "13a")

# The exposure each differential family targets (its estimand is the one
# File 13a recorded for that exposure).
family_exposure <- c(liking = "science_liking14", competence = "science_competence14",
                     favourite = "science_favourite14", maths = "math_liking14")


###############################################################################
# Block 2. Read and check the scenario files
#
# The scenario table is rebuilt exactly as in File 13a.
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

missing_tasks <- scenarios$task[!file.exists(scenarios$file)]
if (length(missing_tasks) > 0) {
  stop(length(missing_tasks), " scenario files are missing (tasks ",
       paste(missing_tasks, collapse = ", "), "). Run them with File 13a.")
}

results  <- data.frame()
task_log <- data.frame()
for (t in scenarios$task) {
  task_result <- readRDS(scenarios$file[t])
  if (task_result$task != t || task_result$family != scenarios$family[t] ||
      task_result$delta != scenarios$delta[t]) {
    stop("Scenario file ", t, " doesn't match the scenario table.")
  }
  task_result$results$family <- task_result$family
  task_result$results$delta  <- task_result$delta
  results <- rbind(results, task_result$results)
  task_log <- rbind(task_log, data.frame(
    task = t, family = task_result$family, delta = task_result$delta,
    imputations = length(unique(task_result$results$imputation)),
    minutes = round(task_result$minutes, 1), computer = task_result$computer,
    r_version = task_result$r_version, mice_version = task_result$mice_version
  ))
}

if (length(unique(task_log$imputations)) != 1) {
  stop("Scenarios have different numbers of imputations (see 13_task_log.csv).")
}
if (length(unique(task_log$r_version)) > 1 || length(unique(task_log$mice_version)) > 1) {
  stop("Scenarios were run with different versions of R or mice.")
}
m <- task_log$imputations[1]


###############################################################################
# Block 3. Rubin's rules for every scenario x estimand
###############################################################################

# One key per scenario x estimand; each key's rows are its values in the
# completed datasets.
key <- paste(results$family, results$delta, results$estimand, sep = " | ")

pooled <- data.frame()
for (g in unique(key)) {
  rows <- results[key == g, ]

  q_bar   <- mean(rows$estimate)
  within  <- mean(rows$variance)
  between <- var(rows$estimate)
  total   <- within + (1 + 1 / m) * between
  lambda  <- (1 + 1 / m) * between / total
  df_complete <- rows$df_complete[1]
  df_observed <- (df_complete + 1) / (df_complete + 3) * df_complete * (1 - lambda)
  if (lambda > 0) {
    df_old <- (m - 1) / lambda^2
    df     <- df_old * df_observed / (df_old + df_observed)
  } else {
    df <- df_observed
  }
  se <- sqrt(total)
  critical <- qt(0.975, df = df)

  pooled <- rbind(pooled, data.frame(
    family = rows$family[1], delta = rows$delta[1], estimand = rows$estimand[1],
    scale = rows$scale[1], exposure = rows$exposure[1], estimate = q_bar, se = se,
    low = q_bar - critical * se, high = q_bar + critical * se, df = df,
    smallest_events = min(rows$events), smallest_non_events = min(rows$non_events)
  ))
}

# The MAR scenario belongs to every family.
mar <- pooled[pooled$family == "uniform" & pooled$delta == 0, ]
for (f in setdiff(families, "uniform")) {
  mar_f <- mar
  mar_f$family <- f
  pooled <- rbind(pooled, mar_f)
}
pooled <- pooled[order(match(pooled$family, families), pooled$estimand, pooled$delta), ]
rownames(pooled) <- NULL

# Odds ratios are shown on the odds ratio scale; prevalences as percentages;
# risk differences as they are.
is_or  <- pooled$scale == "log odds ratio"
is_pct <- pooled$scale == "proportion"
pooled$shown_estimate <- pooled$estimate
pooled$shown_low      <- pooled$low
pooled$shown_high     <- pooled$high
pooled$shown_estimate[is_or]  <- exp(pooled$estimate[is_or])
pooled$shown_low[is_or]       <- exp(pooled$low[is_or])
pooled$shown_high[is_or]      <- exp(pooled$high[is_or])
pooled$shown_estimate[is_pct] <- 100 * pooled$estimate[is_pct]
pooled$shown_low[is_pct]      <- 100 * pooled$low[is_pct]
pooled$shown_high[is_pct]     <- 100 * pooled$high[is_pct]


###############################################################################
# Block 4. Conclusions and tipping points
#
# For associations, "evidence" = the 95% CI excludes no association. For
# prevalences, the reference is the MAR scenario's 95% CI. In each family and
# direction, the tipping point is the delta closest to 0 whose conclusion
# differs from MAR's; "none within +/-2" if there is none.
###############################################################################

pooled$is_association <- pooled$scale != "proportion"
pooled$evidence <- NA                                   # not applicable to prevalences
excludes_zero <- pooled$low > 0 | pooled$high < 0
pooled$evidence[pooled$is_association] <- excludes_zero[pooled$is_association]

tipping <- data.frame()
for (f in families) {
  # The uniform family shifts every estimand; the others shift their own
  # exposure's estimands.
  estimands <- unique(pooled$estimand)
  if (f != "uniform") estimands <- unique(pooled$estimand[pooled$exposure %in% family_exposure[[f]]])
  if (length(estimands) == 0) stop("No estimand found for family ", f, ".")
  for (e in estimands) {
    rows <- pooled[pooled$family == f & pooled$estimand == e, ]
    mar_row <- rows[rows$delta == 0, ]
    # Which deltas change the conclusion: for an association, whether its
    # 95% CI excludes no association; for a prevalence, whether it leaves
    # the MAR 95% CI.
    if (mar_row$is_association) {
      changed <- rows$evidence != mar_row$evidence
      mar_conclusion <- "95% CI includes no association"
      if (mar_row$evidence) mar_conclusion <- "95% CI excludes no association"
    } else {
      changed <- rows$estimate < mar_row$low | rows$estimate > mar_row$high
      mar_conclusion <- "reference: MAR 95% CI"
    }

    # The tipping points: the changed deltas closest to 0 on each side.
    negative <- rows$delta[rows$delta < 0 & changed]
    positive <- rows$delta[rows$delta > 0 & changed]
    tipping_negative <- "none within -2"
    tipping_positive <- "none within +2"
    odds_ratio_negative <- ""
    odds_ratio_positive <- ""
    if (length(negative) > 0) tipping_negative    <- as.character(max(negative))
    if (length(negative) > 0) odds_ratio_negative <- as.character(round(exp(max(negative)), 2))
    if (length(positive) > 0) tipping_positive    <- as.character(min(positive))
    if (length(positive) > 0) odds_ratio_positive <- as.character(round(exp(min(positive)), 2))

    tipping <- rbind(tipping, data.frame(
      family = f, estimand = e, conclusion_under_mar = mar_conclusion,
      tipping_delta_negative = tipping_negative,
      tipping_delta_positive = tipping_positive,
      tipping_odds_ratio_negative = odds_ratio_negative,
      tipping_odds_ratio_positive = odds_ratio_positive,
      smallest_events = min(rows$smallest_events),
      smallest_non_events = min(rows$smallest_non_events)
    ))
  }
}


###############################################################################
# Block 5. Check: the MAR scenario against File 12
#
# The same estimands from File 12's primary analysis. The differences come
# only from drawing the imputed outcome again; they are summarized, not
# shown.
###############################################################################

file12 <- readRDS(file.path(data_dir, "12_pooled_estimates.rds"))
primary <- file12$analyses$analysis[1]
f12_rd  <- file12$risk_differences[file12$risk_differences$analysis == primary &
                                     file12$risk_differences$model == "M0 baseline", ]
mar_rd  <- mar[mar$scale == "risk difference", ]
rd_difference <- max(abs(mar_rd$estimate -
                         f12_rd$estimate[match(sub("^Risk difference: ", "", mar_rd$estimand), f12_rd$term)]))
f12_prev <- file12$prevalence[file12$prevalence$analysis == primary &
                                file12$prevalence$model == "SRET overall", "estimate"]
prevalence_difference <- abs(mar$estimate[mar$estimand == "Prevalence: SRET overall"] - f12_prev)

write.csv(pooled, file.path(output_dir, "private_13_mnar_estimates.csv"), row.names = FALSE)
write.csv(tipping, file.path(output_dir, "private_13_tipping_points.csv"), row.names = FALSE)
saveRDS(list(scenarios = scenarios, results = results, pooled = pooled, tipping = tipping),
        file.path(data_dir, "13_mnar_results.rds"))


###############################################################################
# Block 6. Disclosure checks and review tables
###############################################################################

disclosure_log <- data.frame(table = character(), primary_suppressed = integer())

# --- Disclosure check: estimates ----------------------------------------------
# An estimand is shown only if every completed dataset in every scenario has
# at least 10 events and 10 non-events for it; counts removed.
estimand_small <- unique(pooled$estimand[pooled$smallest_events < disclosure_threshold |
                                           pooled$smallest_non_events < disclosure_threshold])
checked <- pooled
hide <- checked$estimand %in% estimand_small
checked[hide, c("shown_estimate", "shown_low", "shown_high", "evidence")] <- NA
checked <- data.frame(
  family = checked$family, delta = checked$delta,
  odds_ratio_of_delta = round(exp(checked$delta), 2),
  estimand = checked$estimand,
  estimate = round(checked$shown_estimate, 4), low = round(checked$shown_low, 4),
  high = round(checked$shown_high, 4),
  evidence_of_association = checked$evidence
)
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "13_mnar_estimates.csv"), row.names = FALSE, na = "")
disclosure_log <- rbind(disclosure_log, data.frame(table = "13_mnar_estimates",
                                                   primary_suppressed = sum(hide)))
# -----------------------------------------------------------------------------

# --- Disclosure check: tipping points -----------------------------------------
checked <- tipping
hide <- checked$estimand %in% estimand_small
checked[hide, c("tipping_delta_negative", "tipping_delta_positive",
                "tipping_odds_ratio_negative", "tipping_odds_ratio_positive")] <- NA
checked <- checked[, setdiff(names(checked), c("smallest_events", "smallest_non_events"))]
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "13_tipping_points.csv"), row.names = FALSE, na = "")
disclosure_log <- rbind(disclosure_log, data.frame(table = "13_tipping_points",
                                                   primary_suppressed = sum(hide)))
# -----------------------------------------------------------------------------

# Task log: scenarios, run times, versions; no data values.
task_log$disclosure_checked <- TRUE
write.csv(task_log, file.path(review_dir, "13_task_log.csv"), row.names = FALSE)


###############################################################################
# Block 7. Plots
#
# One page per family: each estimand's pooled estimate (points) and 95% CI
# (bars) against delta, with the no-association line (dashed) and the MAR
# estimate (dotted). Only estimands that pass the disclosure check are drawn.
###############################################################################

pdf(file.path(review_dir, "13_mnar_plots.pdf"), width = 11, height = 8.5)
for (f in families) {
  rows_f <- pooled[pooled$family == f & !(pooled$estimand %in% estimand_small), ]
  estimands <- unique(rows_f$estimand)
  par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  for (e in estimands) {
    rows <- rows_f[rows_f$estimand == e, ]
    plot(rows$delta, rows$shown_estimate, pch = 19,
         ylim = range(c(rows$shown_low, rows$shown_high)),
         xlab = "delta (log odds ratio, missing vs observed)",
         ylab = if (rows$scale[1] == "proportion") "Weighted %" else
           if (rows$scale[1] == "log odds ratio") "Odds ratio" else "Risk difference",
         main = sub("^(Risk difference|Odds ratio|Prevalence): ", "", e), cex.main = 0.8)
    arrows(rows$delta, rows$shown_low, rows$delta, rows$shown_high,
           angle = 90, code = 3, length = 0.03)
    if (rows$scale[1] == "log odds ratio") abline(h = 1, lty = 2)
    if (rows$scale[1] == "risk difference") abline(h = 0, lty = 2)
    abline(h = rows$shown_estimate[rows$delta == 0], lty = 3)
  }
  mtext(paste("Family:", f), outer = TRUE, cex = 1)
}
dev.off()

disclosure_log <- rbind(disclosure_log, data.frame(
  table = c("13_mnar_plots", "13_task_log"),
  primary_suppressed = c(length(estimand_small), 0L)))
disclosure_log$threshold <- disclosure_threshold
write.csv(disclosure_log, file.path(review_dir, "13_disclosure_log.csv"), row.names = FALSE)


###############################################################################
# Block 8. Summary
#
# Checks and counts only. No estimates.
###############################################################################

cat("\nFile 13b finished.\n")
cat("  Scenarios combined:                    ", nrow(scenarios), "(", m, "imputations each )\n")
cat("  MAR scenario vs File 12 (Monte Carlo noise only):\n")
cat("    largest difference in a risk difference:", signif(rd_difference, 2), "\n")
cat("    difference in overall prevalence:       ", signif(prevalence_difference, 2), "\n")
cat("  Total run time of 13a (hours):         ", round(sum(task_log$minutes) / 60, 1), "\n")
cat("  Review: 13_tipping_points.csv, 13_mnar_estimates.csv, 13_mnar_plots.pdf,\n")
cat("    13_task_log.csv\n")
