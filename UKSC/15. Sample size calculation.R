###############################################################################
# 15. Sample size calculation
#
# Purpose
#   For each prediction landmark, calculate the minimum sample size needed to
#   develop a prediction model for a binary outcome (Riley et al. 2020), and
#   compare it with the data available. This is TRIPOD+AI item 10 (how the
#   study size was arrived at). No model is fitted and no performance is
#   computed.
#
# The criteria (Riley et al. 2020)
#   The required sample size is the largest of three:
#     1. small overfitting: an expected shrinkage of the coefficients of at
#        least 0.9 (overfitting would shrink them by less than 10%);
#     2. a small difference (0.05 or less) between the apparent and adjusted
#        Nagelkerke R-squared;
#     3. a precise estimate of the overall outcome proportion (margin of
#        error 0.05 or less).
#   The calculations use the pmsampsize package.
#
# Inputs
#   - Outcome proportion: the unweighted SRET proportion in the prediction
#     cohort. The models are fitted unweighted (File 17), so the unweighted
#     proportion is the relevant one.
#   - Candidate parameters at each landmark (File 14): one per numeric
#     feature, and one per category beyond the first for categorical
#     features, before any pooling of rare categories.
#   - Anticipated model performance, fixed before this script was run:
#       Primary:     Nagelkerke R-squared = 0.15 at every landmark. Riley et
#                    al. (2020) suggest this value when there is no existing
#                    model to draw on, and no published model predicts SRET
#                    employment at 32 from these measures.
#       Sensitivity: anticipated C-statistics of 0.70, 0.75 and 0.80,
#                    converted to Cox-Snell R-squared as described by Riley,
#                    Van Calster and Collins (2021).
#
# Data available
#   The models are developed on outer training sets (File 17), each about
#   four-fifths of the cohort. The smallest outer training set, and the
#   fewest events in any outer training set, are compared with the
#   requirement. The whole cohort is also shown, because the final models
#   for reporting coefficients are fitted to it (File 19).
#   For each landmark the table also gives:
#     - the most parameters the available data support at a shrinkage of
#       0.9 (criterion 1 solved for the number of parameters);
#     - the expected shrinkage with the landmark's actual number of
#       candidate parameters (criterion 1 solved for the shrinkage).
#
# If a landmark falls short
#   The landmarks keep their features: each landmark's feature set is the
#   information it represents, and reducing it to meet the calculation would
#   change the question the landmark answers. A shortfall is reported with
#   its expected shrinkage and handled as follows:
#     - the primary model is the elastic net, whose penalty is tuned to
#       shrink the coefficients (File 17);
#     - overfitting that remains is measured directly by the calibration
#       slope in the outer test folds (File 18);
#     - the unpenalized logistic benchmark is expected to overfit more at the
#       later landmarks, which have more parameters.
#   Overfitting lowers performance in new data, so a shortfall at the later
#   landmarks would make their added value look smaller rather than larger.
#   Penalization doesn't fully make up for a small sample (Riley et al.
#   2021), which is stated as a limitation.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   The table contains the outcome proportion of the whole prediction cohort,
#   and the people and events in the smallest outer training set and the
#   whole cohort. These use the same cohort and folds as Files 4 and 14. A
#   row's counts are shown only if the group has at least 10 events and 10
#   non-events.
#
# Reads
#   data_dir/14_prediction_design.rds
#
# Writes
#   output_dir/private_15_sample_size.csv   exact table
#   review_dir/15_sample_size.csv           checked table
#   review_dir/15_disclosure_log.csv        what was suppressed
#
# References
#   Riley RD, et al. (2020). Calculating the sample size required for
#   developing a clinical prediction model. BMJ 368: m441.
#   Riley RD, Van Calster B, Collins GS (2021). A note on estimating the
#   Cox-Snell R2 from a reported C statistic (AUROC) to inform sample size
#   calculations for developing a prediction model with a binary outcome.
#   Statistics in Medicine 40(4): 859-864.
#   Riley RD, et al. (2021). Penalization and shrinkage methods produced
#   unreliable clinical prediction models especially when sample size was
#   small. Journal of Clinical Epidemiology 132: 88-96.
#   Collins GS, et al. (2024). TRIPOD+AI statement. BMJ 385: e078378.
###############################################################################

suppressPackageStartupMessages(library(pmsampsize))


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
target_shrinkage     <- 0.9    # criterion 1 (Riley et al. 2020)

# Anticipated performance, fixed before running (see the header).
assumptions <- data.frame(
  assumption = c("Nagelkerke R-squared 0.15 (primary)",
                 "C-statistic 0.70 (sensitivity)",
                 "C-statistic 0.75 (sensitivity)",
                 "C-statistic 0.80 (sensitivity)"),
  measure    = c("nagelkerke", "cstatistic", "cstatistic", "cstatistic"),
  value      = c(0.15, 0.70, 0.75, 0.80),
  primary    = c(TRUE, FALSE, FALSE, FALSE)
)


###############################################################################
# Block 2. Read the prediction design from File 14
###############################################################################

design           <- readRDS(file.path(data_dir, "14_prediction_design.rds"))
cohort           <- design$cohort
landmark_summary <- design$landmark_summary
settings         <- design$settings

y <- cohort[[settings$outcome]]
outcome_proportion <- mean(y)   # unweighted, whole prediction cohort


###############################################################################
# Block 3. Data available
#
# The smallest outer training set (fewest people), and the fewest events in
# any outer training set, across all repeats and folds. Also the whole
# cohort.
###############################################################################

training_people <- c()
training_events <- c()
training_non_events <- c()
for (r in 1:settings$n_repeats) {
  for (f in 1:settings$n_outer_folds) {
    in_training <- cohort[[paste0("outer_fold_r", r)]] != f
    training_people     <- c(training_people, sum(in_training))
    training_events     <- c(training_events, sum(y[in_training] == 1))
    training_non_events <- c(training_non_events, sum(y[in_training] == 0))
  }
}

available <- data.frame(
  data       = c("Smallest outer training set", "Whole prediction cohort"),
  people     = c(min(training_people), length(y)),
  events     = c(min(training_events), sum(y == 1)),
  non_events = c(min(training_non_events), sum(y == 0))
)


###############################################################################
# Block 4. Required sample size for each landmark and assumption
#
# For each landmark x anticipated performance:
#   a. pmsampsize gives the required people and events (the largest of
#      criteria 1-3), the events per parameter, and the anticipated
#      Cox-Snell R-squared;
#   b. for each data source (smallest training set, whole cohort):
#        - whether the required number of people is available;
#        - the most parameters supported at a shrinkage of 0.9:
#            parameters = n x (S - 1) x ln(1 - R2 / S), with S = 0.9;
#        - the expected shrinkage with the landmark's actual parameters:
#            the S that solves n = parameters / ((S - 1) x ln(1 - R2 / S)).
#      These are Riley et al.'s (2020) criterion 1 solved for the number of
#      parameters and for the shrinkage. The shrinkage is found numerically
#      (uniroot), searching between the anticipated R-squared and 1.
###############################################################################

results <- data.frame()

for (i in seq_len(nrow(landmark_summary))) {
  landmark   <- landmark_summary$landmark[i]
  parameters <- landmark_summary$candidate_parameters[i]

  for (k in seq_len(nrow(assumptions))) {

    # a. Required sample size.
    if (assumptions$measure[k] == "nagelkerke") {
      calculation <- pmsampsize(type = "b", nagrsquared = assumptions$value[k],
                                parameters = parameters,
                                prevalence = outcome_proportion,
                                shrinkage = target_shrinkage)
    } else {
      calculation <- pmsampsize(type = "b", cstatistic = assumptions$value[k],
                                parameters = parameters,
                                prevalence = outcome_proportion,
                                shrinkage = target_shrinkage)
    }
    r2 <- calculation$rsquared   # anticipated Cox-Snell R-squared
    limiting <- rownames(calculation$results_table)[
      which.max(calculation$results_table[c("Criteria 1", "Criteria 2", "Criteria 3"),
                                          "Samp_size"])]

    # b. Compared with each data source.
    for (d in seq_len(nrow(available))) {
      n <- available$people[d]
      max_parameters <- floor(n * (target_shrinkage - 1) *
                                log(1 - r2 / target_shrinkage))
      # The expected shrinkage S solves n = parameters / ((S - 1) log(1 -
      # R2 / S)) (Riley et al. 2020). uniroot() finds S numerically; it needs
      # the equation written as a small function of S that is 0 at the answer.
      expected_shrinkage <- uniroot(
        function(S) n - parameters / ((S - 1) * log(1 - r2 / S)),
        interval = c(r2 * 1.001, 0.99999)
      )$root

      results <- rbind(results, data.frame(
        landmark                 = landmark,
        assumption               = assumptions$assumption[k],
        primary                  = assumptions$primary[k],
        candidate_parameters     = parameters,
        anticipated_cox_snell_r2 = round(r2, 4),
        required_people          = calculation$sample_size,
        required_events          = ceiling(calculation$events),
        required_events_per_parameter = round(calculation$EPP, 2),
        limiting_criterion       = limiting,
        data                     = available$data[d],
        available_people         = n,
        available_events         = available$events[d],
        available_non_events     = available$non_events[d],
        requirement_met          = n >= calculation$sample_size,
        max_parameters_at_0.9    = max_parameters,
        expected_shrinkage       = round(expected_shrinkage, 3)
      ))
    }
  }
}

results$outcome_proportion <- round(outcome_proportion, 4)

write.csv(results, file.path(output_dir, "private_15_sample_size.csv"),
          row.names = FALSE)


###############################################################################
# Block 5. Disclosure check and review table
###############################################################################

# --- Disclosure check: sample size table -------------------------------------
# The available people, events and non-events are shown only if the group
# has at least 10 events and 10 non-events; otherwise they, the outcome
# proportion and the comparisons that depend on them are hidden. The
# required sample sizes are calculations, not counts of people, and are
# always shown.
checked <- results
small <- checked$available_events < disclosure_threshold |
         checked$available_non_events < disclosure_threshold
checked[small, c("available_people", "available_events", "available_non_events",
                 "requirement_met", "max_parameters_at_0.9",
                 "expected_shrinkage", "outcome_proportion")] <- NA
checked <- checked[, setdiff(names(checked), "available_non_events")]
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "15_sample_size.csv"),
          row.names = FALSE, na = "")

disclosure_log <- data.frame(table = "15_sample_size",
                             primary_suppressed = sum(small),
                             secondary_suppressed = 0L,
                             threshold = disclosure_threshold)
write.csv(disclosure_log, file.path(review_dir, "15_disclosure_log.csv"),
          row.names = FALSE)
# -----------------------------------------------------------------------------


###############################################################################
# Block 6. Summary
#
# The primary assumption, compared with the smallest outer training set.
###############################################################################

summary_rows <- results[results$primary & results$data == "Smallest outer training set", ]

cat("\nFile 15 finished.\n")
cat("  Primary assumption: Nagelkerke R-squared 0.15; target shrinkage 0.9\n")
cat("  Compared with the smallest outer training set:\n")
cat("    landmark: parameters / required people / met / most parameters supported / expected shrinkage\n")
for (i in seq_len(nrow(summary_rows))) {
  cat("    ", summary_rows$landmark[i], ": ",
      summary_rows$candidate_parameters[i], " / ",
      summary_rows$required_people[i], " / ",
      ifelse(summary_rows$requirement_met[i], "yes", "no"), " / ",
      summary_rows$max_parameters_at_0.9[i], " / ",
      summary_rows$expected_shrinkage[i], "\n", sep = "")
}
cat("  Full table (all assumptions, both data sources): 15_sample_size.csv\n")
