###############################################################################
# 4. Cohort and missingness description
#
# Purpose
#   Describe the cohort before any modelling:
#     Table 1  response at each sweep (sample flow)
#     Table 2  the age-32 outcome: all Wave-9 respondents by outcome group,
#              and SRET prevalence (unweighted and weighted) with effective
#              sample size in the prediction cohort
#     Table 3  missingness state of every analysis variable, in the
#              imputation population and in the prediction cohort
#     Table 4  prediction cohort and outcome events at each landmark, with
#              complete-case counts (for reference only; prediction doesn't
#              use complete cases)
#     Table 5  survey weight distributions
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   Every table is first saved exactly as computed, to a private file in
#   output_dir, which is never shared. A disclosure-checked copy is then
#   saved to output_dir/review. Each check is written out in full just
#   before the table is saved:
#     - primary suppression: any count from 1 to 9 is replaced by NA, along
#       with its percentage and any weighted estimate built on it;
#     - secondary suppression: if a group has exactly one suppressed cell,
#       it could be recovered from the group total, so the next smallest
#       count in that group is also suppressed;
#     - no minima or maxima are reported, and weighted prevalences are shown
#       only when both events and non-events are at least 10;
#     - tables share one definition of each group (e.g. the prediction
#       cohort), so a small group can't be recovered by subtracting one
#       table from another;
#     - a disclosure_checked column marks the checked tables, and every
#       suppression is recorded in 04_disclosure_log.csv.
#   Zero counts are not suppressed: a zero here describes an empty response
#   category, not a small group of people.
#
# Reads
#   data_dir/03_missingness_and_landmarks.rds
#
# Writes
#   output_dir/private_04_*.csv       exact tables (never shared)
#   review_dir/04_*.csv               disclosure-checked tables
#   review_dir/04_disclosure_log.csv  what was suppressed in each table
###############################################################################

library(haven)   # needed to read the deposited (SPSS-labelled) variables
suppressPackageStartupMessages(library(dplyr))


###############################################################################
# Block 1. Paths and disclosure threshold#
# Folder locations are set once, in File 0, and read here from
# project_paths.rds. Run this script with the working directory set to the
# folder that holds the scripts.
###############################################################################

paths      <- readRDS("project_paths.rds")
output_dir <- paths$output_dir   # exact (private) tables are written here
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs

disclosure_threshold <- 10   # smallest count that may be shown


###############################################################################
# Block 2. Read the cohort and landmark definitions from File 3
###############################################################################

file3 <- readRDS(file.path(data_dir, "03_missingness_and_landmarks.rds"))
cohort            <- file3$cohort
state_levels      <- file3$state_levels
landmark_features <- file3$landmark_features

# The prediction cohort (in_prediction_cohort, defined in File 3): a
# productive Wave-9 interview, an observed age-32 outcome, a positive Wave-9
# final weight, and a known PSU and stratum. Tables 2, 3 and 4 all use this
# one definition, so no small group can be recovered by subtracting one
# table from another.

disclosure_log <- data.frame(table = character(), primary_suppressed = integer(),
                             secondary_suppressed = integer())


###############################################################################
# Block 3. Table 1: response at each sweep
#
# For each sweep, the number of original-cohort members in each response
# outcome, and their percentage of the original cohort.
###############################################################################

sweep_age <- c(14, 15, 16, 17, 18, 19, 20, 25, 32)
response_outcomes <- data.frame(
  code    = c(1, 2, 3, 5, 4, -1),
  outcome = c("Productive", "Refusal", "Non-contact and other unproductive",
              "Untraced", "Ineligible", "Not issued")
)

sample_flow <- data.frame()
for (wave in 1:9) {
  outcome_code <- as.numeric(cohort[[paste0("W", wave, "OUTCOME")]])
  for (i in seq_len(nrow(response_outcomes))) {
    sample_flow <- rbind(sample_flow, data.frame(
      sweep   = wave,
      age     = sweep_age[wave],
      outcome = response_outcomes$outcome[i],
      n       = sum(outcome_code == response_outcomes$code[i], na.rm = TRUE),
      total   = length(outcome_code)
    ))
  }
  # Any code outside the six above is counted separately (expected to be 0).
  sample_flow <- rbind(sample_flow, data.frame(
    sweep = wave, age = sweep_age[wave], outcome = "Other or missing code",
    n = sum(is.na(outcome_code) | !(outcome_code %in% response_outcomes$code)),
    total = length(outcome_code)
  ))
}
sample_flow$percent <- round(100 * sample_flow$n / sample_flow$total, 1)

write.csv(sample_flow, file.path(output_dir, "private_04_sample_flow.csv"),
          row.names = FALSE)

# --- Disclosure check: Table 1 -----------------------------------------------
# Primary: counts 1-9 and their percentages.
checked <- sample_flow
primary <- checked$n > 0 & checked$n < disclosure_threshold
# Secondary: within each sweep, if exactly one count is suppressed, also
# suppress the smallest other non-zero count (the sweep total is known).
# The Wave-9 "Ineligible" count is never chosen as that second cell: the
# target population is the original cohort minus this group, so Table 3's
# total reveals it anyway.
w9_ineligible_row <- checked$sweep == 9 & checked$outcome == "Ineligible"
secondary <- rep(FALSE, nrow(checked))
for (wave in unique(checked$sweep)) {
  rows <- which(checked$sweep == wave)
  if (sum(primary[rows]) == 1) {
    candidates <- rows[!primary[rows] & checked$n[rows] > 0 & !w9_ineligible_row[rows]]
    if (length(candidates) > 0) {
      secondary[candidates[which.min(checked$n[candidates])]] <- TRUE
    }
  }
}
checked$n[primary | secondary]       <- NA
checked$percent[primary | secondary] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "04_sample_flow.csv"), row.names = FALSE)
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "04_sample_flow", primary_suppressed = sum(primary),
  secondary_suppressed = sum(secondary)))
# -----------------------------------------------------------------------------


###############################################################################
# Block 4. Table 2: the age-32 outcome
#
# 2a  every Wave-9 respondent, in exactly one of six groups:
#       no final weight or design information (PSU/stratum);
#       and, among those with them: SRET professional job / other job /
#       not in paid work / in paid work but occupation unresolved /
#       employment status unresolved
#     Because 2a covers all respondents, every small group is counted (and
#     suppressed) directly, rather than being recoverable by subtraction.
# 2b  SRET overall, in the prediction cohort: n, events, unweighted % and
#     W9FINWT-weighted %, plus Kish's effective sample size of the weights
# 2c  SRET among the employed with a valid occupation code (the
#     employed-only sensitivity outcome), with the same measures
###############################################################################

respondents <- cohort[cohort$responded_w9, ]
has_weight_and_design <- !is.na(respondents$w9_final_weight) &
  !is.na(respondents$psu) & !is.na(respondents$stratum)

# 2a. All Wave-9 respondents, in six groups.
outcome_group32 <- case_when(
  !has_weight_and_design                          ~ "No final weight or design information",
  respondents$sret_professional_32 %in% 1         ~ "SRET professional job",
  respondents$sret_among_employed32 %in% 0        ~ "Other job",
  respondents$employed32 %in% "Not in paid work"  ~ "Not in paid work",
  respondents$employed32 %in% "In paid work"      ~ "In paid work, occupation unresolved",
  TRUE                                            ~ "Employment status unresolved"
)
outcome_groups <- c("SRET professional job", "Other job", "Not in paid work",
                    "In paid work, occupation unresolved",
                    "Employment status unresolved",
                    "No final weight or design information")

respondent_groups32 <- data.frame(group = outcome_groups, n = NA_integer_,
                                  total = nrow(respondents))
for (i in seq_along(outcome_groups)) {
  respondent_groups32$n[i] <- sum(outcome_group32 == outcome_groups[i])
}
respondent_groups32$percent <- round(100 * respondent_groups32$n /
                                       respondent_groups32$total, 1)

# 2b and 2c. SRET prevalence in the prediction cohort, overall and among the
# employed.
prediction_cohort <- cohort[cohort$in_prediction_cohort, ]
outcome_summary <- data.frame()
for (outcome_name in c("sret_professional_32", "sret_among_employed32")) {
  observed <- prediction_cohort[!is.na(prediction_cohort[[outcome_name]]), ]
  y <- observed[[outcome_name]]
  w <- observed$w9_final_weight
  outcome_summary <- rbind(outcome_summary, data.frame(
    outcome            = outcome_name,
    n_observed         = nrow(observed),
    events             = sum(y == 1),
    non_events         = sum(y == 0),
    percent_unweighted = round(100 * mean(y), 1),
    percent_weighted   = round(100 * sum(w * y) / sum(w), 1),
    kish_effective_n   = round(sum(w)^2 / sum(w^2))
  ))
}

write.csv(respondent_groups32, file.path(output_dir, "private_04_respondents32.csv"),
          row.names = FALSE)
write.csv(outcome_summary, file.path(output_dir, "private_04_outcome32.csv"),
          row.names = FALSE)

# --- Disclosure check: Table 2a ----------------------------------------------
# Primary: counts 1-9. Secondary: if exactly one group is suppressed, also
# suppress the smallest other non-zero group (the respondent total is known).
checked <- respondent_groups32
primary <- checked$n > 0 & checked$n < disclosure_threshold
secondary <- rep(FALSE, nrow(checked))
if (sum(primary) == 1) {
  candidates <- which(!primary & checked$n > 0)
  if (length(candidates) > 0) secondary[candidates[which.min(checked$n[candidates])]] <- TRUE
}
checked$n[primary | secondary]       <- NA
checked$percent[primary | secondary] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "04_respondents32.csv"), row.names = FALSE)
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "04_respondents32", primary_suppressed = sum(primary),
  secondary_suppressed = sum(secondary)))
# -----------------------------------------------------------------------------

# --- Disclosure check: Table 2b/2c -------------------------------------------
# A row is shown only if both events and non-events are at least 10. Its
# percentages and effective n are built from the same people.
checked <- outcome_summary
too_small <- checked$events < disclosure_threshold |
             checked$non_events < disclosure_threshold
checked[too_small, c("n_observed", "events", "non_events", "percent_unweighted",
                     "percent_weighted", "kish_effective_n")] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "04_outcome32.csv"), row.names = FALSE)
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "04_outcome32", primary_suppressed = sum(too_small),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------


###############################################################################
# Block 5. Table 3: missingness state of every analysis variable
#
# Counts and percentages of the four states for each state_* variable, in:
#   - the target population (the imputation population; File 3);
#   - the prediction cohort (observed age-32 outcome with weight and design
#     information; defined in File 3).
###############################################################################

# --- Disclosure check: the target population's size --------------------------
# The target population is the original cohort minus the Wave-9 "Ineligible"
# group (File 3). Its total, shown in this table, therefore reveals the size
# of that group by subtraction from Table 1. If the group has 1-9 people,
# Table 1 hides it, so the script stops rather than reveal it here. (Table 3
# would then have to be reported for the original cohort instead.)
n_out_of_scope <- sum(!cohort$in_target_population)
if (n_out_of_scope > 0 && n_out_of_scope < disclosure_threshold) {
  stop("Fewer than 10 original-cohort members are out of scope (Wave-9 ",
       "ineligible). Table 3's target-population total would reveal this ",
       "group; report Table 3 for the original cohort instead.")
}
# -----------------------------------------------------------------------------

state_variables <- grep("^state_", names(cohort), value = TRUE)
samples <- list(
  "Target population" = cohort$in_target_population,
  "Prediction cohort" = cohort$in_prediction_cohort
)

missingness <- data.frame()
for (sample_name in names(samples)) {
  in_sample <- samples[[sample_name]]
  for (state_name in state_variables) {
    state_values <- cohort[[state_name]][in_sample]
    for (level in state_levels) {
      missingness <- rbind(missingness, data.frame(
        sample   = sample_name,
        variable = sub("^state_", "", state_name),
        state    = level,
        n        = sum(state_values == level),
        total    = length(state_values)
      ))
    }
  }
}
missingness$percent <- round(100 * missingness$n / missingness$total, 1)

write.csv(missingness, file.path(output_dir, "private_04_missingness.csv"),
          row.names = FALSE)

# --- Disclosure check: Table 3 -----------------------------------------------
# Primary: counts 1-9. Secondary: within each sample x variable, if exactly
# one state is suppressed, also suppress the next smallest non-zero state.
checked <- missingness
primary <- checked$n > 0 & checked$n < disclosure_threshold
secondary <- rep(FALSE, nrow(checked))
group <- paste(checked$sample, checked$variable)
for (g in unique(group)) {
  rows <- which(group == g)
  if (sum(primary[rows]) == 1) {
    candidates <- rows[!primary[rows] & checked$n[rows] > 0]
    if (length(candidates) > 0) {
      secondary[candidates[which.min(checked$n[candidates])]] <- TRUE
    }
  }
}
checked$n[primary | secondary]       <- NA
checked$percent[primary | secondary] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "04_missingness.csv"), row.names = FALSE)
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "04_missingness", primary_suppressed = sum(primary),
  secondary_suppressed = sum(secondary)))
# -----------------------------------------------------------------------------


###############################################################################
# Block 6. Table 4: prediction cohort and events at each landmark
#
# The prediction cohort (File 3) and its events are the same at every
# landmark. Only the features change, cumulatively.
# Complete-case counts (people observed on every feature up to that
# landmark) are reported for reference. The prediction analysis imputes
# missing features instead of dropping people.
###############################################################################

# prediction_cohort was defined in Block 4 from the shared flag in Block 2.
landmark_order <- c("L14", "L16", "L20", "L25a", "L25b")
landmark_summary <- data.frame()
for (i in seq_along(landmark_order)) {
  features <- landmark_features$variable[
    landmark_features$landmark %in% landmark_order[1:i]
  ]

  # A feature counts as observed when its state is "Observed" (features
  # without a state, such as months of observed activity, are always known).
  complete <- rep(TRUE, nrow(prediction_cohort))
  for (feature in features) {
    state_name <- paste0("state_", feature)
    if (state_name %in% names(prediction_cohort)) {
      complete <- complete & prediction_cohort[[state_name]] == "Observed"
    } else {
      complete <- complete & !is.na(prediction_cohort[[feature]])
    }
  }

  landmark_summary <- rbind(landmark_summary, data.frame(
    landmark              = landmark_order[i],
    features              = length(features),
    cohort_n              = nrow(prediction_cohort),
    events                = sum(prediction_cohort$sret_professional_32 == 1),
    non_events            = sum(prediction_cohort$sret_professional_32 == 0),
    complete_case_n       = sum(complete),
    complete_case_events  = sum(complete & prediction_cohort$sret_professional_32 == 1),
    complete_case_non_events = sum(complete & prediction_cohort$sret_professional_32 == 0)
  ))
}

write.csv(landmark_summary, file.path(output_dir, "private_04_landmarks.csv"),
          row.names = FALSE)

# --- Disclosure check: Table 4 -----------------------------------------------
# Cohort and complete-case counts are shown only when events and non-events
# are both at least 10. The complete-case n is suppressed together with its
# events and non-events, because it is their sum.
checked <- landmark_summary
cohort_small <- checked$events < disclosure_threshold |
                checked$non_events < disclosure_threshold
cc_small <- checked$complete_case_events < disclosure_threshold |
            checked$complete_case_non_events < disclosure_threshold
checked[cohort_small, c("cohort_n", "events", "non_events")] <- NA
checked[cc_small, c("complete_case_n", "complete_case_events",
                    "complete_case_non_events")] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "04_landmarks.csv"), row.names = FALSE)
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "04_landmarks", primary_suppressed = sum(cohort_small) + sum(cc_small),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------


###############################################################################
# Block 7. Table 5: survey weight distributions
#
# Design weight (original cohort) and Wave-9 final weight (Wave-9
# respondents with a weight). Reported: n, mean, SD, the 5th, 25th, 50th,
# 75th and 95th percentiles, Kish's effective sample size, and the design
# effect due to weighting (n / effective n).
###############################################################################

weight_summary <- data.frame()
for (weight_name in c("design_weight", "w9_final_weight")) {
  w <- cohort[[weight_name]]
  w <- w[!is.na(w)]
  percentiles <- quantile(w, c(0.05, 0.25, 0.50, 0.75, 0.95), names = FALSE)
  effective_n <- sum(w)^2 / sum(w^2)
  weight_summary <- rbind(weight_summary, data.frame(
    weight = weight_name, n = length(w),
    mean = round(mean(w), 3), sd = round(sd(w), 3),
    p5 = round(percentiles[1], 3), p25 = round(percentiles[2], 3),
    median = round(percentiles[3], 3), p75 = round(percentiles[4], 3),
    p95 = round(percentiles[5], 3),
    kish_effective_n = round(effective_n),
    design_effect_weighting = round(length(w) / effective_n, 2)
  ))
}

write.csv(weight_summary, file.path(output_dir, "private_04_weights.csv"),
          row.names = FALSE)

# --- Disclosure check: Table 5 -----------------------------------------------
# No minima or maxima are computed. Percentiles, means and SDs are shown
# only for weights based on at least 10 people. They are rounded to three
# decimals.
checked <- weight_summary
too_small <- checked$n < disclosure_threshold
checked[too_small, setdiff(names(checked), "weight")] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "04_weights.csv"), row.names = FALSE)
disclosure_log <- rbind(disclosure_log, data.frame(
  table = "04_weights", primary_suppressed = sum(too_small),
  secondary_suppressed = 0L))
# -----------------------------------------------------------------------------


###############################################################################
# Block 8. Disclosure log
#
# The number of cells suppressed in each checked table. It contains no data
# values.
###############################################################################

disclosure_log$threshold <- disclosure_threshold
write.csv(disclosure_log, file.path(review_dir, "04_disclosure_log.csv"),
          row.names = FALSE)


###############################################################################
# Block 9. Summary
#
# Table names and suppression counts only.
###############################################################################

cat("\nFile 4 finished.\n")
cat("  Tables written (exact, private):        ", nrow(disclosure_log), "\n")
cat("  Tables written (checked, in review/):   ", nrow(disclosure_log), "\n")
cat("  Cells suppressed (primary / secondary): ",
    sum(disclosure_log$primary_suppressed), "/",
    sum(disclosure_log$secondary_suppressed), "\n")
cat("  Disclosure log: 04_disclosure_log.csv\n")
