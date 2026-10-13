###############################################################################
# 14. Prediction cohort and fold design
#
# Purpose
#   Set up the prediction analysis (Files 15-19):
#     1. the prediction cohort;
#     2. the variables available at each prediction landmark: the model
#        features, the routed parts imputed for combined features, and the
#        imputation auxiliaries;
#     3. the cross-validation folds: 5 repeats of 5 outer folds, and 5 inner
#        tuning folds within every outer training set.
#   This script fits no model and computes no performance measure.
#
# Study population and aim (REFORMS 1a-c)
#   The prediction analysis asks how well information available at each age
#   distinguishes people who are in professional science, research,
#   engineering and technology (SRET) employment at 32 from those who are
#   not, and how much each developmental stage adds. It is a study of
#   predictability. The models are not intended to screen or select
#   individuals, and none of the associations is causal.
#   Performance can only be assessed where the outcome is known, so the
#   prediction cohort is the people with a productive Wave-9 interview, an
#   observed age-32 SRET outcome, and a Wave-9 final weight, PSU and stratum
#   (the single definition made in File 3). Performance is weighted by the
#   Wave-9 final weight, so it refers to the population that the Wave-9
#   respondents represent after non-response weighting. Imputed outcomes are
#   never used as the truth when scoring predictions. How this cohort
#   relates to the full cohort is shown by the multiple imputation analysis
#   (File 12, observed-outcome sample).
#
# Landmarks (defined in File 3)
#   Cumulative sets: each landmark adds what becomes known at that age and
#   keeps everything from earlier landmarks.
#     L14   background and science/mathematics orientation
#     L16   + course route (15), continuation, HE likelihood and intended
#             subject route (16)
#     L20   + HE pathway and activity histories 16-19
#     L25a  + degree pathway and main activity at 25
#     L25b  + SRET occupation at 25 (a sensitivity analysis)
#   Combined features (intended subject route, HE pathway, degree pathway)
#   are imputed through their routed parts, as in the multiple imputation
#   (File 7). Auxiliaries are the CLS predictors of non-response known by
#   each landmark (File 3).
#   Candidate parameters per landmark are counted for the sample size
#   calculation (File 15): one per numeric feature, and one per category
#   beyond the first for categorical features.
#
# Cross-validation folds (REFORMS 5c, 6; TRIPOD 12a)
#   - Next Steps sampled schools (PSUs), and pupils in the same school share
#     a context. If people from one PSU were split between training and test
#     folds, performance could be overstated, so every PSU is kept whole
#     within one fold.
#   - Outer folds: 5 repeats of 5 folds. Repeating the split with different
#     random assignments reduces the dependence of results on one split.
#   - Assignment, for each repeat: shuffle the PSUs using the repeat's own
#     seed, sort them by their number of SRET events, and deal them to the
#     folds back and forth (1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 1, 2, ...). This
#     balances events across folds. PSUs with the same number of events keep
#     their shuffled order, so the repeats differ.
#   - Inner folds: within each outer training set, the training PSUs are
#     split into 5 inner folds in the same way, for tuning (File 17).
#   - The same outer and inner folds are used at every landmark, so
#     landmarks are compared on the same people (paired comparisons).
#   - Everything estimated from data (imputation, pooling of rare
#     categories, scaling, tuning) is fitted within training folds in
#     File 17, never on the whole cohort.
#
# Disclosure control (UK Data Service output rules; threshold 10)
#   The fold balance table shows a fold only if it has at least 10 events
#   and 10 non-events, with secondary suppression within each repeat (the
#   folds of a repeat add up to the whole cohort). Its groups use the same
#   prediction cohort definition as File 4, so no small group can be
#   recovered by subtraction. The other review tables hold names, category
#   counts of variables (not of people) and pass/fail checks only.
#
# Reads
#   data_dir/03_missingness_and_landmarks.rds
#
# Writes
#   data_dir/14_prediction_design.rds       prediction cohort with folds,
#                                           landmark variables, settings
#                                           (participant level)
#   output_dir/private_14_fold_balance.csv  exact fold table
#   review_dir/14_landmark_design.csv       features, parts, auxiliaries and
#                                           candidate parameters by landmark
#   review_dir/14_landmark_variables.csv    every variable by landmark and
#                                           role (names only)
#   review_dir/14_fold_balance.csv          people, PSUs and events by fold
#                                           (checked)
#   review_dir/14_fold_checks.csv           pass/fail checks
#   review_dir/14_disclosure_log.csv        what was suppressed
#
# References
#   Kapoor S, et al. (2024). REFORMS: Consensus-based recommendations for
#   machine-learning-based science. Science Advances 10(18): eadk3452.
#   Collins GS, et al. (2024). TRIPOD+AI statement. BMJ 385: e078378.
###############################################################################

library(haven)   # needed to read the deposited (SPSS-labelled) variables
suppressPackageStartupMessages(library(dplyr))


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

n_repeats     <- 5      # repeats of the outer cross-validation
n_outer_folds <- 5      # outer folds per repeat
n_inner_folds <- 5      # inner tuning folds per outer training set
fold_seed     <- 1414   # outer folds use fold_seed + repeat;
                        # inner folds use fold_seed + 100 * repeat + fold

landmark_order <- c("L14", "L16", "L20", "L25a", "L25b")
landmark_age   <- c(L14 = 14, L16 = 16, L20 = 20, L25a = 25, L25b = 25)
outcome        <- "sret_professional_32"


###############################################################################
# Block 2. Read the cohort and landmark definitions from File 3
###############################################################################

file3                <- readRDS(file.path(data_dir, "03_missingness_and_landmarks.rds"))
cohort               <- file3$cohort
landmark_features    <- file3$landmark_features
landmark_auxiliaries <- file3$landmark_auxiliaries


###############################################################################
# Block 3. Variables at each landmark
#
# One row per landmark x variable, with its role:
#   Feature            enters the prediction models
#   Part of <feature>  a routed part imputed for a combined feature (the
#                      feature is rebuilt from its parts after imputation)
#   Auxiliary          helps the imputation only; not a model feature
###############################################################################

landmark_variables <- data.frame()

for (k in seq_along(landmark_order)) {
  landmark <- landmark_order[k]

  # Features introduced at this landmark or an earlier one.
  features <- landmark_features[match(landmark_features$landmark, landmark_order) <= k, ]
  landmark_variables <- rbind(landmark_variables, data.frame(
    landmark = landmark, variable = features$variable, role = "Feature",
    block = features$block, available_by_age = features$available_by_age
  ))

  # Routed parts of the combined features.
  combined <- features[features$imputed_through != features$variable, ]
  for (i in seq_len(nrow(combined))) {
    parts <- strsplit(combined$imputed_through[i], " + ", fixed = TRUE)[[1]]
    landmark_variables <- rbind(landmark_variables, data.frame(
      landmark = landmark, variable = parts,
      role = paste("Part of", combined$variable[i]),
      block = combined$block[i], available_by_age = combined$available_by_age[i]
    ))
  }

  # Auxiliaries known by this landmark.
  auxiliaries <- landmark_auxiliaries[landmark_auxiliaries$landmark == landmark, ]
  landmark_variables <- rbind(landmark_variables, data.frame(
    landmark = landmark, variable = auxiliaries$variable, role = "Auxiliary",
    block = "CLS predictor of non-response",
    available_by_age = auxiliaries$known_by_age
  ))
}

# Checks (each stops the script if it fails):
# every variable is known by its landmark's age;
if (any(landmark_variables$available_by_age > landmark_age[landmark_variables$landmark])) {
  stop("A variable is used at a landmark before it is known.")
}
# no age-32 or Wave-9 information enters any landmark (it is the outcome or
# comes after every landmark);
future_information <- c(outcome, "employed32", "sret_among_employed32",
                        "responded_w9", "w9_final_weight", "outcome_observed32",
                        "in_prediction_cohort")
if (any(landmark_variables$variable %in% future_information)) {
  stop("Age-32 or Wave-9 information is used at a landmark.")
}
# every variable exists in the data.
not_found <- setdiff(landmark_variables$variable, names(cohort))
if (length(not_found) > 0) {
  stop("Landmark variables not found: ", paste(not_found, collapse = ", "))
}


###############################################################################
# Block 4. The prediction cohort
#
# People with in_prediction_cohort = TRUE (File 3). Only the variables the
# prediction analysis needs are kept: identifiers and design, the outcome,
# every landmark variable, and the missingness state of each (File 17 uses
# the states to impute item and wave non-response and never structural
# values).
###############################################################################

prediction <- cohort[cohort$in_prediction_cohort, ]

# Checks: the cohort definition holds (each stops the script if it fails).
if (anyNA(prediction[[outcome]]) || any(!(prediction[[outcome]] %in% c(0, 1)))) {
  stop("The outcome is missing or not 0/1 in the prediction cohort.")
}
if (anyNA(prediction$w9_final_weight) || any(prediction$w9_final_weight <= 0)) {
  stop("A Wave-9 final weight is missing or not positive.")
}
if (anyNA(prediction$psu) || anyNA(prediction$stratum)) {
  stop("A PSU or stratum is missing in the prediction cohort.")
}

all_landmark_variables <- unique(landmark_variables$variable)
state_columns <- intersect(paste0("state_", all_landmark_variables), names(prediction))

prediction <- as.data.frame(prediction[, c("NSID", "psu", "stratum", "w9_final_weight",
                                           outcome, all_landmark_variables,
                                           state_columns)])
prediction$psu <- droplevels(prediction$psu)
prediction$stratum <- droplevels(prediction$stratum)


###############################################################################
# Block 5. Landmark summary and candidate parameters
#
# For each landmark: numbers of features, routed parts and auxiliaries, and
# the candidate parameters for File 15 (one per numeric feature; one per
# category beyond the first for categorical features, counting the
# categories present in the prediction cohort).
###############################################################################

landmark_summary <- data.frame()

for (landmark in landmark_order) {
  rows     <- landmark_variables[landmark_variables$landmark == landmark, ]
  features <- rows$variable[rows$role == "Feature"]

  candidate_parameters <- 0
  for (variable in features) {
    values <- prediction[[variable]]
    if (is.factor(values)) {
      candidate_parameters <- candidate_parameters + nlevels(droplevels(values)) - 1
    } else {
      candidate_parameters <- candidate_parameters + 1
    }
  }

  landmark_summary <- rbind(landmark_summary, data.frame(
    landmark             = landmark,
    age                  = landmark_age[[landmark]],
    features             = length(features),
    candidate_parameters = candidate_parameters,
    routed_parts         = sum(startsWith(rows$role, "Part of")),
    auxiliaries          = sum(rows$role == "Auxiliary"),
    feature_names        = paste(features, collapse = " | "),
    auxiliary_names      = paste(rows$variable[rows$role == "Auxiliary"], collapse = " | ")
  ))
}


###############################################################################
# Block 6. Outer folds
#
# One column per repeat (outer_fold_r1 ... outer_fold_r5): the outer fold
# (1-5) each person's PSU is assigned to. See the header for the method.
###############################################################################

# One row per PSU: its people and events.
psu_table <- data.frame(psu = levels(prediction$psu))
psu_table$people <- as.vector(table(prediction$psu))
psu_table$events <- as.vector(table(prediction$psu[prediction[[outcome]] == 1]))

deal_order_outer <- rep(c(1:n_outer_folds, n_outer_folds:1),
                        length.out = nrow(psu_table))

for (r in 1:n_repeats) {
  set.seed(fold_seed + r)
  shuffled <- psu_table[sample(nrow(psu_table)), ]
  dealt    <- shuffled[order(-shuffled$events), ]   # ties keep the shuffled order
  dealt$fold <- deal_order_outer
  prediction[[paste0("outer_fold_r", r)]] <-
    dealt$fold[match(as.character(prediction$psu), dealt$psu)]
}


###############################################################################
# Block 7. Inner folds
#
# One column per repeat x outer fold (inner_fold_r1_f1 ... inner_fold_r5_f5):
# the inner fold (1-5) of each person in that outer training set, and NA for
# the people in the outer test fold. The training PSUs are split with the
# same method as the outer folds.
###############################################################################

for (r in 1:n_repeats) {
  for (f in 1:n_outer_folds) {
    in_training   <- prediction[[paste0("outer_fold_r", r)]] != f
    training_psus <- psu_table[psu_table$psu %in% prediction$psu[in_training], ]

    set.seed(fold_seed + 100 * r + f)
    shuffled <- training_psus[sample(nrow(training_psus)), ]
    dealt    <- shuffled[order(-shuffled$events), ]
    dealt$fold <- rep(c(1:n_inner_folds, n_inner_folds:1), length.out = nrow(dealt))

    inner_fold <- dealt$fold[match(as.character(prediction$psu), dealt$psu)]
    inner_fold[!in_training] <- NA
    prediction[[paste0("inner_fold_r", r, "_f", f)]] <- inner_fold
  }
}


###############################################################################
# Block 8. Fold checks
#
#   a. every person has an outer fold in every repeat;
#   b. every PSU is in exactly one outer fold in each repeat;
#   c. every outer test fold has at least 10 events and 10 non-events;
#   d. inner folds cover exactly the outer training set, keep PSUs whole,
#      and every inner fold has at least 10 events and 10 non-events (the
#      minimum for tuning; the fold's counts are not shown);
#   e. no two repeats have the same assignment.
# The pass/fail table is written first; the script then stops if any check
# failed.
###############################################################################

y <- prediction[[outcome]]
fold_checks <- data.frame()

for (r in 1:n_repeats) {
  outer <- prediction[[paste0("outer_fold_r", r)]]

  # a and b
  every_person_assigned <- !anyNA(outer)
  psu_by_fold <- table(prediction$psu, outer)
  psus_whole  <- all(rowSums(psu_by_fold > 0) == 1)

  # c
  outer_events_ok <- TRUE
  for (f in 1:n_outer_folds) {
    if (sum(y[outer == f] == 1) < disclosure_threshold ||
        sum(y[outer == f] == 0) < disclosure_threshold) {
      outer_events_ok <- FALSE
    }
  }

  # d
  inner_ok <- TRUE
  for (f in 1:n_outer_folds) {
    inner <- prediction[[paste0("inner_fold_r", r, "_f", f)]]
    if (any(is.na(inner) != (outer == f))) inner_ok <- FALSE
    training <- !is.na(inner)
    psu_by_inner <- table(droplevels(prediction$psu[training]), inner[training])
    if (!all(rowSums(psu_by_inner > 0) == 1)) inner_ok <- FALSE
    for (g in 1:n_inner_folds) {
      too_few <- sum(y[training & inner == g] == 1) < disclosure_threshold ||
        sum(y[training & inner == g] == 0) < disclosure_threshold
      if (too_few) inner_ok <- FALSE
    }
  }

  fold_checks <- rbind(fold_checks, data.frame(
    repeat_number = r,
    check = c("Every person has an outer fold",
              "Every PSU is in one outer fold",
              "Every outer fold has at least 10 events and 10 non-events",
              "Inner folds cover the training set, keep PSUs whole, and have at least 10 events and 10 non-events"),
    passed = c(every_person_assigned, psus_whole, outer_events_ok, inner_ok)
  ))
}

# e
repeats_differ <- TRUE
for (r1 in 1:(n_repeats - 1)) {
  for (r2 in (r1 + 1):n_repeats) {
    if (identical(prediction[[paste0("outer_fold_r", r1)]],
                  prediction[[paste0("outer_fold_r", r2)]])) {
      repeats_differ <- FALSE
    }
  }
}
fold_checks <- rbind(fold_checks, data.frame(
  repeat_number = NA, check = "No two repeats have the same assignment",
  passed = repeats_differ
))

# Disclosure check: pass/fail results only; no counts.
fold_checks$disclosure_checked <- TRUE
write.csv(fold_checks, file.path(review_dir, "14_fold_checks.csv"),
          row.names = FALSE, na = "")

if (!all(fold_checks$passed)) {
  stop("A fold check failed; see 14_fold_checks.csv.")
}


###############################################################################
# Block 9. Fold balance table
#
# For each repeat x outer fold (the test fold): people, PSUs, events,
# non-events, and the unweighted and Wave-9-weighted SRET percentages.
###############################################################################

fold_balance <- data.frame()
for (r in 1:n_repeats) {
  outer <- prediction[[paste0("outer_fold_r", r)]]
  for (f in 1:n_outer_folds) {
    in_fold <- outer == f
    w <- prediction$w9_final_weight[in_fold]
    fold_balance <- rbind(fold_balance, data.frame(
      repeat_number = r, fold = f,
      people     = sum(in_fold),
      psus       = length(unique(prediction$psu[in_fold])),
      events     = sum(y[in_fold] == 1),
      non_events = sum(y[in_fold] == 0),
      percent_events          = round(100 * mean(y[in_fold]), 2),
      weighted_percent_events = round(100 * sum(w * y[in_fold]) / sum(w), 2)
    ))
  }
}

write.csv(fold_balance, file.path(output_dir, "private_14_fold_balance.csv"),
          row.names = FALSE)

# --- Disclosure check: fold balance ------------------------------------------
# Primary: a fold's counts and percentages are hidden unless it has at least
# 10 events and 10 non-events. Secondary: within each repeat, if exactly one
# fold is hidden, the fold with the next fewest events is hidden too (the
# folds of a repeat add up to the whole cohort).
checked <- fold_balance
primary <- checked$events < disclosure_threshold |
           checked$non_events < disclosure_threshold
secondary <- rep(FALSE, nrow(checked))
for (r in unique(checked$repeat_number)) {
  rows <- which(checked$repeat_number == r)
  if (sum(primary[rows]) == 1) {
    candidates <- rows[!primary[rows]]
    secondary[candidates[which.min(checked$events[candidates])]] <- TRUE
  }
}
checked[primary | secondary, c("people", "psus", "events", "non_events",
                               "percent_events", "weighted_percent_events")] <- NA
checked$disclosure_checked <- TRUE
write.csv(checked, file.path(review_dir, "14_fold_balance.csv"),
          row.names = FALSE, na = "")

disclosure_log <- data.frame(
  table = c("14_fold_balance", "14_fold_checks", "14_landmark_design",
            "14_landmark_variables"),
  primary_suppressed   = c(sum(primary), 0L, 0L, 0L),
  secondary_suppressed = c(sum(secondary), 0L, 0L, 0L),
  threshold            = disclosure_threshold
)
write.csv(disclosure_log, file.path(review_dir, "14_disclosure_log.csv"),
          row.names = FALSE)
# -----------------------------------------------------------------------------


###############################################################################
# Block 10. Landmark design tables (metadata)
#
# Disclosure check: variable names, roles, ages and parameter counts only.
# Parameter counts are numbers of categories of variables, not of people.
###############################################################################

landmark_summary$disclosure_checked <- TRUE
write.csv(landmark_summary, file.path(review_dir, "14_landmark_design.csv"),
          row.names = FALSE)

landmark_variables$disclosure_checked <- TRUE
write.csv(landmark_variables, file.path(review_dir, "14_landmark_variables.csv"),
          row.names = FALSE)


###############################################################################
# Block 11. Save
###############################################################################

saveRDS(
  list(
    cohort             = prediction,
    landmark_variables = landmark_variables,
    landmark_summary   = landmark_summary,
    settings           = list(outcome = outcome, landmark_order = landmark_order,
                              landmark_age = landmark_age, n_repeats = n_repeats,
                              n_outer_folds = n_outer_folds,
                              n_inner_folds = n_inner_folds, fold_seed = fold_seed)
  ),
  file.path(data_dir, "14_prediction_design.rds"),
  compress = "xz"
)


###############################################################################
# Block 12. Summary
#
# Design counts and checks only. People and events by fold are in the
# checked table 14_fold_balance.csv.
###############################################################################

cat("\nFile 14 finished.\n")
cat("  PSUs in the prediction cohort:        ", nrow(psu_table), "\n")
cat("  Outer folds:", n_repeats, "repeats x", n_outer_folds, "folds;",
    "inner folds per training set:", n_inner_folds, "\n")
cat("  Fold checks:                           all passed\n")
cat("  Folds suppressed (primary/secondary): ", sum(primary), "/", sum(secondary), "\n")
cat("  Landmarks (features / candidate parameters / auxiliaries):\n")
for (i in seq_len(nrow(landmark_summary))) {
  cat("    ", landmark_summary$landmark[i], ": ", landmark_summary$features[i], " / ",
      landmark_summary$candidate_parameters[i], " / ",
      landmark_summary$auxiliaries[i], "\n", sep = "")
}
cat("  Saved:", file.path(data_dir, "14_prediction_design.rds"), "\n")
