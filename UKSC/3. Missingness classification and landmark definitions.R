###############################################################################
# 3. Missingness classification and landmark definitions
#
# Purpose
#   1. Classify every missing value of every analysis variable into one of
#      four states:
#        Observed                     a substantive value is present
#        Item non-response            the question applied and the sweep was
#                                     answered, but no valid value was
#                                     recorded (e.g. don't know, refused,
#                                     unresolved occupation code)
#        Wave non-response            no productive interview at the sweep
#                                     the value comes from
#        Structural (not applicable)  the value is undefined because of an
#                                     earlier answer (routing), e.g. no
#                                     degree subject for someone without a
#                                     degree
#      Item and wave non-response are imputed in Files 7-13 and 17.
#      Structural values are never imputed.
#   2. Define the samples used later: the target population (for multiple
#      imputation) and the prediction cohort (observed age-32 outcome with
#      weight and design information).
#   3. Define the prediction landmarks (ages 14, 16, 20 and 25): which
#      features, and which imputation auxiliaries, are available at each.
#
# Rules
#   - "Structural" is assigned only by the explicit routing rules in Block 5.
#     A deposited "not applicable" code on an ordinary question (e.g. sex,
#     income) is treated as item non-response, so it is imputed rather than
#     left as a hole in the analysis.
#   - Categories that are themselves meaningful count as observed (e.g.
#     family NS-SEC "not applicable", parent education "no applicable
#     resident-parent record").
#   - Deposited missing-value codes mean different things by sweep (Waves
#     1-7: -1 = don't know, -91 = not applicable; Waves 8-9 derived
#     variables: -1 = not applicable, -8 = insufficient information,
#     -9 = refused). The derived variables from File 2 are already NA for
#     every non-substantive code, so the states here are defined from the
#     derived values, the response indicators and the routing rules.
#
# Reads
#   data_dir/02_constructed.rds
#
# Writes
#   data_dir/03_missingness_and_landmarks.rds  list: cohort (with state_*
#       variables and sample flags), state rules, landmark features,
#       landmark auxiliaries (participant level)
#   review_dir/03_missingness_rules.csv        how each variable's state is
#                                              defined (metadata only)
#   review_dir/03_landmark_features.csv        features by landmark (names only)
#   review_dir/03_landmark_auxiliaries.csv     auxiliaries by landmark (names
#                                              only)
#   Counts of each state are produced in File 4, after disclosure checks.
###############################################################################

library(haven)   # needed to read the deposited (SPSS-labelled) variables
suppressPackageStartupMessages(library(dplyr))


###############################################################################
# Block 1. Paths#
# Folder locations are set once, in File 0, and read here from
# project_paths.rds. Run this script with the working directory set to the
# folder that holds the scripts.
###############################################################################

paths      <- readRDS("project_paths.rds")
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs


###############################################################################
# Block 2. Read the constructed cohort from File 2
###############################################################################

cohort <- readRDS(file.path(data_dir, "02_constructed.rds"))

state_levels <- c("Observed", "Item non-response", "Wave non-response",
                  "Structural (not applicable)")


###############################################################################
# Block 3. Samples
#
# in_target_population: members of the original cohort who are in the
#   population the age-32 inferences refer to. This is the population for
#   the full-cohort multiple imputation (Files 5-13).
#   - Out of scope: Wave-9 response outcome 4, "ineligible" (found during
#     fieldwork to have died, or to be in prison or on probation). The Age 32
#     Sweep User Guide (section 7.2.2) defines the target population of the
#     Wave-9 weight as cohort members alive and not in prison or on
#     probation at age 32.
#   - In scope, as Wave-9 non-response: everyone else, including outcome -1,
#     "not issued". That group was not sent to fieldwork because members had
#     asked not to be contacted, were permanently untraced, or were known to
#     have died or to be in prison or on probation (User Guide, section
#     3.1). Refusal and loss of contact are non-response, not reasons to
#     leave the target population (Silverwood et al. 2024). The deposited
#     code doesn't say which reason applies, so the few known deaths and
#     prison cases in this group can't be removed; Silverwood et al. faced
#     the same limitation and expected it to matter little in a young
#     cohort.
# outcome_observed32: a productive Wave-9 interview with an observed age-32
#   SRET outcome.
# in_prediction_cohort: outcome_observed32, plus a positive Wave-9 final
#   weight and a known PSU and stratum. This single definition is used by
#   File 4's tables and by the prediction analysis (File 14 onward).
#
# The Wave-9 outcome labels are checked first, so a different coding scheme
# stops the script instead of silently defining the wrong sample.
###############################################################################

w9_outcome_labels <- attr(cohort$W9OUTCOME, "labels")
label_for_minus_1 <- tolower(names(w9_outcome_labels)[w9_outcome_labels == -1])
label_for_4       <- tolower(names(w9_outcome_labels)[w9_outcome_labels == 4])

if (!identical(label_for_minus_1, "not issued") || !identical(label_for_4, "ineligible")) {
  stop("W9OUTCOME doesn't use the expected codes (-1 = Not issued, ",
       "4 = Ineligible). Check its value labels before defining the sample.")
}

w9_outcome_code <- as.numeric(cohort$W9OUTCOME)
cohort$in_target_population <- !(w9_outcome_code %in% 4)
cohort$outcome_observed32    <- cohort$responded_w9 & !is.na(cohort$sret_professional_32)
cohort$in_prediction_cohort  <- cohort$outcome_observed32 &
  !is.na(cohort$w9_final_weight) & !is.na(cohort$psu) & !is.na(cohort$stratum)


###############################################################################
# Block 4. Variables measured at a single sweep, with no routing
#
# For these, the state is:
#   Observed           if the derived value is present
#   Wave non-response  if there was no productive interview at that sweep
#   Item non-response  otherwise
# One row per variable: the variable and the sweep it comes from.
###############################################################################

single_sweep_variables <- data.frame(
  variable = c(
    # Age 14 (Wave 1): background and orientation
    "sex14", "ethnicity14", "home_language14", "nssec_family14",
    "highest_parent_education14", "household_income14", "single_parent14",
    "independent_school14",
    "science_liking14", "science_competence14", "science_favourite14",
    "math_liking14", "math_competence14",
    # Age 15 (Wave 2) and age 16 (Wave 3)
    "science_course_route15",
    "w3_continue_fte16", "he_likelihood16",
    # Age 25 (Wave 8)
    "degree_attained25", "activity25", "employed25", "sret_professional_25",
    # Age 32 (Wave 9)
    "employed32", "sret_professional_32",
    # CLS predictors of non-response
    "nr_parents_know_where14", "nr_upset_by_names14", "nr_computer_games14",
    "nr_played_instrument14", "nr_housing_tenure14", "nr_internet_at_home14",
    "nr_school_contacted_parents15", "nr_under_strain15",
    "nr_smokes16", "nr_parent_age16",
    "nr_nightclubs17",
    "nr_same_address18", "nr_groups_treated_better18",
    "nr_teacher_expectations18", "nr_main_activity18",
    "nr_spoke_to_teacher19", "nr_willing_sex_questions19",
    "nr_willing_sex_questions20"
  ),
  sweep = c(
    rep(1, 13),
    2,
    3, 3,
    8, 8, 8, 8,
    9, 9,
    1, 1, 1, 1, 1, 1,
    2, 2,
    3, 3,
    4,
    5, 5, 5, 5,
    6, 6,
    7
  )
)

for (i in seq_len(nrow(single_sweep_variables))) {
  variable  <- single_sweep_variables$variable[i]
  responded <- cohort[[paste0("responded_w", single_sweep_variables$sweep[i])]]

  cohort[[paste0("state_", variable)]] <- factor(
    case_when(
      !is.na(cohort[[variable]]) ~ "Observed",
      !responded                 ~ "Wave non-response",
      TRUE                       ~ "Item non-response"
    ),
    levels = state_levels
  )
}


###############################################################################
# Block 5. Routed variables (the only source of structural states)
#
# Each part is defined only where an earlier answer makes it applicable:
#   w3_route_stage16       only for those staying on in full-time education
#                          at Wave 3
#   w3_intended_domain16   only where a subject area was given at Wave 3
#   he_domain20            only for those in HE by age 20
#   degree_domain25        only for degree holders at age 25
#   sret_among_employed32  only for those in paid work at age 32
# If the earlier answer is itself missing, applicability is unknown. The part
# is then item or wave non-response, and the imputation of the earlier
# answer decides it (File 7).
###############################################################################

# 5a. HE route at 16 (Wave 3), for those staying on in full-time education.
# An observed answer of leaving makes it not applicable, whatever the
# interview status, so that rule comes before wave non-response.
cohort$state_w3_route_stage16 <- factor(
  case_when(
    !is.na(cohort$w3_route_stage16)                                     ~ "Observed",
    cohort$w3_continue_fte16 %in% "Leave or not immediately return"     ~ "Structural (not applicable)",
    !cohort$responded_w3                                                ~ "Wave non-response",
    TRUE                                                                ~ "Item non-response"
  ),
  levels = state_levels
)

# 5b. Intended HE subject domain at 16 (Wave 3): not applicable for those
# leaving full-time education, or on a route other than subject given
# (again before wave non-response, for the same reason).
cohort$state_w3_intended_domain16 <- factor(
  case_when(
    !is.na(cohort$w3_intended_domain16)                                 ~ "Observed",
    cohort$w3_continue_fte16 %in% "Leave or not immediately return" |
      (!is.na(cohort$w3_route_stage16) &
         cohort$w3_route_stage16 != "Subject area given")               ~ "Structural (not applicable)",
    !cohort$responded_w3                                                ~ "Wave non-response",
    TRUE                                                                ~ "Item non-response"
  ),
  levels = state_levels
)

# 5c. HE participation (Waves 6 and 7) and HE subject domain at 20.
# Participation needs both sweeps unless HE was reported at either one. If
# it is unresolved and either sweep lacks an interview, that is wave
# non-response.
cohort$state_he_participation20 <- factor(
  case_when(
    !is.na(cohort$he_participation20)             ~ "Observed",
    !cohort$responded_w6 | !cohort$responded_w7   ~ "Wave non-response",
    TRUE                                          ~ "Item non-response"
  ),
  levels = state_levels
)

cohort$state_he_domain20 <- factor(
  case_when(
    !is.na(cohort$he_domain20)                                    ~ "Observed",
    cohort$he_participation20 %in% "Not in HE by age 20"          ~ "Structural (not applicable)",
    !cohort$responded_w6 & !cohort$responded_w7                   ~ "Wave non-response",
    is.na(cohort$he_participation20) &
      (!cohort$responded_w6 | !cohort$responded_w7)               ~ "Wave non-response",
    TRUE                                                          ~ "Item non-response"
  ),
  levels = state_levels
)

# 5d. Degree subject domain at 25 (Wave 8).
cohort$state_degree_domain25 <- factor(
  case_when(
    !is.na(cohort$degree_domain25)                    ~ "Observed",
    !cohort$responded_w8                              ~ "Wave non-response",
    cohort$degree_attained25 %in% "No degree"         ~ "Structural (not applicable)",
    TRUE                                              ~ "Item non-response"
  ),
  levels = state_levels
)

# 5e. SRET occupation among the employed at 32 (Wave 9). Item non-response
# here includes someone in work whose occupation couldn't be coded.
cohort$state_sret_among_employed32 <- factor(
  case_when(
    !is.na(cohort$sret_among_employed32)              ~ "Observed",
    !cohort$responded_w9                              ~ "Wave non-response",
    cohort$employed32 %in% "Not in paid work"         ~ "Structural (not applicable)",
    TRUE                                              ~ "Item non-response"
  ),
  levels = state_levels
)


###############################################################################
# Block 6. Combined pathway variables
#
# Each combination (continuation x route x domain, participation x domain,
# attained x domain) is observed when it has a value. Otherwise it takes the
# state of whichever part is missing, with wave non-response taking
# precedence over item non-response. The combinations are rebuilt from
# their imputed parts in File 7.
###############################################################################

cohort$state_w3_intention_route16 <- factor(
  case_when(
    !is.na(cohort$w3_intention_route16)                                  ~ "Observed",
    cohort$state_w3_continue_fte16 == "Wave non-response" |
      cohort$state_w3_route_stage16 == "Wave non-response" |
      cohort$state_w3_intended_domain16 == "Wave non-response"           ~ "Wave non-response",
    TRUE                                                                 ~ "Item non-response"
  ),
  levels = state_levels
)

cohort$state_he_pathway20 <- factor(
  case_when(
    !is.na(cohort$he_pathway20)                                          ~ "Observed",
    cohort$state_he_participation20 == "Wave non-response" |
      cohort$state_he_domain20 == "Wave non-response"                    ~ "Wave non-response",
    TRUE                                                                 ~ "Item non-response"
  ),
  levels = state_levels
)

cohort$state_degree_pathway25 <- factor(
  case_when(
    !is.na(cohort$degree_pathway25)                                      ~ "Observed",
    cohort$state_degree_attained25 == "Wave non-response" |
      cohort$state_degree_domain25 == "Wave non-response"                ~ "Wave non-response",
    TRUE                                                                 ~ "Item non-response"
  ),
  levels = state_levels
)


###############################################################################
# Block 7. Monthly activity histories, ages 16-19
#
# The activity counts are complete only when every month has a substantive
# code:
#   Observed           all months observed
#   Item non-response  some months observed (the history stops early)
#   Wave non-response  no months observed
# The same state applies to every activity count and to the transition
# count. months_activity_observed_16_19 itself is always known, so it is
# never imputed (it is an L20 prediction feature).
###############################################################################

cohort$state_activity_16_19 <- factor(
  case_when(
    cohort$months_activity_observed_16_19 == cohort$activity_months_expected ~ "Observed",
    cohort$months_activity_observed_16_19 > 0                               ~ "Item non-response",
    TRUE                                                                    ~ "Wave non-response"
  ),
  levels = state_levels
)

activity_count_variables <- c("months_education_16_19", "months_employed_16_19",
                              "months_apprenticeship_16_19", "months_neet_16_19",
                              "activity_transitions_16_19")
for (variable in activity_count_variables) {
  cohort[[paste0("state_", variable)]] <- cohort$state_activity_16_19
}


###############################################################################
# Block 8. Record the state rules (metadata)
###############################################################################

state_rules <- rbind(
  data.frame(
    variable = single_sweep_variables$variable,
    rule = paste0("Single sweep (Wave ", single_sweep_variables$sweep,
                  "): observed / wave non-response / item non-response"),
    structural_when = "never"
  ),
  data.frame(
    variable = c("w3_route_stage16", "w3_intended_domain16", "he_participation20",
                 "he_domain20", "degree_domain25", "sret_among_employed32"),
    rule = c("Wave 3, routed on w3_continue_fte16",
             "Wave 3, routed on w3_continue_fte16 and w3_route_stage16",
             "Waves 6-7: resolved if in HE at either sweep, or not in HE at both",
             "Waves 6-7, routed on he_participation20",
             "Wave 8, routed on degree_attained25",
             "Wave 9, routed on employed32"),
    structural_when = c("w3_continue_fte16 is 'Leave or not immediately return'",
                        "leaving full-time education, or w3_route_stage16 is not 'Subject area given'",
                        "never",
                        "he_participation20 is 'Not in HE by age 20'",
                        "degree_attained25 is 'No degree'",
                        "employed32 is 'Not in paid work'")
  ),
  data.frame(
    variable = c("w3_intention_route16", "he_pathway20", "degree_pathway25"),
    rule = "Combination of its parts: takes the state of the missing part",
    structural_when = "never"
  ),
  data.frame(
    variable = c("activity_16_19", activity_count_variables),
    rule = "Monthly coverage: all months / some months / no months observed",
    structural_when = "never"
  )
)
state_rules$disclosure_checked <- TRUE


###############################################################################
# Block 9. Prediction landmarks: features
#
# Cumulative feature sets. Each landmark adds the information that becomes
# available at that age, and keeps everything from earlier landmarks:
#   L14   age-14 background and science/mathematics orientation
#   L16   + age-15 course route and age-16 continuation, HE likelihood and
#           intended-subject route
#   L20   + HE pathway by 20 and activity histories 16-19
#   L25a  + degree pathway and main activity at 25
#   L25b  + SRET occupation at 25 (a sensitivity analysis: it adds
#           occupational persistence rather than developmental information)
# Combined pathway variables are features. In the imputation they are
# rebuilt from their imputed parts.
###############################################################################

landmark_features <- rbind(
  data.frame(landmark = "L14", available_by_age = 14,
             block = "Background",
             variable = c("sex14", "ethnicity14", "home_language14",
                          "nssec_family14", "highest_parent_education14",
                          "household_income14", "single_parent14",
                          "independent_school14")),
  data.frame(landmark = "L14", available_by_age = 14,
             block = "Science and mathematics orientation",
             variable = c("science_liking14", "science_competence14",
                          "science_favourite14", "math_liking14",
                          "math_competence14")),
  data.frame(landmark = "L16", available_by_age = 16,
             block = "Course route and intentions",
             variable = c("science_course_route15", "w3_continue_fte16",
                          "he_likelihood16", "w3_intention_route16")),
  data.frame(landmark = "L20", available_by_age = 20,
             block = "HE and early activity",
             variable = c("he_pathway20", "months_employed_16_19",
                          "months_apprenticeship_16_19", "months_neet_16_19",
                          "activity_transitions_16_19",
                          "months_activity_observed_16_19")),
  data.frame(landmark = "L25a", available_by_age = 25,
             block = "Attained education and activity",
             variable = c("degree_pathway25", "activity25")),
  data.frame(landmark = "L25b", available_by_age = 25,
             block = "Occupation at 25 (sensitivity)",
             variable = "sret_professional_25")
)

# The parts that are imputed for each combined feature.
landmark_features$imputed_through <- case_when(
  landmark_features$variable == "w3_intention_route16" ~ "w3_continue_fte16 + w3_route_stage16 + w3_intended_domain16",
  landmark_features$variable == "he_pathway20"         ~ "he_participation20 + he_domain20",
  landmark_features$variable == "degree_pathway25"     ~ "degree_attained25 + degree_domain25",
  TRUE                                                 ~ landmark_features$variable
)
landmark_features$disclosure_checked <- TRUE


###############################################################################
# Block 10. Prediction landmarks: imputation auxiliaries
#
# Variables that may help impute missing features but are not model
# features. They are the CLS predictors of non-response (user guide
# Appendix C, Table C1), used only from the age each is known:
#   - the individual predictors, by the age they were measured;
#   - previous non-response (one or more sweeps without an interview),
#     defined as in Appendix C but over the sweeps known by the landmark:
#       L16   Waves 2-3
#       L20   Waves 2-7 (the published Appendix C item)
#       L25   Waves 2-8
#     Each landmark uses only its own span, which replaces the shorter one.
# Wave-9 response is never used: it is future information at every
# landmark. Sex is a model feature, not an auxiliary. Candidates are
# screened inside each training fold (File 17).
###############################################################################

landmark_order <- c("L14", "L16", "L20", "L25a", "L25b")
landmark_age   <- c(L14 = 14, L16 = 16, L20 = 20, L25a = 25, L25b = 25)

# known_by_age:    the age from which the auxiliary is available.
# replaced_by_age: for the previous non-response items, the age from which a
#                  longer span replaces it (NA = never replaced).
auxiliary_candidates <- data.frame(
  variable = c(
    "nr_parents_know_where14", "nr_upset_by_names14", "nr_computer_games14",
    "nr_played_instrument14", "nr_housing_tenure14", "nr_internet_at_home14",
    "nr_school_contacted_parents15", "nr_under_strain15",
    "nr_smokes16", "nr_parent_age16",
    "nr_nightclubs17",
    "nr_same_address18", "nr_groups_treated_better18",
    "nr_teacher_expectations18", "nr_main_activity18",
    "nr_spoke_to_teacher19", "nr_willing_sex_questions19",
    "nr_willing_sex_questions20",
    "nr_previous_nonresponse_w2_w3",
    "nr_previous_nonresponse_w2_w7",
    "nr_previous_nonresponse_w2_w8"
  ),
  known_by_age = c(
    14, 14, 14, 14, 14, 14,
    15, 15,
    16, 16,
    17,
    18, 18, 18, 18,
    19, 19,
    20,
    16, 20, 25
  ),
  replaced_by_age = c(
    rep(NA, 18),
    20, 25, NA
  )
)

landmark_auxiliaries <- data.frame()
for (landmark in landmark_order) {
  age <- landmark_age[[landmark]]
  available <- auxiliary_candidates$known_by_age <= age &
    (is.na(auxiliary_candidates$replaced_by_age) |
       age < auxiliary_candidates$replaced_by_age)
  landmark_auxiliaries <- rbind(landmark_auxiliaries, data.frame(
    landmark     = landmark,
    variable     = auxiliary_candidates$variable[available],
    known_by_age = auxiliary_candidates$known_by_age[available]
  ))
}
landmark_auxiliaries$disclosure_checked <- TRUE


###############################################################################
# Block 11. Checks
#
# Each check stops the script if it fails. No data values are printed.
###############################################################################

state_variables <- grep("^state_", names(cohort), value = TRUE)

# Every state variable is complete (no NA states).
for (state_name in state_variables) {
  if (anyNA(cohort[[state_name]])) stop(state_name, " has missing states.")
}

# "Observed" exactly when the value is present, for every classified
# variable except the activity summary state.
for (state_name in setdiff(state_variables, "state_activity_16_19")) {
  variable <- sub("^state_", "", state_name)
  if (variable %in% activity_count_variables) next
  value_present <- !is.na(cohort[[variable]])
  if (any(value_present != (cohort[[state_name]] == "Observed"))) {
    stop(state_name, " disagrees with whether ", variable, " is present.")
  }
}

# Structural states only occur in the five routed parts.
routed_parts <- c("state_w3_route_stage16", "state_w3_intended_domain16",
                  "state_he_domain20", "state_degree_domain25",
                  "state_sret_among_employed32")
for (state_name in setdiff(state_variables, routed_parts)) {
  if (any(cohort[[state_name]] == "Structural (not applicable)")) {
    stop(state_name, " has structural states but isn't a routed variable.")
  }
}

# Every landmark feature and auxiliary exists in the data.
missing_names <- setdiff(c(landmark_features$variable, landmark_auxiliaries$variable),
                         names(cohort))
if (length(missing_names) > 0) {
  stop("Landmark variables not found: ", paste(missing_names, collapse = ", "))
}

# No auxiliary may be used at a landmark before the age it becomes known.
if (any(landmark_auxiliaries$known_by_age > landmark_age[landmark_auxiliaries$landmark])) {
  stop("An auxiliary is used before it is known.")
}


###############################################################################
# Block 12. Save
###############################################################################

saveRDS(
  list(
    cohort               = cohort,
    state_levels         = state_levels,
    state_rules          = state_rules,
    landmark_features    = landmark_features,
    landmark_auxiliaries = landmark_auxiliaries
  ),
  file.path(data_dir, "03_missingness_and_landmarks.rds"),
  compress = "xz"
)

# Review tables. Disclosure check: variable names, rules and landmark
# assignments only; no participant data and no counts.
write.csv(state_rules, file.path(review_dir, "03_missingness_rules.csv"),
          row.names = FALSE)
write.csv(landmark_features, file.path(review_dir, "03_landmark_features.csv"),
          row.names = FALSE)
write.csv(landmark_auxiliaries, file.path(review_dir, "03_landmark_auxiliaries.csv"),
          row.names = FALSE)


###############################################################################
# Block 13. Summary
#
# Counts of variables and checks only. State counts are in File 4.
###############################################################################

cat("\nFile 3 finished.\n")
cat("  Variables with a missingness state:  ", length(state_variables), "\n")
cat("  Routed variables (structural states):", length(routed_parts), "\n")
cat("  Landmark features (all landmarks):   ", length(unique(landmark_features$variable)), "\n")
cat("  Auxiliary candidates:                ", nrow(auxiliary_candidates), "\n")
cat("  Checks:                               all passed\n")
cat("  Saved:", file.path(data_dir, "03_missingness_and_landmarks.rds"), "\n")
