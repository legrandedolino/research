###############################################################################
# 2. Variable construction
#
# Purpose
#   Build every analysis variable from the deposited Next Steps variables.
#   Each derived variable has its own commented block showing its source
#   variables and coding.
#
# Principles
#   - Deposited variables are kept unchanged. Derived variables get new,
#     descriptive names ending in the age they refer to (e.g. ...14 = age
#     14).
#   - A derived variable is set to NA wherever its source isn't a
#     substantive answer (deposited negative codes such as refused, don't
#     know, not applicable, or absence from a sweep). The one exception is a
#     category that is itself meaningful (e.g. "no applicable resident-parent
#     record"). File 3 classifies every NA as item non-response, wave
#     non-response or structural (not applicable), using the response
#     indicators and explicit routing rules.
#   - Recodes applied through value labels (subject groups, and the
#     collapsed categories of five CLS non-response predictors) are written
#     to review tables, so they can be checked against the deposited labels.
#
# Contents
#   Block 3   survey design and response indicators
#   Block 4   age-14 background (analysis covariates)
#   Block 5   age-14 science and mathematics orientation
#   Block 6   age-15 science course route
#   Block 7   age-16 continuation, HE likelihood and intended HE subject
#   Block 8   subject-group domains (HE at 19-20, degree at 25)
#   Block 9   HE participation and subject at ages 19-20
#   Block 10  degree attained by age 25
#   Block 11  monthly activity, ages 16-19
#   Block 12  activity, employment and occupation at ages 25 and 32
#   Block 13  CLS predictors of non-response (Appendix C, Table C1)
#   Block 14  consistency checks
#   Block 15  review tables (metadata only)
#   Block 16  save
#   Block 17  summary
#
# Reads
#   data_dir/01_merged.rds
#
# Writes
#   data_dir/02_constructed.rds          deposited + derived variables
#                                        (participant level)
#   review_dir/02_derived_variables.csv  name, source and categories of each
#                                        derived variable (metadata only)
#   review_dir/02_value_labels.csv       deposited value labels relied on for
#                                        code-based recodes (metadata only)
#   review_dir/02_label_mappings.csv     deposited label -> derived category
#                                        for label-based recodes (metadata
#                                        only)
#   review_dir/02_w3_routing_check.csv   Wave-3 plan code x HE question codes,
#                                        banded counts (0 / fewer than 10 /
#                                        10 or more)
#
# References
#   Silverwood R, et al. (2026). Handling missing data in the CLS cohort
#   studies: User guide (Version 2), Appendix C. UCL Centre for Longitudinal
#   Studies.
#   Silverwood RJ, et al. (2024). Longitudinal and Life Course Studies 15(2):
#   227-250.
#   Centre for Longitudinal Studies. Stata code deriving the Next Steps
#   predictors of non-response: github.com/CLS-Data/variables-predicting-
#   non-response-in-the-cohorts (MIT licence).
###############################################################################

library(haven)
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
# Block 2. Read the merged cohort from File 1
###############################################################################

cohort <- readRDS(file.path(data_dir, "01_merged.rds"))

# Record the deposited value labels of every source variable whose numeric
# codes are used directly below. These are written to a review table in
# Block 15, so each assumed code can be checked against its label.
code_based_sources <- c(
  "W1OUTCOME", "W1sexYP", "W1ethgrpYP", "W1englangYP", "W1nssecfam",
  "W1hiqualgmum", "W1hiqualgdad", "W1inc1est", "W1famtyp2", "IndSchool",
  "W1lsciYP", "W1gsciYP", "W1favesYP", "W1lmathYP", "W1gmathYP",
  "W2GenSciYP", "W2ExSciYP", "W3plan16YP", "W3heposs9YP", "W3hesciYP",
  "W3hesubYP1", "W3HESUBYP21", "W6HEFlag", "W7HEFlag", "W8DDEGP", "W8DHANVQH",
  "W8DWRK", "W9DWRK", "W8DACTIVITYC",
  grep("^W7FinAct[0-9]+_B11$", names(cohort), value = TRUE)[1]
)


###############################################################################
# Block 3. Survey design and response indicators
#
# responded_w1 ... responded_w9: TRUE when the sweep's outcome code is 1
# (productive interview), FALSE otherwise. These drive the missingness
# classification in File 3 and are imputation auxiliaries.
#
# Weights are converted from SPSS format to plain numbers. A weight that is
# zero, negative or missing is set to NA (e.g. no Wave-9 final weight for
# Wave-9 non-respondents).
###############################################################################

for (wave in 1:9) {
  outcome_code <- as.numeric(cohort[[paste0("W", wave, "OUTCOME")]])
  cohort[[paste0("responded_w", wave)]] <- !is.na(outcome_code) & outcome_code == 1
}

cohort$design_weight <- as.numeric(zap_labels(zap_missing(cohort$DESIGNWEIGHT)))
cohort$design_weight[!(cohort$design_weight > 0)] <- NA

cohort$w9_final_weight <- as.numeric(zap_labels(zap_missing(cohort$W9FINWT)))
cohort$w9_final_weight[!(cohort$w9_final_weight > 0)] <- NA

cohort$psu     <- factor(as.numeric(cohort$SAMPPSU))
cohort$stratum <- factor(as.numeric(cohort$SAMPSTRATUM))


###############################################################################
# Block 4. Age-14 background (analysis covariates)
#
# Categorical variables keep their deposited value labels as categories.
# Negative codes become NA, except where noted.
###############################################################################

# Sex, ethnic group and main home language: deposited categories.
sex_code <- as.numeric(cohort$W1sexYP)
cohort$sex14 <- as_factor(cohort$W1sexYP, levels = "labels")
cohort$sex14[is.na(sex_code) | sex_code < 0] <- NA
cohort$sex14 <- droplevels(cohort$sex14)

ethnicity_code <- as.numeric(cohort$W1ethgrpYP)
cohort$ethnicity14 <- as_factor(cohort$W1ethgrpYP, levels = "labels")
cohort$ethnicity14[is.na(ethnicity_code) | ethnicity_code < 0] <- NA
cohort$ethnicity14 <- droplevels(cohort$ethnicity14)

language_code <- as.numeric(cohort$W1englangYP)
cohort$home_language14 <- as_factor(cohort$W1englangYP, levels = "labels")
cohort$home_language14[is.na(language_code) | language_code < 0] <- NA
cohort$home_language14 <- droplevels(cohort$home_language14)

# Family NS-SEC. Code -91 ("not applicable") is a genuine category (no
# classifiable family occupation), not missing, so it is kept as a level.
# Other negative codes become NA.
nssec_code   <- as.numeric(cohort$W1nssecfam)
nssec_labels <- attr(cohort$W1nssecfam, "labels")
if (any(!is.na(nssec_code) & nssec_code >= 0 & !(nssec_code %in% nssec_labels))) {
  stop("W1nssecfam has substantive codes without value labels.")
}
nssec_levels <- c(names(sort(nssec_labels[nssec_labels >= 0])), "Not applicable")

cohort$nssec_family14 <- as.character(as_factor(cohort$W1nssecfam, levels = "labels"))
cohort$nssec_family14[nssec_code %in% -91] <- "Not applicable"
cohort$nssec_family14[is.na(nssec_code) | (nssec_code < 0 & nssec_code != -91)] <- NA
cohort$nssec_family14 <- factor(cohort$nssec_family14, levels = nssec_levels)

# Highest resident-parent qualification: the higher (lower code) of the two
# parents' qualifications, codes 1 (degree) to 7 (no qualification).
# -98 and -91 mean there is no applicable parent record, which is structural.
# If one parent's record is unresolved (e.g. refused), the value stays NA,
# rather than assuming the observed parent has the higher qualification.
# If neither parent has an applicable record, the value is category 8, "No
# applicable resident-parent record".
mother_code <- as.numeric(cohort$W1hiqualgmum)
father_code <- as.numeric(cohort$W1hiqualgdad)
mother_observed   <- mother_code %in% 1:7
father_observed   <- father_code %in% 1:7
mother_structural <- mother_code %in% c(-98, -91)
father_structural <- father_code %in% c(-98, -91)
both_resolved <- (mother_observed | mother_structural) &
                 (father_observed | father_structural)

parent_education_code <- case_when(
  both_resolved & mother_observed & father_observed ~ pmin(mother_code, father_code),
  both_resolved & mother_observed                   ~ mother_code,
  both_resolved & father_observed                   ~ father_code,
  mother_structural & father_structural             ~ 8,
  TRUE                                              ~ NA_real_
)
cohort$highest_parent_education14 <- factor(
  parent_education_code,
  levels = 1:8,
  labels = c("Degree", "Higher education below degree", "A level",
             "GCSE grades A-C", "Level 1 or below", "Other qualification",
             "No qualification", "No applicable resident-parent record")
)

# Household income: banded income code 0-32, kept as a number (higher code =
# higher income band). Negative codes become NA.
income_code <- as.numeric(cohort$W1inc1est)
cohort$household_income14 <- ifelse(income_code %in% 0:32, income_code, NA)

# Single-parent household. The deposited family-type variable has two codes
# (0 and 2) that both mean "not a single-parent household"; 1 means single
# parent.
family_type_code <- as.numeric(cohort$W1famtyp2)
cohort$single_parent14 <- factor(
  case_when(
    family_type_code == 1          ~ "Single parent",
    family_type_code %in% c(0, 2)  ~ "Two parents",
    TRUE                           ~ NA_character_
  ),
  levels = c("Two parents", "Single parent")
)

# Independent (fee-paying) school at age 14: 1 = yes, 0 = no.
school_code <- as.numeric(cohort$IndSchool)
cohort$independent_school14 <- factor(
  case_when(
    school_code == 1 ~ "Independent school",
    school_code == 0 ~ "Maintained school",
    TRUE             ~ NA_character_
  ),
  levels = c("Maintained school", "Independent school")
)


###############################################################################
# Block 5. Age-14 science and mathematics orientation
#
# Liking and perceived competence are 4-point items where 1 is the most
# positive answer. They are reversed, so higher = more positive:
#   score = 5 - deposited code, giving 1 (least) to 4 (most).
# Science as favourite subject: deposited code 2 = science is a favourite
# subject (1), any other positive code = not (0).
# These are separate observed indicators, not a combined scale.
###############################################################################

science_liking_code     <- as.numeric(cohort$W1lsciYP)
science_competence_code <- as.numeric(cohort$W1gsciYP)
math_liking_code        <- as.numeric(cohort$W1lmathYP)
math_competence_code    <- as.numeric(cohort$W1gmathYP)
favourite_code          <- as.numeric(cohort$W1favesYP)

cohort$science_liking14 <- ifelse(science_liking_code %in% 1:4,
                                  5 - science_liking_code, NA)
cohort$science_competence14 <- ifelse(science_competence_code %in% 1:4,
                                      5 - science_competence_code, NA)
cohort$math_liking14 <- ifelse(math_liking_code %in% 1:4,
                               5 - math_liking_code, NA)
cohort$math_competence14 <- ifelse(math_competence_code %in% 1:4,
                                   5 - math_competence_code, NA)

cohort$science_favourite14 <- factor(
  case_when(
    favourite_code == 2 ~ "Science is a favourite subject",
    favourite_code >= 1 ~ "Science is not a favourite subject",
    TRUE                ~ NA_character_
  ),
  levels = c("Science is not a favourite subject",
             "Science is a favourite subject")
)


###############################################################################
# Block 6. Age-15 science course route (Wave 2)
#
# W2GenSciYP: 2 = separate sciences, 1 = general (double/single) science,
# 3 = no science course. W2ExSciYP (asked only when W2GenSciYP = 1):
# 1 = leads to a GCSE, 2 = doesn't.
#
# Three categories. A general science course that doesn't lead to a GCSE and
# no science course are both small, and share a meaning (no GCSE science),
# so they form one category (CLS user guide: collapse sparse categories).
###############################################################################

general_science_code <- as.numeric(cohort$W2GenSciYP)
general_gcse_code    <- as.numeric(cohort$W2ExSciYP)

cohort$science_course_route15 <- factor(
  case_when(
    general_science_code == 2                          ~ "Separate sciences",
    general_science_code == 1 & general_gcse_code == 1 ~ "General science GCSE",
    general_science_code == 1 & general_gcse_code == 2 ~ "No GCSE science course",
    general_science_code == 3                          ~ "No GCSE science course",
    TRUE                                               ~ NA_character_
  ),
  levels = c("Separate sciences", "General science GCSE", "No GCSE science course")
)


###############################################################################
# Block 7. Age-16 continuation, HE likelihood and intended HE subject (Wave 3)
#
# The Wave-3 questions:
#   W3plan16YP (asked of everyone interviewed): plan after Year 11.
#     1 = staying on in, or returning to, full-time education; 2 = leaving;
#     3 = leaving now but returning later (spontaneous); -1 = don't know;
#     -99 = young person not interviewed.
#   W3heposs9YP (asked of everyone interviewed): how likely to apply to
#     university.
#   W3hesciYP: whether the young person would like to study for a science
#     degree. Codes 0-1 (no/yes) mean a subject answer was recorded, so the
#     subject flags apply; -93 = not sure whether to go to university, or
#     subject not decided; -95 = not planning to go to university; -91 =
#     routed past the HE questions; -94 = insufficient information (item
#     non-response).
#   W3hesubYP1-20: the classified subject areas mentioned (0/1).
#   W3HESUBYP21: subject answer "don't know".
#
# How the questions were routed (checked in 02_w3_routing_check.csv):
#   - everyone leaving full-time education (plan 2 or 3) was routed past the
#     HE questions (W3hesciYP = -91), as were people whose plan was "don't
#     know";
#   - people staying on were asked, although some of them also have -91;
#   - in the deposited data, W3hesciYP is 0 ("no"), not -99, for young
#     people who weren't interviewed at Wave 3 (plan -99). Those values are
#     treated as missing.
#
# The intended-subject route follows this routing in three parts, which
# File 7 imputes and then recombines:
#   w3_continue_fte16     staying on in full-time education or leaving
#                         (the first part; also an analysis variable)
#   w3_route_stage16      for those staying on: subject area given / unsure
#                         or undecided / no university route. Not applicable
#                         for those leaving.
#   w3_intended_domain16  the subject domain, only where a subject area was
#                         given
# w3_intention_route16 (the analysis variable) is their combination:
#   leaving, or no university route  -> No university or full-time education
#                                       route
#   unsure or undecided              -> Unsure about university or subject
#                                       undecided
#   subject area given               -> Subject given: <domain>
# A "no university route" answer fixes the combination even when the plan
# itself is unknown. Conversely, anyone asked the HE questions was staying
# on, so a missing plan is completed as staying on for them (deductive
# completion, not imputation).
#
# People with a subject answer (W3hesciYP 0-1) but no classified subject:
#   - "don't know" only: they haven't decided their subject, which is what
#     code -93 records, so they are on the undecided route;
#   - refused, couldn't be classified, or no flag: they are on the
#     subject-given route, but their subject domain is missing (item
#     non-response), to be imputed.
###############################################################################

# Continuation in full-time education: 1 = continue or return; 2 or 3 =
# leave or not immediately return.
plan_code <- as.numeric(cohort$W3plan16YP)
cohort$w3_continue_fte16 <- factor(
  case_when(
    plan_code == 1         ~ "Continue full-time education",
    plan_code %in% c(2, 3) ~ "Leave or not immediately return",
    TRUE                   ~ NA_character_
  ),
  levels = c("Leave or not immediately return", "Continue full-time education")
)

# Likelihood of applying to university: 4 deposited categories, kept as an
# ordered factor in the deposited order.
he_likelihood_code <- as.numeric(cohort$W3heposs9YP)
cohort$he_likelihood16 <- as_factor(cohort$W3heposs9YP, levels = "labels")
cohort$he_likelihood16[!(he_likelihood_code %in% 1:4)] <- NA
cohort$he_likelihood16 <- factor(cohort$he_likelihood16, ordered = TRUE)
cohort$he_likelihood16 <- droplevels(cohort$he_likelihood16)

# Subject domains from the classified subject flags (1 = mentioned):
#   flags 1-2 = health (medicine, allied to medicine)
#   flags 3-9 = core STEM (sciences, mathematics, computing, engineering)
#   flag 10   = built environment
#   flags 11-20 = other subjects (including combined or general courses and
#                 other courses)
# Built environment is grouped with other subjects, because its cell is
# small.
subject_variables <- paste0("W3hesubYP", 1:20)
subject_flags <- matrix(NA_real_, nrow = nrow(cohort), ncol = 20,
                        dimnames = list(NULL, subject_variables))
for (k in 1:20) subject_flags[, k] <- as.numeric(cohort[[subject_variables[k]]])
chose_health    <- rowSums(subject_flags[, 1:2]   == 1, na.rm = TRUE) > 0
chose_core_stem <- rowSums(subject_flags[, 3:9]   == 1, na.rm = TRUE) > 0
chose_other     <- rowSums(subject_flags[, 10:20] == 1, na.rm = TRUE) > 0
domains_chosen  <- chose_health + chose_core_stem + chose_other
subject_dont_know <- as.numeric(cohort$W3HESUBYP21) %in% 1

# Route stage (second part), for those staying on in full-time education.
# Not applicable (NA) for those leaving. A subject answer of "don't know"
# alone (no classified subject) joins the undecided route (see the block
# header). W3hesciYP is set to missing for young people not interviewed at
# Wave 3 (deposited as 0).
subject_given_code <- as.numeric(cohort$W3hesciYP)
subject_given_code[plan_code %in% -99] <- NA

# Deductive completion of the plan: anyone asked the HE questions
# (W3hesciYP 0-1, -93 or -95) was staying on in full-time education, since
# everyone leaving was routed past them (see the block header). Where the
# plan itself is missing, it is completed accordingly. A check in Block 14
# stops the script if anyone leaving was asked these questions.
asked_he_questions <- subject_given_code %in% c(0, 1, -93, -95)
cohort$w3_continue_fte16[is.na(cohort$w3_continue_fte16) & asked_he_questions] <-
  "Continue full-time education"

dont_know_only <- subject_given_code %in% 0:1 & domains_chosen == 0 & subject_dont_know
leaving_fte    <- cohort$w3_continue_fte16 %in% "Leave or not immediately return"

cohort$w3_route_stage16 <- factor(
  case_when(
    leaving_fte                          ~ NA_character_,   # routed: not applicable
    dont_know_only                       ~ "Unsure about university or subject undecided",
    subject_given_code %in% 0:1          ~ "Subject area given",
    subject_given_code == -93            ~ "Unsure about university or subject undecided",
    subject_given_code %in% c(-91, -95)  ~ "No university route",
    TRUE                                 ~ NA_character_
  ),
  levels = c("Subject area given", "Unsure about university or subject undecided",
             "No university route")
)

# Intended subject domain, only on the subject-given route. With no
# classified subject (refused, couldn't be classified, or no flag), the
# domain is missing: item non-response.
cohort$w3_intended_domain16 <- factor(
  case_when(
    !(cohort$w3_route_stage16 %in% "Subject area given") ~ NA_character_,  # routed
    domains_chosen > 1             ~ "More than one domain",
    chose_core_stem                ~ "Core STEM",
    chose_health                   ~ "Health",
    chose_other                    ~ "Other or built environment",
    TRUE                           ~ NA_character_   # no classified subject
  ),
  levels = c("Core STEM", "Health", "Other or built environment",
             "More than one domain")
)

# The combination of the three parts, for observed data (recombined after
# imputation in File 7, with the same rules). Missing when a part it needs
# is missing.
cohort$w3_intention_route16 <- factor(
  case_when(
    leaving_fte | cohort$w3_route_stage16 %in% "No university route" ~
      "No university or full-time education route",
    cohort$w3_route_stage16 %in% "Unsure about university or subject undecided" ~
      "Unsure about university or subject undecided",
    cohort$w3_route_stage16 %in% "Subject area given" &
      !is.na(cohort$w3_intended_domain16) ~
      paste("Subject given:", cohort$w3_intended_domain16),
    TRUE ~ NA_character_
  ),
  levels = c(paste("Subject given:", levels(cohort$w3_intended_domain16)),
             "Unsure about university or subject undecided",
             "No university or full-time education route")
)


###############################################################################
# Block 8. Subject-group domains (HE at 19-20 and degree at 25)
#
# The Wave 6, 7 and 8 subject-group variables use different numeric codes,
# so each code is classified from its value-label wording, with rules
# applied in this order:
#   1. architecture, building or planning       -> Other or built environment
#   2. biology, veterinary, agriculture, physical sciences, mathematics,
#      computing, engineering or technology     -> Core STEM
#   3. medicine or allied to medicine           -> Health
#   4. any other substantive subject            -> Other or built environment
# Negative codes (not applicable, refused, etc.) -> NA.
# Every label and its assigned domain is written to 02_label_mappings.csv.
# A substantive code without a value label stops the script.
###############################################################################

domain_levels <- c("Core STEM", "Health", "Other or built environment")
subject_group_sources <- c(he_domain_w6 = "W6HESubGroup2",
                           he_domain_w7 = "W7HESubGroup",
                           degree_domain_w8 = "W8DHESUBGROUP")

label_mappings <- data.frame()

for (new_name in names(subject_group_sources)) {
  source_name  <- subject_group_sources[[new_name]]
  value_labels <- attr(cohort[[source_name]], "labels")
  label_text   <- tolower(names(value_labels))

  domain_for_label <- case_when(
    value_labels < 0 ~ NA_character_,
    grepl("architect|building|planning", label_text) ~ "Other or built environment",
    grepl("biolog|veterinary|agricultur|physical|mathemat|computer|engineering|technolog",
          label_text) ~ "Core STEM",
    grepl("medicine|allied", label_text) ~ "Health",
    TRUE ~ "Other or built environment"
  )

  label_mappings <- rbind(label_mappings, data.frame(
    source_variable  = source_name,
    code             = as.numeric(value_labels),
    deposited_label  = names(value_labels),
    derived_category = domain_for_label
  ))

  code <- as.numeric(cohort[[source_name]])
  unlabelled <- !is.na(code) & code >= 0 & !(code %in% value_labels)
  if (any(unlabelled)) {
    stop(source_name, " has substantive codes without value labels; ",
         "the subject-domain rules can't classify them.")
  }
  cohort[[new_name]] <- factor(domain_for_label[match(code, value_labels)],
                               levels = domain_levels)
}


###############################################################################
# Block 9. HE participation and subject at ages 19-20 (Waves 6 and 7)
#
# W6HEFlag / W7HEFlag: 1 = in higher education, 2 = not.
#   he_participation20: in HE at Wave 6 or 7 -> "In HE by age 20"; not in HE
#     at both waves -> "Not in HE by age 20"; otherwise NA (e.g. not in HE at
#     one wave and absent at the other).
#   he_domain20: the subject domain, only for those in HE. Core STEM takes
#     precedence, then health, then other, across the two waves.
#   he_pathway20: their combination, for observed data (recombined after
#     imputation in File 7).
###############################################################################

w6_he_code <- as.numeric(cohort$W6HEFlag)
w7_he_code <- as.numeric(cohort$W7HEFlag)

cohort$he_participation20 <- factor(
  case_when(
    w6_he_code == 1 | w7_he_code == 1 ~ "In HE by age 20",
    w6_he_code == 2 & w7_he_code == 2 ~ "Not in HE by age 20",
    TRUE                              ~ NA_character_
  ),
  levels = c("Not in HE by age 20", "In HE by age 20")
)

w6_domain <- as.character(cohort$he_domain_w6)
w7_domain <- as.character(cohort$he_domain_w7)

cohort$he_domain20 <- factor(
  case_when(
    !(cohort$he_participation20 %in% "In HE by age 20")              ~ NA_character_,
    w6_domain %in% "Core STEM" | w7_domain %in% "Core STEM"          ~ "Core STEM",
    w6_domain %in% "Health" | w7_domain %in% "Health"                ~ "Health",
    w6_domain %in% "Other or built environment" |
      w7_domain %in% "Other or built environment"                    ~ "Other or built environment",
    TRUE                                                             ~ NA_character_
  ),
  levels = domain_levels
)

cohort$he_pathway20 <- factor(
  case_when(
    cohort$he_participation20 == "Not in HE by age 20" ~ "Not in HE by age 20",
    !is.na(cohort$he_domain20) ~ paste("In HE:", cohort$he_domain20),
    TRUE ~ NA_character_
  ),
  levels = c("Not in HE by age 20", paste("In HE:", domain_levels))
)


###############################################################################
# Block 10. Degree attained by age 25 (Wave 8)
#
# W8DDEGP: 1 = has a degree, 0 = doesn't. Some respondents are routed to -1
# although their highest academic qualification (W8DHANVQH codes 1-3 or 96:
# below degree or none) shows they have no degree; they are coded "No
# degree".
#   degree_attained25: degree / no degree.
#   degree_domain25: the degree's subject domain, only for degree holders.
#   degree_pathway25: their combination, for observed data (recombined after
#     imputation in File 7).
###############################################################################

degree_code        <- as.numeric(cohort$W8DDEGP)
highest_qualification_code <- as.numeric(cohort$W8DHANVQH)
confirmed_no_degree <- degree_code == 0 |
  (degree_code == -1 & highest_qualification_code %in% c(1:3, 96))

cohort$degree_attained25 <- factor(
  case_when(
    degree_code == 1    ~ "Degree",
    confirmed_no_degree ~ "No degree",
    TRUE                ~ NA_character_
  ),
  levels = c("No degree", "Degree")
)

cohort$degree_domain25 <- cohort$degree_domain_w8
cohort$degree_domain25[cohort$degree_attained25 %in% "No degree"] <- NA
cohort$degree_domain25[is.na(cohort$degree_attained25)] <- NA

cohort$degree_pathway25 <- factor(
  case_when(
    cohort$degree_attained25 == "No degree" ~ "No degree by age 25",
    !is.na(cohort$degree_domain25) ~ paste("Degree:", cohort$degree_domain25),
    TRUE ~ NA_character_
  ),
  levels = c("No degree by age 25", paste("Degree:", domain_levels))
)


###############################################################################
# Block 11. Monthly activity, ages 16-19 (Waves 4-7 activity histories)
#
# One variable per month. Substantive activity codes:
#   1 = education, 4 = employment, 5 = apprenticeship or training,
#   6 = not in education, employment or training (NEET).
# Any other value (including negative codes) is treated as "month not
# observed".
#
#   months_activity_observed_16_19  months with a substantive code (also an
#                                   imputation auxiliary: it separates a
#                                   true zero from incomplete coverage)
#   months_education_16_19 / months_employed_16_19 /
#   months_apprenticeship_16_19 / months_neet_16_19
#                                   months in each activity
#   activity_transitions_16_19      changes of activity between adjacent
#                                   observed months
#   activity_months_expected        number of monthly variables (the
#                                   maximum possible)
###############################################################################

activity_variables <- grep("^W7FinAct[0-9]+_B11$", names(cohort), value = TRUE)
activity <- matrix(NA_real_, nrow = nrow(cohort), ncol = length(activity_variables),
                   dimnames = list(NULL, activity_variables))      # one column per month
for (k in seq_along(activity_variables)) activity[, k] <- as.numeric(cohort[[activity_variables[k]]])

activity_substantive <- activity == 1 | activity == 4 | activity == 5 | activity == 6

cohort$activity_months_expected       <- length(activity_variables)
cohort$months_activity_observed_16_19 <- rowSums(activity_substantive, na.rm = TRUE)
cohort$months_education_16_19      <- rowSums(activity == 1, na.rm = TRUE)
cohort$months_employed_16_19       <- rowSums(activity == 4, na.rm = TRUE)
cohort$months_apprenticeship_16_19 <- rowSums(activity == 5, na.rm = TRUE)
cohort$months_neet_16_19           <- rowSums(activity == 6, na.rm = TRUE)

# Transitions: compare each month with the next, counting only pairs where
# both months have a substantive code and the activity differs.
this_month <- activity[, -ncol(activity), drop = FALSE]
next_month <- activity[, -1, drop = FALSE]
both_substantive <- activity_substantive[, -ncol(activity), drop = FALSE] &
                    activity_substantive[, -1, drop = FALSE]
cohort$activity_transitions_16_19 <- rowSums(both_substantive & this_month != next_month,
                                             na.rm = TRUE)


###############################################################################
# Block 12. Activity, employment and occupation at ages 25 and 32
#
# W8DWRK / W9DWRK: 1 = in paid work, 2 = not in paid work.
# W8SOC3 / W9CURRJOB2010: current occupation, SOC2010 3-digit minor group.
# A valid code is a whole number from 100 to 999.
#
# Professional science, research, engineering and technology (SRET)
# employment = SOC2010 minor groups 211-215 (sub-major group 21).
#
# Age 25:
#   activity25            main economic activity (deposited categories)
#   employed25            in paid work / not
#   sret_professional_25  1 = SRET job; 0 = other job or not in work;
#                         NA = unresolved (incl. in work with no valid SOC)
#
# Age 32 (the outcome), in two routed parts plus their combination:
#   employed32              in paid work / not
#   sret_among_employed32   among those in work with a valid SOC: 1 = SRET,
#                           0 = other occupation. NA for those not in work
#                           (structural) and for unresolved cases.
#   sret_professional_32    overall outcome: 1 = SRET job; 0 = other job or
#                           not in work; NA = unresolved
###############################################################################

# Age 25
activity25_code <- as.numeric(cohort$W8DACTIVITYC)
cohort$activity25 <- as_factor(cohort$W8DACTIVITYC, levels = "labels")
cohort$activity25[is.na(activity25_code) | activity25_code < 0] <- NA
cohort$activity25 <- droplevels(cohort$activity25)

work25_code <- as.numeric(cohort$W8DWRK)
soc25       <- as.numeric(cohort$W8SOC3)
valid_soc25 <- !is.na(soc25) & soc25 >= 100 & soc25 <= 999 & soc25 == round(soc25)

cohort$employed25 <- factor(
  case_when(
    work25_code == 1 ~ "In paid work",
    work25_code == 2 ~ "Not in paid work",
    TRUE             ~ NA_character_
  ),
  levels = c("Not in paid work", "In paid work")
)

cohort$sret_professional_25 <- case_when(
  work25_code == 1 & valid_soc25 ~ as.integer(soc25 %in% 211:215),
  work25_code == 2               ~ 0L,
  TRUE                           ~ NA_integer_
)

# Age 32
work32_code <- as.numeric(cohort$W9DWRK)
soc32       <- as.numeric(cohort$W9CURRJOB2010)
valid_soc32 <- !is.na(soc32) & soc32 >= 100 & soc32 <= 999 & soc32 == round(soc32)

cohort$employed32 <- factor(
  case_when(
    work32_code == 1 ~ "In paid work",
    work32_code == 2 ~ "Not in paid work",
    TRUE             ~ NA_character_
  ),
  levels = c("Not in paid work", "In paid work")
)

cohort$sret_among_employed32 <- case_when(
  work32_code == 1 & valid_soc32 ~ as.integer(soc32 %in% 211:215),
  TRUE                           ~ NA_integer_
)

cohort$sret_professional_32 <- case_when(
  work32_code == 1 & valid_soc32 ~ as.integer(soc32 %in% 211:215),
  work32_code == 2               ~ 0L,
  TRUE                           ~ NA_integer_
)


###############################################################################
# Block 13. CLS predictors of non-response (Appendix C, Table C1)
#
# Candidate auxiliary variables for multiple imputation, screened in File 6.
# The codings follow the CLS user guide (Appendix C) and the CLS Stata code
# that derives these predictors ("Next Steps - final non-response variables
# dataset", github.com/CLS-Data/variables-predicting-non-response-in-the-
# cohorts):
#   - negative codes -> NA;
#   - 12 predictors keep their deposited categories (13a);
#   - parent's age is kept as a number (13b);
#   - 5 predictors have deposited codes collapsed exactly as in the CLS
#     code (13c);
#   - previous non-response is derived from the response outcomes (13d):
#     at Waves 2-7 as published (the sweeps before Wave 8), and with the
#     same definition over Waves 2-3 and Waves 2-8, for uses where the
#     sweeps known are those up to age 16 or up to age 25.
# W1sexYP is also a CLS predictor; it enters as sex14 from Block 4.
# W4DWPLnk2YP is also listed in Appendix C but isn't in the End User
# Licence release. The CLS code notes it came from an earlier Wave-4 file
# that is no longer in the accessible data.
###############################################################################

# 13a. Predictors kept in their deposited categories.
nonresponse_as_deposited <- c(
  nr_upset_by_names14           = "W1namesYP",
  nr_computer_games14           = "W1HcomGYP",
  nr_played_instrument14        = "W1alei1YP0h",
  nr_internet_at_home14         = "W1condur6MP",
  nr_school_contacted_parents15 = "W2SchcontMP",
  nr_under_strain15             = "W2strainYP",
  nr_smokes16                   = "W3cignowYP",
  nr_groups_treated_better18    = "W5FairTreat2YP",
  nr_teacher_expectations18     = "W5TeaExpecYP",
  nr_spoke_to_teacher19         = "W6IagIntroYP0b",
  nr_willing_sex_questions19    = "W6SexIntroYP",
  nr_willing_sex_questions20    = "W7SexIntroYP"
)

for (new_name in names(nonresponse_as_deposited)) {
  source_name <- nonresponse_as_deposited[[new_name]]
  code <- as.numeric(cohort[[source_name]])
  cohort[[new_name]] <- as_factor(cohort[[source_name]], levels = "labels")
  cohort[[new_name]][is.na(code) | code < 0] <- NA
  cohort[[new_name]] <- droplevels(cohort[[new_name]])
  code_based_sources <- c(code_based_sources, source_name)
}

# 13b. Parent's age at Wave 3, in years.
parent_age_code <- as.numeric(cohort$W3ageMP)
cohort$nr_parent_age16 <- ifelse(!is.na(parent_age_code) & parent_age_code > 0,
                                 parent_age_code, NA)

# 13c. Predictors whose deposited codes are collapsed, as in the CLS code:
#   W1gowhereYP  codes 3-5 -> "Sometimes-never"
#   W1hous12HH   codes 2-3 -> being bought on a mortgage/bank loan (code 3,
#                             shared ownership, is grouped with code 2)
#                codes 4-8 -> "Rented/other"
#   W4ClubberYP  codes 1-3 -> "Once a week or more"
#   W5SameaddHH  code 3    -> the same category as code 1 (yes)
#                code 4 ("not to be asked", household grid) -> NA: it is a
#                positive code that records the question wasn't asked, so
#                the variable keeps Table C1's Yes/No categories
#   W5actYP      codes 2-3 -> "Working or part working and part college"
#                codes 4-7 -> "Other"
# Every other substantive code keeps its deposited label. Each rule is
# applied to the variable's value labels first, giving a code -> category
# table (written to 02_label_mappings.csv), which is then applied to the
# data.

gowhere_labels <- attr(cohort$W1gowhereYP, "labels")
gowhere_category <- case_when(
  gowhere_labels < 0      ~ NA_character_,
  gowhere_labels %in% 3:5 ~ "Sometimes-never",
  TRUE                    ~ names(gowhere_labels)
)

tenure_labels <- attr(cohort$W1hous12HH, "labels")
tenure_category <- case_when(
  tenure_labels < 0      ~ NA_character_,
  tenure_labels %in% 2:3 ~ "Being bought on a mortgage/bank loan (incl. shared ownership)",
  tenure_labels %in% 4:8 ~ "Rented/other",
  TRUE                   ~ names(tenure_labels)
)

clubbing_labels <- attr(cohort$W4ClubberYP, "labels")
clubbing_category <- case_when(
  clubbing_labels < 0      ~ NA_character_,
  clubbing_labels %in% 1:3 ~ "Once a week or more",
  TRUE                     ~ names(clubbing_labels)
)

address_labels <- attr(cohort$W5SameaddHH, "labels")
address_category <- case_when(
  address_labels < 0  ~ NA_character_,
  address_labels == 4 ~ NA_character_,
  address_labels == 3 ~ names(address_labels)[address_labels == 1],
  TRUE                ~ names(address_labels)
)

activity18_labels <- attr(cohort$W5actYP, "labels")
activity18_category <- case_when(
  activity18_labels < 0      ~ NA_character_,
  activity18_labels %in% 2:3 ~ "Working or part working and part college",
  activity18_labels %in% 4:7 ~ "Other",
  TRUE                       ~ names(activity18_labels)
)

collapsed_predictors <- list(
  nr_parents_know_where14 = list(source = "W1gowhereYP", labels = gowhere_labels,
                                 category = gowhere_category),
  nr_housing_tenure14     = list(source = "W1hous12HH", labels = tenure_labels,
                                 category = tenure_category),
  nr_nightclubs17         = list(source = "W4ClubberYP", labels = clubbing_labels,
                                 category = clubbing_category),
  nr_same_address18       = list(source = "W5SameaddHH", labels = address_labels,
                                 category = address_category),
  nr_main_activity18      = list(source = "W5actYP", labels = activity18_labels,
                                 category = activity18_category)
)

for (new_name in names(collapsed_predictors)) {
  item <- collapsed_predictors[[new_name]]

  label_mappings <- rbind(label_mappings, data.frame(
    source_variable  = item$source,
    code             = as.numeric(item$labels),
    deposited_label  = names(item$labels),
    derived_category = item$category
  ))

  code <- as.numeric(cohort[[item$source]])
  if (any(!is.na(code) & code >= 0 & !(code %in% item$labels))) {
    stop(item$source, " has substantive codes without value labels.")
  }

  # Categories are ordered by their deposited codes.
  categories_by_code <- item$category[order(item$labels)]
  category_levels    <- unique(categories_by_code[!is.na(categories_by_code)])

  cohort[[new_name]] <- factor(item$category[match(code, item$labels)],
                               levels = category_levels)
}

# 13d. Previous non-response: one or more sweeps without a productive
# interview, as in the CLS code. Appendix C defines it over the sweeps
# before its Wave-8 target (Waves 2-7). The same definition is also built
# over other spans of sweeps:
#   Waves 2-3: sweeps known by age 16 (prediction landmark L16)
#   Waves 2-7: the published item; sweeps known by age 20 (landmark L20)
#   Waves 2-8: sweeps before the Wave-9 outcome (multiple imputation), and
#              sweeps known by age 25 (landmarks L25a and L25b)
responded_w2_to_w3 <- as.matrix(cohort[paste0("responded_w", 2:3)])
cohort$nr_previous_nonresponse_w2_w3 <- factor(
  ifelse(rowSums(!responded_w2_to_w3) > 0,
         "One or more instances of non-response",
         "Complete response up to Wave 3"),
  levels = c("Complete response up to Wave 3",
             "One or more instances of non-response")
)

responded_w2_to_w7 <- as.matrix(cohort[paste0("responded_w", 2:7)])
cohort$nr_previous_nonresponse_w2_w7 <- factor(
  ifelse(rowSums(!responded_w2_to_w7) > 0,
         "One or more instances of non-response",
         "Complete response up to Wave 7"),
  levels = c("Complete response up to Wave 7",
             "One or more instances of non-response")
)

responded_w2_to_w8 <- as.matrix(cohort[paste0("responded_w", 2:8)])
cohort$nr_previous_nonresponse_w2_w8 <- factor(
  ifelse(rowSums(!responded_w2_to_w8) > 0,
         "One or more instances of non-response",
         "Complete response up to Wave 8"),
  levels = c("Complete response up to Wave 8",
             "One or more instances of non-response")
)


###############################################################################
# Block 14. Consistency checks
#
# Each check stops the script if it fails. No data values are printed.
###############################################################################

# The outcome's routed parts must agree with the overall outcome.
if (any(cohort$sret_professional_32 == 1 & cohort$employed32 != "In paid work",
        na.rm = TRUE)) {
  stop("An SRET outcome of 1 at age 32 occurs without paid work.")
}
if (any(!is.na(cohort$sret_among_employed32) &
        !(cohort$employed32 %in% "In paid work"))) {
  stop("sret_among_employed32 is defined for someone not in paid work.")
}
if (any(cohort$sret_professional_32 != cohort$sret_among_employed32, na.rm = TRUE)) {
  stop("Overall and among-employed SRET outcomes disagree for workers.")
}

# Routed parts must be defined only where they apply.
if (any(cohort$w3_continue_fte16 %in% "Leave or not immediately return" &
        asked_he_questions)) {
  stop("Someone leaving full-time education was asked the Wave-3 HE ",
       "questions, so the routing assumed in Block 7 doesn't hold.")
}
if (any(!is.na(cohort$w3_route_stage16) &
        cohort$w3_continue_fte16 %in% "Leave or not immediately return")) {
  stop("A Wave-3 HE route exists for someone leaving full-time education.")
}
if (any(!is.na(cohort$w3_intended_domain16) &
        !(cohort$w3_route_stage16 %in% "Subject area given"))) {
  stop("An intended subject domain exists without a subject being given.")
}
if (any(!is.na(cohort$he_domain20) &
        !(cohort$he_participation20 %in% "In HE by age 20"))) {
  stop("An HE subject domain exists for someone not in HE.")
}
if (any(!is.na(cohort$degree_domain25) &
        !(cohort$degree_attained25 %in% "Degree"))) {
  stop("A degree subject domain exists for someone without a degree.")
}

# Activity months can't exceed the number of monthly variables.
if (any(cohort$months_activity_observed_16_19 > length(activity_variables))) {
  stop("More observed activity months than monthly variables.")
}

# Response indicators have no missing values.
if (anyNA(cohort[paste0("responded_w", 1:9)])) {
  stop("A response indicator is missing.")
}


###############################################################################
# Block 15. Review tables (metadata only)
###############################################################################

# 15a. Derived variables: name, source and categories.
derived_variables <- data.frame(
  variable = c(
    "responded_w1..w9", "design_weight", "w9_final_weight", "psu", "stratum",
    "sex14", "ethnicity14", "home_language14", "nssec_family14",
    "highest_parent_education14", "household_income14", "single_parent14",
    "independent_school14",
    "science_liking14", "science_competence14", "math_liking14",
    "math_competence14", "science_favourite14",
    "science_course_route15",
    "w3_continue_fte16", "he_likelihood16", "w3_route_stage16",
    "w3_intended_domain16", "w3_intention_route16",
    "he_participation20", "he_domain20", "he_pathway20",
    "degree_attained25", "degree_domain25", "degree_pathway25",
    "months_activity_observed_16_19", "months_education_16_19",
    "months_employed_16_19", "months_apprenticeship_16_19",
    "months_neet_16_19", "activity_transitions_16_19",
    "activity25", "employed25", "sret_professional_25",
    "employed32", "sret_among_employed32", "sret_professional_32",
    names(nonresponse_as_deposited), "nr_parent_age16",
    names(collapsed_predictors), "nr_previous_nonresponse_w2_w3",
    "nr_previous_nonresponse_w2_w7", "nr_previous_nonresponse_w2_w8"
  ),
  source = c(
    "W1OUTCOME..W9OUTCOME", "DESIGNWEIGHT", "W9FINWT", "SAMPPSU", "SAMPSTRATUM",
    "W1sexYP", "W1ethgrpYP", "W1englangYP", "W1nssecfam",
    "W1hiqualgmum + W1hiqualgdad", "W1inc1est", "W1famtyp2", "IndSchool",
    "W1lsciYP", "W1gsciYP", "W1lmathYP", "W1gmathYP", "W1favesYP",
    "W2GenSciYP + W2ExSciYP",
    "W3plan16YP", "W3heposs9YP",
    "W3plan16YP + W3hesciYP + W3hesubYP1-20 + W3HESUBYP21",
    "W3hesciYP + W3hesubYP1-20 + W3HESUBYP21",
    "w3_continue_fte16 + w3_route_stage16 + w3_intended_domain16",
    "W6HEFlag + W7HEFlag", "W6HESubGroup2 + W7HESubGroup",
    "he_participation20 + he_domain20",
    "W8DDEGP + W8DHANVQH", "W8DHESUBGROUP", "degree_attained25 + degree_domain25",
    rep("W7FinAct*_B11 (monthly)", 6),
    "W8DACTIVITYC", "W8DWRK", "W8DWRK + W8SOC3",
    "W9DWRK", "W9DWRK + W9CURRJOB2010", "W9DWRK + W9CURRJOB2010",
    unname(nonresponse_as_deposited), "W3ageMP",
    "W1gowhereYP", "W1hous12HH", "W4ClubberYP", "W5SameaddHH", "W5actYP",
    "W2OUTCOME..W3OUTCOME", "W2OUTCOME..W7OUTCOME", "W2OUTCOME..W8OUTCOME"
  )
)

derived_variables$categories <- NA_character_
for (i in seq_len(nrow(derived_variables))) {
  name <- derived_variables$variable[i]
  if (name %in% names(cohort) && is.factor(cohort[[name]])) {
    derived_variables$categories[i] <- paste(levels(cohort[[name]]), collapse = " | ")
  } else if (name %in% names(cohort)) {
    derived_variables$categories[i] <- paste("numeric:", class(cohort[[name]])[1])
  }
}
derived_variables$disclosure_checked <- TRUE
write.csv(derived_variables, file.path(review_dir, "02_derived_variables.csv"),
          row.names = FALSE, na = "")

# 15b. Deposited value labels of every source variable used through its codes.
value_labels_table <- data.frame()
for (source_name in unique(code_based_sources)) {
  value_labels <- attr(cohort[[source_name]], "labels")
  if (length(value_labels) > 0) {
    value_labels_table <- rbind(value_labels_table, data.frame(
      source_variable = source_name,
      code            = as.numeric(value_labels),
      deposited_label = names(value_labels)
    ))
  }
}
value_labels_table$disclosure_checked <- TRUE
write.csv(value_labels_table, file.path(review_dir, "02_value_labels.csv"),
          row.names = FALSE, na = "")

# 15c. Label-based recodes: deposited label -> derived category.
label_mappings$disclosure_checked <- TRUE
write.csv(label_mappings, file.path(review_dir, "02_label_mappings.csv"),
          row.names = FALSE, na = "")

# 15d. Wave-3 routing check: which codes of the HE-subject question
# (W3hesciYP) and the HE-likelihood question (W3heposs9YP) occur with each
# post-16 plan code (W3plan16YP). It shows how the HE questions were routed,
# e.g. whether everyone with plan code 2 or 3 was routed past them (code -91).
# --- Disclosure check: counts are replaced by bands (0 / fewer than 10 / 10
# or more), so no exact count is shown. ---------------------------------------
plan_code_all   <- as.numeric(cohort$W3plan16YP)
hesci_code_all  <- as.numeric(cohort$W3hesciYP)
heposs_code_all <- as.numeric(cohort$W3heposs9YP)

routing_check <- rbind(
  data.frame(question = "W3hesciYP",
             as.data.frame(table(plan_code = plan_code_all, answer_code = hesci_code_all,
                                 useNA = "ifany"))),
  data.frame(question = "W3heposs9YP",
             as.data.frame(table(plan_code = plan_code_all, answer_code = heposs_code_all,
                                 useNA = "ifany")))
)
routing_check$people <- ifelse(routing_check$Freq == 0, "0",
                               ifelse(routing_check$Freq < 10, "fewer than 10",
                                      "10 or more"))
routing_check$Freq <- NULL
routing_check$disclosure_checked <- TRUE
write.csv(routing_check, file.path(review_dir, "02_w3_routing_check.csv"),
          row.names = FALSE, na = "")
# -----------------------------------------------------------------------------


###############################################################################
# Block 16. Save
###############################################################################

saveRDS(cohort, file.path(data_dir, "02_constructed.rds"), compress = "xz")


###############################################################################
# Block 17. Summary
#
# Counts of variables and checks only. No data values.
###############################################################################

cat("\nFile 2 finished.\n")
cat("  Derived variables documented:   ", nrow(derived_variables), "\n")
cat("  Monthly activity variables used: ", length(activity_variables), "\n")
cat("  Consistency checks:               all passed\n")
cat("  Review tables written (metadata only):\n")
cat("    02_derived_variables.csv, 02_value_labels.csv, 02_label_mappings.csv\n")
cat("  Wave-3 routing check (banded counts): 02_w3_routing_check.csv\n")
cat("  Saved:", file.path(data_dir, "02_constructed.rds"), "\n")
