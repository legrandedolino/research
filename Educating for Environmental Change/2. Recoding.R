################################################################################
#
# 2. RECODING - EDUCATING FOR ENVIRONMENTAL CHANGE
# Latest Version: Jun 11, 2026
# 
# This file processes items into workable values for analysis in R.
#
################################################################################

################################################################################
# LOADING IN DATA AND INITIAL DEFINITIONS
################################################################################
setwd("~/Documents/EfEC")
data <- read.csv("data.csv")

library(dplyr)
library(tidyr)
library(stringr)

################################################################################
# RECODING DATA
################################################################################

# Define scale for ccUnder
data_clean <- data %>%
  mutate(
    across(
      starts_with("ccUnder_"), 
      ~ case_match(
        .x,
        "Strongly disagree" ~ 1,
        "Disagree"          ~ 2,
        "Somewhat disagree" ~ 3,
        "Somewhat agree"    ~ 4,
        "Agree"             ~ 5,
        "Strongly agree"    ~ 6,
        .default = NA 
      )
    )
  )

# Define scale for Efficacy
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("Efficacy_"), 
      ~ case_match(
        .x,
        "Strongly disagree" ~ 1,
        "Disagree"          ~ 2,
        "Slightly disagree" ~ 3,
        "Slightly agree"    ~ 4,
        "Agree"             ~ 5,
        "Strongly agree"    ~ 6,
        .default = NA 
      )
    )
  )

# Define scale for kse
# NOTE: Since kse in 2022 is in a different scale, we cannot report them under
# the same construct. Instead, it will make sense to report them as background.
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("kse"),
      ~ case_match(
        .x,
        # The "Agree" Scale
        "Strongly disagree"         ~ 1,
        "Disagree"                  ~ 2,
        "Somewhat disagree"         ~ 3,
        "Somewhat agree"            ~ 4,
        "Agree"                     ~ 5,
        "Strongly agree"            ~ 6,
        
        # The "Importance" Scale
        "Not very important at all" ~ 1, 
        "Slightly important"        ~ 2, 
        "Moderately important"      ~ 3, 
        "Very important"            ~ 4, 
        "Extremely important"       ~ 5,
        
        # Anything else gets wiped
        .default = NA_real_ 
      )
    )
  )

# Define scale for schoolRate
# Not converting these to numbers since there's no pre-post comparison
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("schoolRate_"), 
      ~ factor(
        .x,
        levels = c("Poor", "Fair", "Average", "Good", "Excellent"),
        ordered = TRUE
      )
    )
  )

# Define scale for studentUnderstand
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("studentUnderstand_"), 
      ~ factor(
        .x,
        levels = c("Strongly disagree", "Disagree", "Somewhat disagree",
                   "Somewhat agree", "Agree", "Strongly agree"),
        ordered = TRUE
      )
    )
  )

# Define scale for topic
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("topic_"), 
      ~ case_match(
        .x,
        "I am unfamiliar with this topic" ~ 1, 
        "I know a little about this topic" ~ 2, 
        "I know a moderate amount about this topic" ~ 3, 
        "I know a lot about this topic" ~ 4,
        .default = NA 
      )
    )
  )

# Define scale for challenge (pre)
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("challenge_"), 
      ~ if_else(
        test == "pre", 
        as.character(
          case_match(
            .x,
            "Not at all challenging" ~ "1", 
            "Somewhat challenging" ~ "2", 
            "Moderately challenging" ~ "3", 
            "Challenging" ~ "4", 
            "Extremely challenging" ~ "5",
            .default = .x
          )
        ),
        .x
      )
    )
  )


# Define scale for challenge (post)
# NOTE: Same issue with kse (2022), cannot make pre-post comparisons
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("challenge_"), 
      ~ if_else(
        test == "post", 
        case_match(
          .x,
          "Not at all prepared" ~ "1", 
          "Somewhat prepared"   ~ "2", 
          "Moderately prepared" ~ "3", 
          "Completely prepared" ~ "4",
          .default = .x
        ),
        .x
      )
    )
  )

# Defining subjectsTeach (pre)
data_clean <- data_clean %>%
  mutate(
    temp_subj = replace_na(subjectsTeach, ""),
    subElem = as.integer(str_detect(subjectsTeach, fixed("Elementary"))),
    subMidLife = as.integer(str_detect(subjectsTeach, fixed("MS Life Science"))),
    subMidEar = as.integer(str_detect(subjectsTeach, fixed("MS Earth Science"))),
    subHiLife = as.integer(str_detect(subjectsTeach, fixed("HS Life Science"))),
    subHiEar = as.integer(str_detect(subjectsTeach, fixed("HS Earth Science"))),
    subOther = as.integer(str_detect(subjectsTeach, fixed("Other (please specify)"))),
    subUnc = as.integer(is.na(subjectsTeach) | subjectsTeach == "NA" | subjectsTeach == "")
  ) %>%
  mutate(
    across(subElem:subUnc, ~ if_else(test == "pre", .x, NA_integer_))
  ) %>%
  select(-RecordedDate, -temp_subj, -subjectsTeach, -subjectsTeach_6_TEXT)

# Defining classType (pre)
data_clean <- data_clean %>%
  mutate(
    indGen = as.integer(str_detect(classType, fixed("General / Regular"))),
    indInc = as.integer(str_detect(classType, fixed("Inclusion"))),
    indSpE = as.integer(str_detect(classType, fixed("Special Education (e.g., resource room, self-contained)"))),
    indESL = as.integer(str_detect(classType, fixed("ESL (English as the Second Language)"))),
    indAdv = as.integer(str_detect(classType, fixed("Advanced / Gifted and Talented")))
  ) %>%
  mutate(across(starts_with("ind_"), ~ tidyr::replace_na(., 0))) %>%
  select(-classType)

# Define scale for yearsTeaching
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("yearsTeaching_"), 
      ~ factor(
        .x,
        levels = c("1 - 5 years", "6 - 10 years", "11 - 15 years",
                   "16 - 20 years", "More than 20 years"),
        ordered = TRUE
      )
    )
  )

# Define scale for education
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("education_"), 
      ~ factor(
        .x,
        levels = c("Bachelor's degree", "Graduate work but no advanced degree",
                   "Master's degree", "Post-master's work but no doctorate",
                   "Ed.D. or PhD."),
        ordered = TRUE
      )
    )
  )

# Define scale for degrees
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("degrees_"), 
      ~ factor(
        .x,
        levels = c("No Degree", "Bachelor's Degree",
                   "Master's Degree", "Ed.D. or Ph.D."),
        ordered = TRUE
      )
    )
  )

# Define scale for studentSocioStatus
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("studentSocioStatus_"), 
      ~ factor(
        .x,
        levels = c("Lower income class", "Lower-middle income class",
                   "Middle income class", "Middle-upper income class",
                   "Upper income class"),
        ordered = TRUE
      )
    )
  )

# Define scale for ccControversial
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("ccControversial"), 
      ~ factor(
        .x,
        levels = c("Not at all controversial", "Not very controversial",
                   "Somewhat controversial", "Very controversial"),
        ordered = TRUE
      )
    )
  )

# gender, race, previousParticipant do not need to be ordered

################################################################################
# SETTING LEVELS FOR EVALUATIONS
################################################################################
# Again, no need to convert this to numbers as sophisticated quantitative
# analysis is not appropriate.

# Define scale for workshopObj
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("workshopObj_"), 
      ~ factor(
        .x,
        levels = c("Strongly disagree", "Disagree", "Somewhat disagree",
                   "Somewhat agree", "Agree", "Strongly agree"),
        ordered = TRUE
      )
    )
  )

# Define scale for workshopRating
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("workshopRating_"), 
      ~ factor(
        .x,
        levels = c("Poor", "Fair", "Average",
                   "Very good", "Excellent"),
        ordered = TRUE
      )
    )
  )

# Define scale for workshopAct
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("workshopAct_"), 
      ~ factor(
        .x,
        levels = c("Not at all helpful", "Somewhat helpful",
                   "Moderately helpful", "Very helpful"),
        ordered = TRUE
      )
    )
  )

################################################################################
# SAVING DATA
################################################################################
save(data_clean, file="data_clean.RData")
