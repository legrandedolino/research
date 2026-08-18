################################################################################
#
# 2. RECODING - EDUCATING FOR ENVIRONMENTAL CHANGE
# Latest Version: August 18, 2026
#
# Recodes survey items into values and factors for analysis.
#
################################################################################

################################################################################
# LOADING IN DATA AND INITIAL DEFINITIONS
################################################################################

setwd("~/Documents/EfEC")

library(dplyr)
library(tidyr)
library(stringr)

data <- read.csv("data.csv")

# Identify repeat participants
repeated <- data %>%
  add_count(ID) %>%
  filter(n == 4) %>%
  select(year, test, ID)

################################################################################
# RECODING DATA
################################################################################

# Recode climate-change understanding
data_clean <- data %>%
  mutate(
    across(
      starts_with("ccUnder_"),
      ~ case_match(
        .x,
        "Strongly disagree" ~ 1,
        "Disagree" ~ 2,
        "Somewhat disagree" ~ 3,
        "Somewhat agree" ~ 4,
        "Agree" ~ 5,
        "Strongly agree" ~ 6,
        .default = NA
      )
    )
  )

# Recode teaching efficacy
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("Efficacy_"),
      ~ case_match(
        .x,
        "Strongly disagree" ~ 1,
        "Disagree" ~ 2,
        "Slightly disagree" ~ 3,
        "Slightly agree" ~ 4,
        "Agree" ~ 5,
        "Strongly agree" ~ 6,
        .default = NA
      )
    )
  )

# Recode KSE; 2022 used a different scale and can only be reported as background
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("kse"),
      ~ case_match(
        .x,
        # Agreement scale
        "Strongly disagree" ~ 1,
        "Disagree" ~ 2,
        "Somewhat disagree" ~ 3,
        "Somewhat agree" ~ 4,
        "Agree" ~ 5,
        "Strongly agree" ~ 6,

        # Importance scale
        "Not very important at all" ~ 1,
        "Slightly important" ~ 2,
        "Moderately important" ~ 3,
        "Very important" ~ 4,
        "Extremely important" ~ 5,

        .default = NA_real_
      )
    )
  )

# Order school ratings without converting them to numeric values
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

# Order student-understanding responses
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

# Recode topic familiarity
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

# Recode pretest challenges
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

# Recode posttest preparedness; the 2022 scale is not comparable
data_clean <- data_clean %>%
  mutate(
    across(
      starts_with("challenge_"),
      ~ if_else(
        test == "post",
        case_match(
          .x,
          "Not at all prepared" ~ "1",
          "Somewhat prepared" ~ "2",
          "Moderately prepared" ~ "3",
          "Completely prepared" ~ "4",
          .default = .x
        ),
        .x
      )
    )
  )

# Create pretest teaching-subject indicators
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
  select(-RecordedDate, -temp_subj, -subjectsTeach)

# Create pretest class-type indicators
data_clean <- data_clean %>%
  mutate(
    indGen = as.integer(str_detect(classType, fixed("General / Regular"))),
    indInc = as.integer(str_detect(classType, fixed("Inclusion"))),
    indSpE = as.integer(str_detect(
      classType,
      fixed("Special Education (e.g., resource room, self-contained)")
    )),
    indESL = as.integer(str_detect(classType, fixed("ESL (English as the Second Language)"))),
    indAdv = as.integer(str_detect(classType, fixed("Advanced / Gifted and Talented")))
  ) %>%
  mutate(across(starts_with("ind_"), ~ tidyr::replace_na(., 0))) %>%
  select(-classType)

# Order years-teaching categories
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

# Order education categories
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

# Order degree categories
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

# Order student socioeconomic-status categories
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

# Order perceived-controversy categories
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

################################################################################
# SETTING LEVELS FOR EVALUATIONS
################################################################################

# Order workshop-objective responses
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

# Order workshop ratings
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

# Order workshop-activity responses
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

save(data_clean, file = "data_clean.RData")
