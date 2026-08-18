################################################################################
#
# 1. DATA CLEANING - EDUCATING FOR ENVIRONMENTAL CHANGE
# Latest Version: August 18, 2026
#
# Processes the pre- and post-test Excel forms for analysis.
#
################################################################################

################################################################################
# LOADING IN DATA AND INITIAL CLEANING
################################################################################

setwd("~/Documents/EfEC") # Change this for your own working directory

library(readxl)
library(dplyr)
library(tidyr)

pre <- read_xlsx("pre.xlsx")
post <- read_xlsx("post.xlsx")

# Identify metadata and signpost columns for deletion
meta <- c("StartDate", "EndDate", "Status", "IPAddress", "Progress",
          "Duration (in seconds)", "Finished", "ResponseId", "consent",
          "RecipientLastName", "RecipientFirstName", "RecipientEmail",
          "ExternalReference", "LocationLatitude", "LocationLongitude",
          "DistributionChannel", "UserLanguage",
          "birthMonth", "middleInitial", "pets", "motherInitial",
          "ccUnderstanding_8", "workshopTopics_11", # Signposts start here
          "subjectsTeach_6_TEXT", "degreeOther")

# Generate IDs and remove question text, invalid responses, and metadata
pre <- pre %>%
  slice(-1) %>%
  mutate(
    ID = toupper(paste0(
      substr(birthMonth, 1, 3),
      substr(middleInitial, 1, 1),
      substr(motherInitial, 1, 1)
    )),
    test = "pre",
    year = format(
      as.Date(as.numeric(EndDate), origin = "1899-12-30"),
      "%Y"
    )
  ) %>%
  relocate(ID, test, year) %>%
  filter(
    Status != "Survey Preview" &
      ID != "NANANA" &
      consent == "I agree"
  ) %>%
  select(-any_of(meta))

# Process post-test data and correct mislabeled columns
post <- post %>%
  slice(-(1:4)) %>%
  mutate(
    ID = toupper(paste0(
      substr(birthMonth, 1, 3),
      substr(middleInitial, 1, 1),
      substr(motherInitial, 1, 1)
    )),
    test = "post",
    year = format(
      as.Date(as.numeric(EndDate), origin = "1899-12-30"),
      "%Y"
    )
  ) %>%
  rename(
    topic_climateCauses = topic_climateCause,
    kse_understand = kse_understanding
  ) %>%
  relocate(ID, test) %>%
  filter(Status != "Survey Preview" & ID != "NANANA") %>%
  select(-any_of(meta))

# Merge datasets and sort by ID
data <- bind_rows(pre, post) %>%
  mutate(
    race = recode(
      race,
      "White,American Indian or Alaska Native" =
        "American Indian or Alaska Native"
    )
  ) %>%
  arrange(ID)

# Identify IDs without one pretest and one posttest in the same year
unpaired_ids <- data %>%
  group_by(ID) %>%
  filter(
    n() != 2 |
      sum(test == "pre") != 1 |
      sum(test == "post") != 1 |
      n_distinct(year) != 1
  ) %>%
  pull(ID) %>%
  unique()
unpaired_summary <- data %>%
  filter(ID %in% unpaired_ids) %>%
  select(ID, test, year) %>%
  arrange(ID)
View(unpaired_summary) # Manually verify the flagged records

# Valid repeat-participant IDs
reps <- c("NOVJP", "NOVKA", "SEPAN", "SEPDE", "SEPLC")

# Distinguish the unrelated 2025 and 2026 SEPMO respondents
data <- data %>%
  mutate(ID = if_else(ID == "SEPMO" & year == 2026, "SEPMX", ID))

# Remove invalid repeat-participant records
data <- data %>%
  group_by(ID, year, test) %>%
  mutate(occurrence = row_number()) %>%
  ungroup() %>%
  filter(
    !(ID == "SEPMA" & year == 2026),
    !(ID == "MARAD" & year == 2026),
    !(ID == "MARAD" & year == 2024 & test == "post" & occurrence > 1),
    !(ID == "JUNDS" & year == 2023)
  ) %>%
  select(-occurrence)

# Correct IDs with strong evidence of typographical errors
data <- data %>%
  mutate(ID = replace_values(
    ID,
    "SEPMB" ~ "SEPMM",
    "SEPMD" ~ "SEPLD",
    "FEBTW" ~ "FEBTM",
    "AUGAJ" ~ "AUGAH",
    "AUGME" ~ "AUGMB",
    "FEBBW" ~ "FEBBB",
    "JUNRY" ~ "JUNRD",
    "JUNDF" ~ "JUNDS"
  ))

# Remove records without a matching pretest or posttest
orphs <- c("SEPTY",
           "NOVHL",
           "MAYHC",
           "MAYKT",
           "MAYML",
           "APRIE",
           "APRRM",
           "DECRC",
           "DECWK",
           "DECZC",
           "MARAR",
           "MARGJ",
           "FEBKK",
           "FEBJM",
           "JANAE",
           "JANAM",
           "JULLJ",
           "SEPAM")

data <- data %>%
  filter(!ID %in% orphs)

# Verify the remaining flagged records
unpaired_summary <- data %>%
  filter(ID %in% unpaired_ids) %>%
  select(ID, test, year) %>%
  arrange(ID)
View(unpaired_summary) # Remaining records are valid repeated pairs

################################################################################
# PRELIMINARY DESCRIPTIVES
################################################################################

yearly_counts <- data %>%
  filter(test == "pre") %>%
  group_by(year) %>%
  summarize(Total_Participants = n())

# Summarize gender by year in wide format
gender_by_year <- data %>%
  filter(test == "pre") %>%
  group_by(year, gender) %>%
  summarize(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = gender, values_from = Count, values_fill = 0)

# Summarize race by year in wide format
race_by_year <- data %>%
  filter(test == "pre") %>%
  group_by(year, race) %>%
  summarize(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = race, values_from = Count, values_fill = 0)

write.csv(data, "data.csv")
