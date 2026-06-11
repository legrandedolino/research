################################################################################
#
# 1. DATA CLEANING - EDUCATING FOR ENVIRONMENTAL CHANGE
# Latest Version: Jun 10, 2026
# 
# This file processes both the pre- and post-test excel forms for further
# data analysis.
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
          "ccUnderstanding_8", "workshopTopics_11", # These onward are signposts
          "subjectsTeach_6_TEXT", "degreeOther") 

# Generate unique IDs from birth month, middle initial, and mother's initial
# Removes question text, and rows without IDs or consent
# Removes metadata columns
pre <- pre %>%
  slice(-1) %>% # Removes question text
  mutate(ID = toupper(paste0(substr(birthMonth, 1, 3), 
                             substr(middleInitial, 1,1), 
                             substr(motherInitial,1,1))),
         test = "pre",
         year = format(as.Date(as.numeric(EndDate), origin = "1899-12-30"), 
                       "%Y")) %>% 
  relocate(ID, test, year) %>%
  filter(Status != "Survey Preview" & ID != "NANANA" & 
           consent == "I agree") %>%
  select(-any_of(meta))

# Do the same for post, but also remove the first three entries
post <- post %>%
  slice(-(1:4)) %>% 
  mutate(ID = toupper(paste0(substr(birthMonth, 1, 3), 
                             substr(middleInitial, 1,1), 
                             substr(motherInitial,1,1))),
         test = "post",
         year = format(as.Date(as.numeric(EndDate), origin = "1899-12-30"), 
                       "%Y")) %>% 
  relocate(ID, test) %>%
  filter(Status != "Survey Preview" & ID != "NANANA") %>%
  select(-any_of(meta))
  
# Merges into one large data frame and sort by ID
data <- bind_rows(pre, post) %>%
  arrange(ID)

# Extracts a list of IDs that do not have exactly 2 rows
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
View(unpaired_summary) # Manual verification was done after this step.

# Repeat teacher IDs - No issues
reps <- c("NOVJP", "NOVKA", "SEPAN", "SEPDE", "SEPLC", "SEPMO")

# Repeat teacher IDs - With issues and for deletion
data <- data %>%
  group_by(ID, year, test) %>%
  mutate(occurrence = row_number()) %>%
  ungroup() %>%
  filter(
    !(ID == "SEPMA" & year == 2026),
    !(ID == "MARAD" & year == 2026),
    !(ID == "MARAD" & year == 2024 & test == "post" & occurrence > 1)
  ) %>%
  select(-occurrence)

# Conversions - Strong evidence of a typo
data <- data %>%
  mutate(ID = case_match(
    ID,
    "SEPMB" ~ "SEPMM",
    "SEPMD" ~ "SEPLD",
    "FEBTW" ~ "FEBTM",
    "AUGAJ" ~ "AUGAH",
    "AUGME" ~ "AUGMB",
    "FEBBW" ~ "FEBBB",
    "JUNRY" ~ "JUNRD",
    .default = ID  # Keeps all other IDs exactly as they are
  ))

# Orphaned - For deletion due to lack of pairing
orphs <- c("SEPTY",
           "JUNDS",
           "SEPBI",
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
           "JUNDF",
           "SEPAM"
           )

data <- data %>%
  filter(!ID %in% orphs)

# Verify that data is cleaned
unpaired_summary <- data %>%
  filter(ID %in% unpaired_ids) %>%
  select(ID, test, year) %>%
  arrange(ID)
View(unpaired_summary) # What shows are repeated but paired rows

################################################################################
# PRELIMINARY DESCRIPTIVES
################################################################################
yearly_counts <- data %>%
  filter(test == "pre") %>%
  group_by(year) %>%
  summarize(Total_Participants = n())

# Gender Aggregated by Year (Wide Format)
gender_by_year <- data %>%
  filter(test == "pre") %>%
  group_by(year, gender) %>%
  summarize(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = gender, values_from = Count, values_fill = 0)

# Race Aggregated by Year (Wide Format)
race_by_year <- data %>%
  filter(test == "pre") %>%
  group_by(year, race) %>%
  summarize(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = race, values_from = Count, values_fill = 0)

write.csv(data, "data.csv")
