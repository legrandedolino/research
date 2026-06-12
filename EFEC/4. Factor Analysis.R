################################################################################
#
# 4. FACTOR ANALYSIS - EDUCATING FOR ENVIRONMENTAL CHANGE
# Latest Version: Jun 11, 2026
# 
# Generates exploratory factor analysis and preliminary pre/post results.
#
################################################################################

################################################################################
# LOADING IN DATA
################################################################################

library(psych)
library(dplyr)
library(readxl)
library(tidyr)
library(broom)
library(effectsize)

load("data_clean.RData")

################################################################################
# EFFICACY
################################################################################

pre_eff_num <- data_clean %>% 
  filter(test == "pre") %>%
  select(starts_with("Efficacy_")) %>%
  # REVERSE CODE the negatively worded items
  mutate(across(c(Efficacy_2, Efficacy_4, Efficacy_5, Efficacy_7, 
                  Efficacy_8, Efficacy_9, Efficacy_11), ~ 7 - .)) %>%
  na.omit()

## Check for assumptions
# Kaiser-Meyer-Olkin (KMO) Measure of Sampling Adequacy
# Values > 0.60 are acceptable, > 0.80 is great.
KMO(pre_eff_num)

# Bartlett's Test of Sphericity
# We want the p-value to be < 0.05 (meaning items are correlated enough to factor)
cortest.bartlett(cor(pre_eff_num), n = nrow(pre_eff_num))

## Extract factors and model
# Run parallel analysis to determine the optimal number of factors
fa.parallel(pre_eff_num, fm = "pa", fa = "fa", main = "Scree Plot")

# Run the Exploratory Factor Analysis (EFA) model
# NOTE: Replace 'nfactors = 1' came from parallel analysis
efa_model <- fa(pre_eff_num, nfactors = 1, rotate = "oblimin", fm = "pa") 
print(efa_model, cut = 0.3, sort = TRUE)

# We drop Efficacy 10 due to low loading
final_eff_data <- pre_eff_num %>%
  select(-Efficacy_10)

# Re-run the final Factor Analysis to confirm the structure is stable
final_efa <- fa(final_eff_data, nfactors = 1, rotate = "oblimin", fm = "pa")
print(final_efa, cut = 0.3, sort = TRUE)

# Check Cronbach's Alpha
# Above 0,80 is great
efficacy_alpha <- alpha(final_eff_data)
print(summary(efficacy_alpha))

## Calculate composite score to both pre and post items
data_scored <- data_clean %>%
  mutate(across(c(Efficacy_2, Efficacy_4, Efficacy_5, Efficacy_7, 
                  Efficacy_8, Efficacy_9, Efficacy_11), ~ 7 - .)) %>%
  mutate(
    Efficacy_Score = rowMeans(
      select(., starts_with("Efficacy_"), -Efficacy_10), # Exclude item 10
      na.rm = TRUE
    )
  )

# Match and run pre and post items
paired_data <- data_scored %>%
  select(ID, year, test, Efficacy_Score) %>% 
  pivot_wider(
    names_from = test, 
    values_from = Efficacy_Score
  ) %>%
  drop_na(pre, post)

# Run t-test
paired_data %>%
  group_by(year) %>%
  summarize(
    tidy(t.test(post, pre, paired = TRUE)), 
    N_Pairs = n(),
    .groups = "drop"
  ) %>%
  select(
    Year = year, 
    N_Pairs,
    `Mean Difference (Post - Pre)` = estimate, 
    `p-value` = p.value, 
    `t-statistic` = statistic,
    `Degrees of Freedom` = parameter
  )

# Check effect size
# Above 0.80 is large
paired_data %>%
  group_by(year) %>%
  summarize(
    d_stats = list(cohens_d(post, pre, paired = TRUE))
  ) %>%
  tidyr::unnest(d_stats)

################################################################################
# UNDERSTANDING
################################################################################
pre_ccUnder_num <- data_clean %>% 
  filter(test == "pre") %>%
  select(starts_with("ccUnder_")) %>%
  na.omit()

# Check for assumptions
KMO(pre_ccUnder_num)
cortest.bartlett(cor(pre_ccUnder_num), n = nrow(pre_ccUnder_num))

# Extract factors and model
fa.parallel(pre_ccUnder_num, fm = "pa", fa = "fa", main = "Scree Plot")
efa_model_cc <- fa(pre_ccUnder_num, nfactors = 1, rotate = "oblimin", fm = "pa") 
print(efa_model_cc, cut = 0.3, sort = TRUE)
final_ccUnder_data <- pre_ccUnder_num

# Check reliability
ccUnder_alpha <- alpha(final_ccUnder_data)

# Calculate scores and match items
data_scored_cc <- data_clean %>%
  mutate(
    ccUnder_Score = rowMeans(
      select(., starts_with("ccUnder_")), 
      na.rm = TRUE
    )
  )

paired_data_cc <- data_scored_cc %>%
  select(ID, year, test, ccUnder_Score) %>% 
  pivot_wider(
    names_from = test, 
    values_from = ccUnder_Score
  ) %>%
  drop_na(pre, post)

# Run the paired t-test
paired_data_cc %>%
  group_by(year) %>%
  summarize(
    tidy(t.test(post, pre, paired = TRUE)), 
    N_Pairs = n(),
    .groups = "drop"
  ) %>%
  select(
    Year = year, 
    N_Pairs,
    `Mean Difference (Post - Pre)` = estimate, 
    `p-value` = p.value, 
    `t-statistic` = statistic,
    `Degrees of Freedom` = parameter
  )

paired_data_cc %>%
  group_by(year) %>%
  summarize(
    d_stats = list(cohens_d(post, pre, paired = TRUE))
  ) %>%
  tidyr::unnest(d_stats)
