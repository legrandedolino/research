################################################################################
#
# 5. LINEAR MIXED-EFFECTS MODELS
# Latest Version: August 18, 2026
#
# Estimates cohort-adjusted pre/post change with participant random intercepts.
#
################################################################################

################################################################################
# PACKAGES AND SHARED SETTINGS
################################################################################

library(lmerTest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

options(contrasts = c("contr.treatment", "contr.poly"))

outcome_labels <- c(
  Efficacy_Score = "Teaching efficacy",
  ccUnder_Score = "Self-reported climate understanding"
)

theme_efec <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 11),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

################################################################################
# HELPERS
################################################################################

fit_outcome_model <- function(outcome, data) {
  model_formula <- as.formula(
    paste(outcome, "~ test + year + (1 | ID)")
  )

  lmer(
    model_formula,
    data = data,
    REML = TRUE,
    na.action = na.exclude,
    control = lmerControl(optimizer = "bobyqa")
  )
}

tidy_fixed_effects <- function(model, outcome) {
  coefficient_table <- as.data.frame(coef(summary(model)))
  coefficient_table$Term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL

  coefficient_table %>%
    transmute(
      Outcome = outcome,
      Term,
      Estimate = Estimate,
      SE = `Std. Error`,
      df = df,
      t = `t value`,
      p = `Pr(>|t|)`,
      CI_low = Estimate - qt(0.975, df) * SE,
      CI_high = Estimate + qt(0.975, df) * SE
    )
}

get_omnibus_tests <- function(model, outcome) {
  test_table <- as.data.frame(anova(model, ddf = "Satterthwaite"))

  data.frame(
    Outcome = outcome,
    Term = rownames(test_table),
    Num_df = test_table$NumDF,
    Den_df = test_table$DenDF,
    F = test_table$`F value`,
    p = test_table$`Pr(>F)`,
    row.names = NULL,
    check.names = FALSE
  )
}

# Calculate an equal-cohort-weighted adjusted mean for one test occasion
get_adjusted_mean <- function(model, outcome_code, test_level, data) {
  reference_grid <- expand.grid(
    test = levels(data$test),
    year = levels(data$year)
  )

  fixed_formula <- as.formula(paste(outcome_code, "~ test + year"))
  design_matrix <- model.matrix(
    delete.response(terms(fixed_formula)),
    reference_grid
  )
  design_matrix <- design_matrix[, names(fixef(model)), drop = FALSE]

  contrast_vector <- colMeans(
    design_matrix[reference_grid$test == test_level, , drop = FALSE]
  )

  result <- as.data.frame(
    contest1D(
      model,
      contrast_vector,
      ddf = "Satterthwaite",
      confint = TRUE,
      level = 0.95
    )
  )

  data.frame(
    Test = test_level,
    Adjusted_mean = result$Estimate,
    SE = result$`Std. Error`,
    df = result$df,
    CI_low = result$lower,
    CI_high = result$upper
  )
}

get_primary_result <- function(model, outcome_code, outcome, data) {
  pre_result <- get_adjusted_mean(model, outcome_code, "pre", data)
  post_result <- get_adjusted_mean(model, outcome_code, "post", data)

  fixed_formula <- as.formula(paste(outcome_code, "~ test + year"))
  reference_grid <- expand.grid(
    test = levels(data$test),
    year = levels(data$year)
  )
  design_matrix <- model.matrix(
    delete.response(terms(fixed_formula)),
    reference_grid
  )
  design_matrix <- design_matrix[, names(fixef(model)), drop = FALSE]

  pre_vector <- colMeans(
    design_matrix[reference_grid$test == "pre", , drop = FALSE]
  )
  post_vector <- colMeans(
    design_matrix[reference_grid$test == "post", , drop = FALSE]
  )

  change_result <- as.data.frame(
    contest1D(
      model,
      post_vector - pre_vector,
      ddf = "Satterthwaite",
      confint = TRUE,
      level = 0.95
    )
  )

  data.frame(
    Outcome = outcome,
    Observations = nobs(model),
    Unique_teachers = n_distinct(model.frame(model)$ID),
    Program_participations = n_distinct(
      interaction(model.frame(model)$ID, model.frame(model)$year, drop = TRUE)
    ),
    Adjusted_pre = pre_result$Adjusted_mean,
    Adjusted_pre_CI_low = pre_result$CI_low,
    Adjusted_pre_CI_high = pre_result$CI_high,
    Adjusted_post = post_result$Adjusted_mean,
    Adjusted_post_CI_low = post_result$CI_low,
    Adjusted_post_CI_high = post_result$CI_high,
    Adjusted_change = change_result$Estimate,
    SE = change_result$`Std. Error`,
    df = change_result$df,
    t = change_result$`t value`,
    p = change_result$`Pr(>|t|)`,
    CI_low = change_result$lower,
    CI_high = change_result$upper,
    Cohort_weighting = "Equal",
    row.names = NULL
  )
}

get_model_diagnostics <- function(model, outcome) {
  variance_components <- as.data.frame(VarCorr(model))
  participant_variance <- variance_components %>%
    filter(grp == "ID") %>%
    pull(vcov)
  residual_variance <- sigma(model)^2

  model_r2 <- suppressWarnings(performance::r2_nakagawa(model))
  standardized_residuals <- residuals(model) / sigma(model)

  data.frame(
    Outcome = outcome,
    Participant_intercept_SD = sqrt(participant_variance),
    Residual_SD = sigma(model),
    ICC = participant_variance / (participant_variance + residual_variance),
    Marginal_R2 = unname(model_r2$R2_marginal),
    Conditional_R2 = unname(model_r2$R2_conditional),
    Singular_fit = performance::check_singularity(model),
    Maximum_absolute_standardized_residual = max(
      abs(standardized_residuals),
      na.rm = TRUE
    ),
    Standardized_residuals_over_3 = sum(
      abs(standardized_residuals) > 3,
      na.rm = TRUE
    )
  )
}

################################################################################
# LOAD AND VALIDATE DATA
################################################################################

load("data_scored.RData")

required_columns <- c(
  "ID", "year", "test", "Efficacy_Score", "ccUnder_Score"
)

missing_columns <- setdiff(required_columns, names(data_scored))
if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

analysis_data <- data_scored %>%
  transmute(
    ID = factor(ID),
    year = factor(year, levels = sort(unique(as.character(year)))),
    test = factor(test, levels = c("pre", "post")),
    Efficacy_Score = as.numeric(Efficacy_Score),
    ccUnder_Score = as.numeric(ccUnder_Score)
  ) %>%
  arrange(year, ID, test)

duplicate_rows <- analysis_data %>%
  count(ID, year, test) %>%
  filter(n > 1)

if (nrow(duplicate_rows) > 0) {
  stop("Duplicate ID-year-test observations must be resolved before modeling.")
}

participation_check <- analysis_data %>%
  count(ID, year, name = "observations")

if (any(participation_check$observations != 2)) {
  stop("Every ID-year participation must have exactly one pretest and one posttest.")
}

sample_summary <- data.frame(
  Observations = nrow(analysis_data),
  Unique_teachers = n_distinct(analysis_data$ID),
  Program_participations = nrow(participation_check),
  Returning_teachers = sum(table(participation_check$ID) > 1)
)
print(sample_summary)

################################################################################
# DESCRIPTIVE STATISTICS
################################################################################

descriptive_statistics <- analysis_data %>%
  pivot_longer(
    cols = all_of(names(outcome_labels)),
    names_to = "Outcome_code",
    values_to = "Score"
  ) %>%
  mutate(Outcome = unname(outcome_labels[Outcome_code])) %>%
  group_by(Outcome, year, test) %>%
  summarise(
    n = sum(!is.na(Score)),
    Mean = mean(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE),
    .groups = "drop"
  )

print(descriptive_statistics)

################################################################################
# PRIMARY MIXED-EFFECTS MODELS
################################################################################

# Fit additive cohort-adjusted models without unsupported interactions or slopes
models <- setNames(
  lapply(names(outcome_labels), fit_outcome_model, data = analysis_data),
  names(outcome_labels)
)

for (outcome_code in names(models)) {
  cat("\n", outcome_labels[[outcome_code]], "\n", sep = "")
  print(summary(models[[outcome_code]]))
}

fixed_effects <- bind_rows(lapply(names(models), function(outcome_code) {
  tidy_fixed_effects(
    models[[outcome_code]],
    outcome_labels[[outcome_code]]
  )
}))

omnibus_tests <- bind_rows(lapply(names(models), function(outcome_code) {
  get_omnibus_tests(
    models[[outcome_code]],
    outcome_labels[[outcome_code]]
  )
}))

primary_results <- bind_rows(lapply(names(models), function(outcome_code) {
  get_primary_result(
    models[[outcome_code]],
    outcome_code,
    outcome_labels[[outcome_code]],
    analysis_data
  )
}))

model_diagnostics <- bind_rows(lapply(names(models), function(outcome_code) {
  get_model_diagnostics(
    models[[outcome_code]],
    outcome_labels[[outcome_code]]
  )
}))

manuscript_results <- primary_results %>%
  left_join(
    omnibus_tests %>%
      filter(Term == "year") %>%
      select(Outcome, Cohort_F = F, Cohort_Num_df = Num_df,
             Cohort_Den_df = Den_df, Cohort_p = p),
    by = "Outcome"
  ) %>%
  transmute(
    Outcome,
    `Adjusted pre, M [95% CI]` = sprintf(
      "%.2f [%.2f, %.2f]",
      Adjusted_pre,
      Adjusted_pre_CI_low,
      Adjusted_pre_CI_high
    ),
    `Adjusted post, M [95% CI]` = sprintf(
      "%.2f [%.2f, %.2f]",
      Adjusted_post,
      Adjusted_post_CI_low,
      Adjusted_post_CI_high
    ),
    `Adjusted change [95% CI]` = sprintf(
      "%.2f [%.2f, %.2f]",
      Adjusted_change,
      CI_low,
      CI_high
    ),
    `t (df)` = sprintf("%.2f (%.1f)", t, df),
    `Pre/post p` = ifelse(p < 0.001, "< .001", sprintf("= %.3f", p)),
    `Cohort F (df1, df2)` = sprintf(
      "%.2f (%.0f, %.1f)",
      Cohort_F,
      Cohort_Num_df,
      Cohort_Den_df
    ),
    `Cohort p` = ifelse(
      Cohort_p < 0.001,
      "< .001",
      sprintf("= %.3f", Cohort_p)
    )
  )

cat("\nPrimary adjusted pre/post results\n")
print(primary_results)

cat("\nOmnibus fixed-effect tests\n")
print(omnibus_tests)

cat("\nModel diagnostics\n")
print(model_diagnostics)

cat("\nManuscript-ready results\n")
print(manuscript_results)

################################################################################
# DIAGNOSTIC FIGURE
################################################################################

diagnostic_data <- bind_rows(lapply(names(models), function(outcome_code) {
  model <- models[[outcome_code]]
  data.frame(
    Outcome = outcome_labels[[outcome_code]],
    Fitted = fitted(model),
    Standardized_residual = residuals(model) / sigma(model)
  )
}))

residual_plot <- ggplot(
  diagnostic_data,
  aes(x = Fitted, y = Standardized_residual)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(color = "grey40", alpha = 0.55, size = 1.5) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    color = "black",
    linewidth = 0.6
  ) +
  facet_wrap(~ Outcome) +
  labs(
    title = "Residuals versus fitted values",
    x = "Fitted composite score",
    y = "Standardized residual"
  ) +
  theme_efec

qq_plot <- ggplot(
  diagnostic_data,
  aes(sample = Standardized_residual)
) +
  stat_qq(color = "grey40", alpha = 0.65, size = 1.5) +
  stat_qq_line(color = "black", linewidth = 0.6) +
  facet_wrap(~ Outcome) +
  labs(
    title = "Normal Q-Q plots",
    x = "Theoretical quantile",
    y = "Observed quantile"
  ) +
  theme_efec

diagnostic_figure <- residual_plot / qq_plot
print(diagnostic_figure)

# No strong evidence of violation of assumptions to test

ggsave(
  "Figure_HLM_Diagnostics.pdf",
  diagnostic_figure,
  width = 7,
  height = 7,
  units = "in"
)

################################################################################
# EXPORT RESULTS
################################################################################

write.csv(
  descriptive_statistics,
  "Table_HLM_Descriptive_Statistics.csv",
  row.names = FALSE
)
write.csv(
  primary_results,
  "Table_HLM_Primary_Results.csv",
  row.names = FALSE
)
write.csv(
  fixed_effects,
  "Table_HLM_Fixed_Effects.csv",
  row.names = FALSE
)
write.csv(
  omnibus_tests,
  "Table_HLM_Omnibus_Tests.csv",
  row.names = FALSE
)
write.csv(
  model_diagnostics,
  "Table_HLM_Diagnostics.csv",
  row.names = FALSE
)
write.csv(
  manuscript_results,
  "Table_HLM_Manuscript_Results.csv",
  row.names = FALSE
)

save(
  models,
  sample_summary,
  descriptive_statistics,
  primary_results,
  fixed_effects,
  omnibus_tests,
  model_diagnostics,
  manuscript_results,
  file = "HLM_Results.RData"
)
