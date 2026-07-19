################################################################################
#
# 1. MODEL FITTING
# Latest Version: July 19, 2026
# 
# This file processes desalination data for model fits
#
################################################################################

############################################
# Preliminaries 
############################################
setwd("~/Documents/Desalination")
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(patchwork)
library(knitr)

data <- read.csv("data.csv")

############################################
# Coding the factor levels
############################################
alpha_practical <- sqrt(2)

# Identify how much each axial point should have
coding <- tibble::tribble(
  ~actual, ~coded, ~center, ~one_coded_unit,
  "temp",  "A",     50,     10 / alpha_practical,
  "wind",  "B",      2,      2 / alpha_practical,
  "conc",  "C",      3,      3 / alpha_practical,
  "irrad", "D",    400,    400 / alpha_practical
)

# Create the axial points
## General function
add_coded_factors <- function(x) {
  x %>%
    mutate(
      A = (temp - 50) / (10 / alpha_practical),
      B = (wind - 2) / (2 / alpha_practical),
      C = (conc - 3) / (3 / alpha_practical),
      D = (irrad - 400) / (400 / alpha_practical)
    )
}

data <- add_coded_factors(data) %>%
  mutate(
    setting_id = interaction(
      sprintf("%.6f", A),
      sprintf("%.6f", B),
      sprintf("%.6f", C),
      sprintf("%.6f", D),
      drop = TRUE,
      lex.order = TRUE
    )
  )

############################################
# Model fits
############################################

# For rigorous comparison, test for linear and quadratic fits
# A full cubic model cannot be fit because there are only 78 points, i.e. not enough degrees of freedom
# Partial cubic models may be fit, but I have no substantive or theoretical justification to pursue it for now

# "Interaction" model contains linear main effects and all linear-by-linear interactions
# The purpose is to examine and argue for relationships between the experimental factors
# To be rigorous, explore squares in the interaction model as well ("Interaction-with-squares")

# These candidate models use coded variables
# For natural variables, rerun chosen model on temp, wind, conc and irrad instead of A,B,C,D
candidate_formulas <- list(
  linear = yield ~ A + B + C + D,
  quadratic = yield ~ (A + B + C + D)^2 +
    I(A^2) + I(B^2) + I(C^2) + I(D^2),
  interaction = yield ~ (A + B + C + D)^2 + A:B:C + B:C:D + A:B:D + A:C:D + A:B:C:D,
  interaction_with_squares = yield ~ (A + B + C + D)^2 +
    I(A^2) + I(B^2) + I(C^2) + I(D^2) + A:B:C + B:C:D + A:B:D + A:C:D + A:B:C:D
)

models <- map(candidate_formulas, ~ lm(.x, data))
summary(models$linear)
summary(models$quadratic)
summary(models$interaction) # Since ABCD is significant, retain ALL lower-order terms by hierarchy
summary(models$interaction_with_squares) # Additional support to not account for main effect squares

anova(models$interaction, models$interaction_with_squares) 
# Adding squares to the interaction model does NOT improve fit


############################################
# Establishing model adequacy
############################################
# PRESS evaluates prediction of new data
## Better models have lower PRESS
press_statistic <- function(model) {
  leverages <- hatvalues(model)
  if (any(leverages >= 1 - sqrt(.Machine$double.eps))) {
    return(NA_real_)
  }
  sum((residuals(model) / (1 - leverages))^2)
}

# There is a need to do leave-one-out cross-validation to generate relevant statistics
## Better models have lower grouped CV errors (MAE, RMSE)
## Better models have higher grouped CV predictive R-squared
grouped_cv <- function(formula, data) {
  groups <- levels(data$setting_id)
  predictions <- rep(NA_real_, nrow(data))
  
  for (group in groups) {
    test <- data$setting_id == group
    training_model <- lm(formula, data = data[!test, , drop = FALSE])
    
    if (training_model$rank < length(coef(training_model))) {
      warning(
        "Rank deficiency occurred while withholding setting ", group,
        call. = FALSE
      )
    }
    
    predictions[test] <- predict(
      training_model,
      newdata = data[test, , drop = FALSE]
    )
  }
  
  errors <- data$yield - predictions
  sst <- sum((data$yield - mean(data$yield))^2)
  
  list(
    predictions = predictions,
    rmse = sqrt(mean(errors^2)),
    mae = mean(abs(errors)),
    predictive_r_squared = 1 - sum(errors^2) / sst
  )
}

# Errors can also be in the form of sum of squares
## Better models have lower error sum of squares
pure_error_components <- function(data) {
  group_means <- ave(data$yield, data$setting_id, FUN = mean)
  pure_error_ss <- sum((data$yield - group_means)^2)
  pure_error_df <- nrow(data) - nlevels(data$setting_id)
  list(ss = pure_error_ss, df = pure_error_df)
}

# Evaluate lack-of-fit
# Do NOT choose models with a significant lack of fit, since the model does not fit the data
lack_of_fit <- function(model, data) {
  pure <- pure_error_components(data)
  residual_ss <- deviance(model)
  residual_df <- df.residual(model)
  lof_ss <- residual_ss - pure$ss
  lof_df <- residual_df - pure$df
  
  if (lof_df <= 0 || pure$df <= 0) {
    return(tibble(
      lack_of_fit_ss = NA_real_, lack_of_fit_df = NA_integer_,
      pure_error_ss = pure$ss, pure_error_df = pure$df,
      lack_of_fit_f = NA_real_, lack_of_fit_p = NA_real_
    ))
  }
  
  f_value <- (lof_ss / lof_df) / (pure$ss / pure$df)
  tibble(
    lack_of_fit_ss = lof_ss,
    lack_of_fit_df = lof_df,
    pure_error_ss = pure$ss,
    pure_error_df = pure$df,
    lack_of_fit_f = f_value,
    lack_of_fit_p = pf(f_value, lof_df, pure$df, lower.tail = FALSE)
  )
}

# For completeness, evaluate total sum-of-squares
sst <- sum((data$yield - mean(data$yield))^2)

# Generate grouped CV results
grouped_cv_results <- imap(
  candidate_formulas,
  ~ grouped_cv(.x, data)
)

# Tabulate all the relevant statistics
model_comparison <- imap_dfr(models, function(fit_object, model_name) {
  fit <- glance(fit_object)
  press <- press_statistic(fit_object)
  grouped <- grouped_cv_results[[model_name]]
  tibble(
    model = model_name,
    coefficients = fit_object$rank,
    residual_df = df.residual(fit_object),
    residual_standard_error = sigma(fit_object),
    r_squared = fit$r.squared,
    adjusted_r_squared = fit$adj.r.squared,
    AIC = AIC(fit_object),
    BIC = BIC(fit_object),
    PRESS = press,
    observation_press_r_squared = 1 - press / sst,
    grouped_cv_RMSE = grouped$rmse,
    grouped_cv_MAE = grouped$mae,
    grouped_cv_predictive_r_squared = grouped$predictive_r_squared
  ) %>%
    bind_cols(lack_of_fit(fit_object, data))
})

# This transposes the table for HTML export
comparison_display <- model_comparison |>
  mutate(
    model = recode(
      model,
      linear = "Linear",
      quadratic = "Quadratic",
      interaction = "Interaction",
      interaction_with_squares = "Interaction + squares"
    ),
    across(where(is.numeric), \(x) round(x, 3))
  ) |>
  pivot_longer(
    cols = -model,
    names_to = "Statistic",
    values_to = "Value"
  ) |>
  pivot_wider(
    names_from = model,
    values_from = Value
  ) |>
  mutate(
    Statistic = gsub("_", " ", Statistic)
  )

kable(
  comparison_display,
  align = c("l", "r", "r", "r", "r"),
  caption = "Model comparison"
)
# It seems reasonable to choose interaction model without dropping anything

############################################
# Diagnostic plots
############################################
# This provides residuals and fitted values to original design
diagnostics <- augment(models$interaction, data) %>%
  arrange(run.order)

# This generates confidence intervals for coefficients
coefficients <- tidy(models$interaction, conf.int = TRUE, conf.level = 0.95)

# Observed vs Fitted plot
# Observed data should roughly follow trend of model-Fitted data
diagnostic_observed <- ggplot(diagnostics, aes(.fitted, yield)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_point(
    size = 2.2,
    alpha = 0.85,
    colour = "#2C7FB8"
  ) +
  coord_equal() +
  labs(x = "Predicted yield (g)", y = "Observed yield (g)")

# Residuals vs Fitted plot
# Residuals should not have a trend across Fitted values (homoscedasticity)
diagnostic_residuals <- ggplot(diagnostics, aes(.fitted, .std.resid)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_hline(
    yintercept = c(-2, 2), linetype = "dashed", linewidth = 0.4,
    colour = "grey55"
  ) +
  geom_point(size = 2.1, alpha = 0.8, colour = "#2C7FB8") +
  labs(x = "Fitted yield (g)", y = "Standardized residual")

# Normal Q-Q Plot
# Points should roughly follow the Q-Q line (normality)
diagnostic_qq <- ggplot(diagnostics, aes(sample = .std.resid)) +
  stat_qq(size = 2, alpha = 0.8, colour = "#2C7FB8") +
  stat_qq_line(linewidth = 0.6, colour = "grey35") +
  labs(x = "Theoretical quantile", y = "Standardized residual")

# Run Order Plot
# There should be no trend across run order (independence)
diagnostic_run <- ggplot(diagnostics, aes(run.order, .std.resid)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_line(linewidth = 0.35, colour = "grey65") +
  geom_point(size = 2, colour = "#2C7FB8") +
  labs(x = "Randomized run order", y = "Standardized residual")

# Just some theme options for ggplot
diagnostic_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    legend.position = "right"
  )

# Generates 2x2 plot for assumption checking
diagnostic_figure <- (
  diagnostic_observed + diagnostic_residuals +
    diagnostic_qq + diagnostic_run
) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste("Model diagnostics:", models$interaction),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  ) & diagnostic_theme

diagnostic_figure
# No strong evidence to violations of assumptions of test
