################################################################################
#
# 3. CONDITIONAL SLOPE PLOTS
# Latest Version: July 19, 2026
# 
# This file generates plots to more clearly display the effect of CaCO3
#
################################################################################

############################################
# Preliminaries 
############################################
library(ggplot2)
library(dplyr)

data <- read.csv("full_data.csv")

# We need to respecify this for the plots
alpha_practical <- sqrt(2)
add_coded_factors <- function(x) {
  x %>%
    mutate(
      A = (temp - 50) / (10 / alpha_practical),
      B = (wind - 2) / (2 / alpha_practical),
      C = (conc - 3) / (3 / alpha_practical),
      D = (irrad - 400) / (400 / alpha_practical)
    )
}

# Specify model
model_formula <- yield ~
  (A + B + C + D)^2 +
  A:B:C + A:B:D + A:C:D + B:C:D + A:B:C:D

model <- lm(model_formula, data = data)

############################################
# CaCO3 plots
############################################
# Specify low/high factorial slices plus center settings
wind_slices <- c("Low wind = 0.59 m/s" = 0.59,
                 "Centered wind = 2.00 m/s" = 2.00,
                 "High wind = 3.41 m/s" = 3.41)
irrad_slices <- c("Low IR = 117 W/m^2" = 117.16,
                  "Centered IR = 400 W/m^2" = 400.00,
                  "High IR = 683 W/m^2" = 682.84)

effect_grid <- expand.grid(
  temp = seq(42.93, 57.07, length.out = 241),
  conc = 3,
  wind_label = names(wind_slices),
  irrad_label = names(irrad_slices),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    wind = unname(wind_slices[wind_label]),
    irrad = unname(irrad_slices[irrad_label]),
    wind_label = factor(wind_label, levels = names(wind_slices)),
    irrad_label = factor(irrad_label, levels = names(irrad_slices))
  ) %>%
  add_coded_factors()

conditional_slope <- function(model, newdata, variable = "conc", epsilon = 1e-5) {
  plus <- newdata
  minus <- newdata
  plus[[variable]] <- plus[[variable]] + epsilon
  minus[[variable]] <- minus[[variable]] - epsilon
  plus <- add_coded_factors(plus)
  minus <- add_coded_factors(minus)
  model_terms <- delete.response(terms(model))
  x_plus <- model.matrix(model_terms, data = plus)
  x_minus <- model.matrix(model_terms, data = minus)
  contrast <- (x_plus - x_minus) / (2 * epsilon)
  estimate <- drop(contrast %*% coef(model))
  covariance <- vcov(model)
  standard_error <- sqrt(rowSums((contrast %*% covariance) * contrast))
  data.frame(
    caco3_slope = estimate,
    caco3_slope_se = standard_error
  )
}

effect_grid <- bind_cols(
  effect_grid,
  conditional_slope(model, effect_grid)
)

# I'm assuming we're using 95% confidence level
t_critical <- qt(0.975, df.residual(model))
effect_grid <- effect_grid %>%
  mutate(
    lower_95 = caco3_slope - t_critical * caco3_slope_se,
    upper_95 = caco3_slope + t_critical * caco3_slope_se,
    support = case_when(
      lower_95 > 0 ~ "Positive",
      upper_95 < 0 ~ "Negative",
      TRUE ~ "Not distinguished from zero"
    )
  )

plot_conditional <- ggplot(
  effect_grid,
  aes(x = temp, y = caco3_slope)
) +
  geom_ribbon(
    aes(ymin = lower_95, ymax = upper_95),
    fill = "grey80",
    colour = NA
  ) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_line(colour = "black", linewidth = 0.75) +
  facet_grid(
    rows = vars(wind_label),
    cols = vars(irrad_label),
    switch = "y"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(breaks = c(43, 47, 50, 53, 57)) +
  scale_y_continuous(
    breaks = seq(-80, 40, 20),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  labs(
    x = expression("Temperature (" * degree * "C)"),
    y = expression("Conditional CaCO"[3] * " slope (g yield per g CaCO"[3] * ")"),
    caption = 
      "Line: estimated conditional CaCO3 slope; grey band: pointwise 95% confidence interval.\n"
    ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.spacing = grid::unit(1.2, "lines"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    axis.text.x = element_text(margin = margin(t = 4)),
    axis.text.y = element_text(margin = margin(r = 4)),
    plot.caption = element_text(hjust = 0, size = 9, margin = margin(t = 8))
  )

plot_conditional
