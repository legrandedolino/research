################################################################################
#
# 2. CONTOUR PLOTS
# Latest Version: July 19, 2026
# 
# This file generates contour plots from the appropriate model
#
################################################################################

############################################
# Preliminaries 
############################################
library(ggplot2)
library(dplyr)
library(metR)

data <- read.csv("full_data.csv")

# We need to respecify this for the contour plots
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

# This populates values for the contour plot
surface_grid <- expand.grid(
  conc = seq(0.88, 5.12, length.out = 161),
  temp = seq(42.93, 57.07, length.out = 161),
  wind_label = names(wind_slices),
  irrad_label = names(irrad_slices),
  KEEP.OUT.ATTRS = FALSE
) %>%
  mutate(
    wind = unname(wind_slices[wind_label]),
    irrad = unname(irrad_slices[irrad_label]),
    wind_label = factor(wind_label, levels = names(wind_slices)),
    irrad_label = factor(irrad_label, levels = names(irrad_slices))
  ) %>%
  add_coded_factors()

surface_grid$predicted_yield <- predict(model, newdata = surface_grid)

#These specify the breaks of the plot
yield_range <- range(surface_grid$predicted_yield, finite = TRUE)
contour_step <- 25
solid_breaks <- seq(
  floor(yield_range[1] / contour_step) * contour_step,
  ceiling(yield_range[2] / contour_step) * contour_step,
  by = contour_step
)
label_breaks <- solid_breaks[solid_breaks %% 50 == 0]
fill_breaks <- seq(
  floor(yield_range[1] / 50) * 50,
  ceiling(yield_range[2] / 50) * 50,
  by = 50
)

plot_caco3 <- ggplot(surface_grid, aes(x = conc, y = temp)) +
  geom_contour_filled(
    aes(z = predicted_yield),
    breaks = fill_breaks,
    show.legend = FALSE
  ) +
  scale_fill_grey(start = 0.98, end = 0.25) +
  geom_contour(
    aes(z = predicted_yield),
    breaks = solid_breaks,
    colour = "black",
    linewidth = 0.20
  ) +
  geom_label_contour(
    aes(
      z = predicted_yield,
      label = after_stat(sprintf("%.0f", level))
    ),
    breaks = label_breaks,
    label.placer = label_placer_fraction(0.50),
    fill = "white",
    label.padding = grid::unit(0.10, "lines"),
    size = 3.1,
    label.size = 0.15,
    skip = 0
  ) +
  facet_grid(
    rows = vars(wind_label),
    cols = vars(irrad_label),
    switch = "y"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(breaks = seq(1, 5, 1)) +
  scale_y_continuous(breaks = seq(43, 57, 2)) +
  labs(
    x = expression(CaCO[3] ~ "concentration (g)"),
    y = expression("Temperature (" * degree * "C)"),
    caption = "Contours and shading: predicted water yield (g per 6 h); higher yield is darker.\n"
    ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.spacing = grid::unit(1.2, "lines"),
    axis.text.x = element_text(margin = margin(t = 4)),
    axis.text.y = element_text(margin = margin(r = 4)),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    plot.caption = element_text(hjust = 0, size = 9, margin = margin(t = 8))
  )

plot_caco3
