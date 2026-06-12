##################################################################
#
# Various plots for simulated data
# Use and reuse as needed
#
##################################################################

library(ggplot2)
library(readxl)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(forcats)
library(patchwork)
library(metR)

setwd("~/Documents/Longitudinal/New")
data <- read_excel("outputs.xlsx", sheet = "Sheet 1")

# non-coverage
noncov <- data[data$ICC %in% c(0.1, 0.2, 0.3) & data$N %in% c(10, 20, 30) & data$T %in% c(3, 4, 5) & 
                 data$pE == 0.4, 
               c("var_epsilon", "b_11_prop","ICC", "N", "T", "b_11_noncov", "b_10_noncov", "b_00_noncov", "tau_00_sq_noncov", "tau_11_sq_noncov", "sigma_sq_noncov")]
noncov <- noncov[order(noncov$ICC, noncov$N, noncov$T), ]
write.csv(noncov,"noncov_lowvar_highprop.csv")

# inadequacy
inadequacy <- data[data$ICC %in% c(0.1, 0.2, 0.3) & data$N %in% c(10, 20, 30) & data$T %in% c(3, 4, 5) & 
                     data$pE == 0.4, c("var_epsilon", "b_11_prop","ICC", "N", "T", "conv_problem")]
inadequacy <- inadequacy[order(inadequacy$ICC, inadequacy$N, inadequacy$T), ]
write.csv(inadequacy,"inadequacy_lowvar_highprop.csv")

# b_11 Power Plots
facet_labels <- c("10" = "N=10", "20" = "N=20", "30" = "N=30")
ggplot(data[data$N %in% c(10, 20, 30) & data$T %in% c(3, 4, 5) & data$ICC %in% c(0.1, 0.5, 0.9), ], 
       aes(x = pE, y = b_11_power, color = factor(ICC), shape = factor(T), linetype = factor(T))) +
  geom_smooth(aes(group = interaction(ICC, T)), method = "loess", se = FALSE, size = 0.6) + 
  geom_point(size = 2, stroke = 0.5, fill = "white") +  
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +  
  scale_shape_manual(values = c(NA, NA, NA)) + 
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +  
  labs(y = "power", color = "ICC", shape = "T", linetype = "T") +
  theme_minimal() +
  theme(
    legend.position = c(0.93, 0.15),  
    legend.justification = c(1, 0),
    legend.box = "horizontal", 
    legend.box.just = "right",  
    legend.spacing.x = unit(.5, "cm"),  
    legend.background = element_rect(fill = "white", color = "black"), 
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  coord_cartesian(ylim = c(0.5, 1)) +
  facet_wrap(~ N, ncol = 3, labeller = as_labeller(facet_labels)) +
  guides(
    color = guide_legend(ncol = 1, title.position = "top"),  # ICC legend in a single column
    shape = guide_legend(ncol = 1, title.position = "top")   # T legend in a single column next to ICC
  )

# Power contour plots
facet_N <- c(`10` = "N = 10", `20` = "N = 20", `30` = "N = 30")
facet_T <- c(`3`  = "T = 3",  `4`  = "T = 4",  `5`  = "T = 5")

solid_breaks  <- setdiff(seq(0.10, 0.90, 0.10), 0.80)   
label_breaks  <- seq(0.50, 0.90, 0.10)                

ggplot(data, aes(pE, ICC, z = b_11_power)) +
  geom_contour_filled(
    breaks      = c(-Inf, 0.80, Inf),   
    show.legend = FALSE,
    alpha       = 1                     
  ) +
  scale_fill_manual(
    values = c(          
      "white",        
      "grey75"          
    )
  ) +
  geom_contour(breaks = solid_breaks,
               colour = "black", linewidth = 0.15) +
  geom_label_contour(
    breaks        = label_breaks,
    label.placer  = label_placer_flattest(),
    fill          = "white",
    label.padding = unit(0.12, "lines"),
    aes(label = after_stat(sprintf("%.2f", level)))
  ) +
  geom_contour(breaks = 0.80,
               colour = "black",
               linetype = "dashed",
               linewidth = 0.6) +
  facet_grid(rows = vars(T), cols = vars(N),
             labeller = labeller(N = facet_N, T = facet_T),
             switch   = "y") +           
  coord_cartesian(expand = 0) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  labs(x = "pE", y = "ICC") +
  theme_minimal(base_size = 11) +
  theme(
    panel.spacing     = unit(1.2, "lines"),
    axis.text.x       = element_text(margin = margin(t = 4)),
    axis.text.y       = element_text(margin = margin(r = 4)),
    panel.grid        = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 0),
    legend.position   = "none"
  )

# RMSEs
rmse <- data[data$ICC %in% c(0.1, 0.2, 0.3) & data$N %in% c(10, 20, 30) & data$T %in% c(3, 4, 5) & 
        data$pE == 0.4, 
        c("ICC", "N", "T", "b_11_rmse", "b_10_rmse", "b_00_rmse", "tau_00_sq_rmse", "tau_11_sq_rmse", "sigma_sq_rmse")]
rmse <- rmse[order(rmse$ICC, rmse$N, rmse$T), ]
write.csv(rmse,"rmse.csv")

# Function to extract RMSE values based on input parameters
extract_rmse_values <- function(df, N_val, T_val, pE_val, ICC_val) {
  df %>%
    dplyr::filter(
      N == N_val,
      T == T_val,
      pE == pE_val,
      ICC == ICC_val
    ) %>%
    dplyr::select(
      b_11_power,
      b_00_rmse,
      b_10_rmse,
      b_11_rmse,
      tau_00_sq_rmse,
      tau_11_sq_rmse,
      sigma_sq_rmse
    )
}

# Baseline
extract_rmse_values(data, N_val=30, T_val=5, pE_val = 0.5, ICC_val = 0.5)

# Small N
extract_rmse_values(data, N_val=10, T_val=5, pE_val = 0.3, ICC_val = 0.5)

# Small T
extract_rmse_values(data, N_val=30, T_val=3, pE_val = 0.3, ICC_val = 0.6)

# Max ICC, N = 20
extract_rmse_values(data, N_val=20, T_val=5, pE_val = 0.5, ICC_val = 0.9)

# Max ICC, N = 30
extract_rmse_values(data, N_val=30, T_val=3, pE_val = 0.4, ICC_val = 0.9)

# RMSE heatmap
rmse_columns <- c("b_11_rmse", "b_10_rmse", "b_00_rmse", "tau_00_sq_rmse", "tau_11_sq_rmse", "sigma_sq_rmse")
data_filtered <- data[data$pE == 0.3, c("N", "T", "ICC", "b_11_rmse", "b_10_rmse", "b_00_rmse", "tau_00_sq_rmse", "tau_11_sq_rmse", "sigma_sq_rmse")]
data_long <- data_filtered %>%
  pivot_longer(cols = all_of(rmse_columns), names_to = "RMSE_Type", values_to = "RMSE_Value")

data_long <- data_long %>%
  mutate(Combined = paste("ICC:", ICC, "| N:", N, "| T:", T))

data_long <- data_long %>%
  arrange(ICC, N, T) %>%  # Sort by ICC, N, and T in ascending order
  mutate(Combined = factor(Combined, levels = unique(Combined)))

data_long$RMSE_Type <- factor(data_long$RMSE_Type, levels = c(
  "b_00_rmse", "b_10_rmse", "b_11_rmse", 
  "sigma_sq_rmse", "tau_00_sq_rmse", "tau_11_sq_rmse"),
  labels = c(
    "b[00]", "b[10]", "b[11]",
    "sigma^2", "tau[00]^2", "tau[11]^2"
  )
)

data_long$Combined <- fct_rev(factor(data_long$Combined))

ggplot(data_long, aes(x = RMSE_Type, y = Combined, fill = RMSE_Value)) +
  geom_tile(color = "white") +  
  geom_text(aes(label = sprintf("%.4f", RMSE_Value)), size = 3, hjust = 0.5, vjust = 0.5) + 
  scale_fill_gradient(low = "lightblue", high = "red", name = "RMSE") +  
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), 
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 12, face = "italic"), 
    legend.position = "bottom", 
    panel.spacing = unit(-0.5, "lines") 
  ) +
  facet_wrap(~ RMSE_Type, ncol = 6, scales = "free_x", labeller = label_parsed)


########################
# Checking if var or prop of b_11 affects
# the convergence problems
########################

inadequacy_lowvar_lowprop <- read.csv("inadequacy_lowvar_lowprop.csv")
inadequacy_lowvar_highprop <- read.csv("inadequacy_lowvar_highprop.csv")
inadequacy_highvar_lowprop <- read.csv("inadequacy_highvar_lowprop.csv")
inadequacy_highvar_highprop <- read.csv("inadequacy_highvar_highprop.csv")
inadequacy_results <- rbind(inadequacy_lowvar_lowprop, inadequacy_lowvar_highprop, inadequacy_highvar_lowprop, inadequacy_highvar_highprop)
View(inadequacy_results)

noncov_lowvar_lowprop <- read.csv("noncov_lowvar_lowprop.csv")
noncov_lowvar_highprop <- read.csv("noncov_lowvar_highprop.csv")
noncov_highvar_lowprop <- read.csv("noncov_highvar_lowprop.csv")
noncov_highvar_highprop <- read.csv("noncov_highvar_highprop.csv")
noncov_results <- rbind(noncov_lowvar_lowprop, noncov_lowvar_highprop, noncov_highvar_lowprop, noncov_highvar_highprop)
View(noncov_results)

########################
# Plotting attrition
########################

library(ggplot2)
library(dplyr)
library(openxlsx)
library(gridExtra)

attr <- read.xlsx("outputs_phase2.new.xlsx")

filtered_data_1 <- attr %>%
  dplyr::filter(name %in% c("input_1", "input_2", "input_3", "input_4", "input_5",
                            "input_17", "input_18", "input_19", "input_20", "input_21"))

filtered_data_2 <- attr %>%
  dplyr::filter(name %in% c("input_33", "input_34", "input_35", "input_36", "input_37",
                            "input_49", "input_50", "input_51", "input_52", "input_53"))

plot_1 <- ggplot(filtered_data_1, aes(x = attrition, y = b_11_power, color = as.factor(ICC), group = ICC)) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +
  labs(
    x = "Attrition",
    y = "Power",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_color_manual(values = c("0.1" = "#1b9e77", "0.3" = "#7570b3")) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.2, 1)) +
  annotate("text", x = 0.05, y = 0.22, label = "N = 20, T = 5", hjust = 0, vjust = 0, size = 4, color = "black")

plot_2 <- ggplot(filtered_data_2, aes(x = attrition, y = b_11_power, color = as.factor(ICC), group = ICC)) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +
  labs(
    x = "Attrition",
    y = "",
    color = "ICC"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.85),
    text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_color_manual(values = c("0.1" = "#1b9e77", "0.3" = "#7570b3")) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.2, 1)) +
  annotate("text", x = 0.05, y = 0.22, label = "N = 30, T = 3", hjust = 0, vjust = 0, size = 4, color = "black")

grid.arrange(plot_1, plot_2, ncol = 2)


########################
# Plotting nonsampling errors
########################

library(ggplot2)
library(dplyr)
library(openxlsx)
library(gridExtra)

nonsamp <- read.xlsx("outputs_final.xlsx", sheet = "Phase 2")

# Filter data for N = 20, T = 5
filtered_data_1 <- attr %>%
  dplyr::filter(
    (name %in% c("input_nonsamp_a","input_nonsamp_b","input_nonsamp_c",
                "input_1_nonsamp", "input_2_nonsamp", "input_3_nonsamp", "input_4_nonsamp", 
                 "input_5_nonsamp", "input_6_nonsamp", "input_7_nonsamp", "input_8_nonsamp", 
                 "input_9_nonsamp", "input_10_nonsamp", "input_11_nonsamp", "input_12_nonsamp") & ICC == 0.1) |
      (name %in% c("input_nonsamp_d","input_nonsamp_e","input_nonsamp_f",
                   "input_13_nonsamp", "input_14_nonsamp", "input_15_nonsamp", "input_16_nonsamp", 
                   "input_17_nonsamp", "input_18_nonsamp", "input_19_nonsamp", "input_20_nonsamp", 
                   "input_21_nonsamp", "input_22_nonsamp", "input_23_nonsamp", "input_24_nonsamp") & ICC == 0.3)
  )

# Filter data for N = 30, T = 3
filtered_data_2 <- attr %>%
  dplyr::filter(
    (name %in% c("input_nonsamp_A","input_nonsamp_B","input_nonsamp_C",
                 "input_25_nonsamp", "input_26_nonsamp", "input_27_nonsamp", "input_28_nonsamp", 
                 "input_29_nonsamp", "input_30_nonsamp", "input_31_nonsamp", "input_32_nonsamp", 
                 "input_33_nonsamp", "input_34_nonsamp", "input_35_nonsamp", "input_36_nonsamp") & ICC == 0.1) |
      (name %in% c("input_nonsamp_D","input_nonsamp_E","input_nonsamp_F",
                   "input_37_nonsamp", "input_38_nonsamp", "input_39_nonsamp", "input_40_nonsamp", 
                   "input_41_nonsamp", "input_42_nonsamp", "input_43_nonsamp", "input_44_nonsamp", 
                   "input_45_nonsamp", "input_46_nonsamp", "input_47_nonsamp", "input_48_nonsamp") & ICC == 0.3)
  )

# Create the first plot (N = 20, T = 5)
plot_1 <- ggplot(filtered_data_1, aes(x = nonsamp_p, y = b_11_power, 
                                      color = as.factor(ICC), shape = as.factor(nonsamp_q), group = interaction(ICC, nonsamp_q))) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +
  labs(
    x = "Non-sampling p",
    y = "Power",
    color = "ICC",
    shape = "Non-sampling q"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_color_manual(values = c("0.1" = "#1b9e77", "0.3" = "#7570b3")) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.6, 1)) +
  annotate("text", x = 0.05, y = 0.6, label = "N = 20, T = 5", hjust = 0, vjust = 0, size = 4, color = "black")

# Create the second plot (N = 30, T = 3)
# Create the second plot (N = 30, T = 3) with a side-by-side legend
plot_2 <- ggplot(filtered_data_2, aes(x = nonsamp_p, y = b_11_power, 
                                      color = as.factor(ICC), shape = as.factor(nonsamp_q), group = interaction(ICC, nonsamp_q))) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +
  labs(
    x = "Non-sampling p",
    y = "",
    color = "ICC",
    shape = "Non-sampling q"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.35, 0.2),
    legend.box = "horizontal",  # Arrange legends side by side
    text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_color_manual(values = c("0.1" = "#1b9e77", "0.3" = "#7570b3")) +
  guides(
    color = guide_legend(nrow = 1),  # Arrange legend items in a single row
    shape = guide_legend(nrow = 1)   # Arrange legend items in a single row
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.6, 1)) +
  annotate("text", x = 0.05, y = 0.6, label = "N = 30, T = 3", hjust = 0, vjust = 0, size = 4, color = "black")

grid.arrange(plot_1, plot_2, ncol = 2)


########################
# Plotting inadmissibility
########################

inadmiss <- read.xlsx("outputs_final.xlsx", sheet = "Inadmissibility")

plot <- ggplot(inadmiss, aes(x = T, y = conv_problem, color = as.factor(N), group = as.factor(N))) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = c(0.05, 0.25), linetype = "dashed", color = "black") +
  annotate("text", x = 3, y = 0.05, label = "0.05", vjust = -1, size = 4) +
  annotate("text", x = 3, y = 0.25, label = "0.25", vjust = -1, size = 4) +
  scale_x_continuous(breaks = seq(3, 10, by = 1)) +  # Ensures x-axis ticks from 3 to 10
  labs(
    x = "T",
    y = "Proportion",
    color = "N"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = c(0.75, 0.9),  # Position legend in the top-right of the plot
    legend.direction = "horizontal",  # Make the legend horizontal
    legend.box = "horizontal"  # Align legend items horizontally
  )

# Print the plot
print(plot)

