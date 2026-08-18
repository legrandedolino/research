################################################################################
#
# 4. FACTOR ANALYSIS - EDUCATING FOR ENVIRONMENTAL CHANGE
# Latest Version: August 18, 2026
#
# Validates scales, calculates composite scores, and creates pre/post figures.
#
################################################################################

################################################################################
# PACKAGES AND SHARED SETTINGS
################################################################################

library(psych)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
library(patchwork)
library(purrr)

reverse_efficacy_items <- c(
  "Efficacy_2", "Efficacy_4", "Efficacy_5", "Efficacy_7",
  "Efficacy_8", "Efficacy_9", "Efficacy_11"
)

efficacy_items <- paste0("Efficacy_", 1:11)
efficacy_items_scored <- setdiff(efficacy_items, "Efficacy_10")
understanding_items <- c(
  "ccUnder_climChange",
  "ccUnder_humanCaused",
  "ccUnder_policy",
  "ccUnder_stepsPeople",
  "ccUnder_mitigateEffects"
)

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
    plot.subtitle = element_text(color = "grey30"),
    plot.margin = margin(6, 10, 6, 6),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

################################################################################
# HELPERS
################################################################################

# Preserve all-missing row means as NA instead of NaN
row_mean_or_na <- function(data) {
  result <- rowMeans(data, na.rm = TRUE)
  result[rowSums(!is.na(data)) == 0] <- NA_real_
  result
}

run_factor_checks <- function(data, label, factors = 1, run_parallel = TRUE) {
  message("\n", label)
  print(psych::KMO(data))
  print(psych::cortest.bartlett(cor(data), n = nrow(data)))

  if (run_parallel) {
    psych::fa.parallel(
      data,
      fm = "pa",
      fa = "fa",
      main = paste(label, "parallel analysis")
    )
  }

  model <- psych::fa(
    data,
    nfactors = factors,
    rotate = "oblimin",
    fm = "pa"
  )
  print(model, cut = 0.3, sort = TRUE)

  reliability <- psych::alpha(data)
  print(reliability)

  invisible(list(model = model, reliability = reliability))
}

make_paired_items <- function(data, variables, excluded_years = NULL) {
  paired <- data

  if (!is.null(excluded_years)) {
    paired <- paired %>% filter(!year %in% excluded_years)
  }

  paired %>%
    select(ID, year, test, all_of(variables)) %>%
    pivot_longer(
      cols = all_of(variables),
      names_to = "Topic",
      values_to = "Score"
    ) %>%
    mutate(Score = as.numeric(as.character(Score))) %>%
    pivot_wider(names_from = test, values_from = Score) %>%
    drop_na(pre, post)
}

summarize_increases <- function(data, labels = NULL) {
  result <- data %>%
    mutate(increased = post > pre) %>%
    group_by(Topic, year) %>%
    summarise(
      total_pairs = n(),
      increased_n = sum(increased),
      increased_pct = increased_n / total_pairs,
      .groups = "drop"
    )

  if (!is.null(labels)) {
    result <- result %>% mutate(Topic = unname(labels[Topic]))
  }

  result %>%
    mutate(increased_pct = scales::percent(increased_pct, accuracy = 0.1)) %>%
    arrange(Topic, year)
}

# Calculate paired Wilcoxon tests and p-values by item and cohort
summarize_wilcoxon <- function(data, labels = NULL) {
  result <- data %>%
    group_by(Topic, year) %>%
    group_modify(~ {
      differences <- .x$post - .x$pre
      nonzero <- differences[differences != 0]

      if (length(nonzero) == 0) {
        return(tibble(
          n = nrow(.x),
          V = 0,
          p_raw = 1,
          rank_biserial = 0
        ))
      }

      test <- wilcox.test(
        .x$post,
        .x$pre,
        paired = TRUE,
        exact = FALSE,
        correct = TRUE
      )
      signed_ranks <- rank(abs(nonzero), ties.method = "average")
      positive_sum <- sum(signed_ranks[nonzero > 0])
      negative_sum <- sum(signed_ranks[nonzero < 0])

      tibble(
        n = nrow(.x),
        V = unname(test$statistic),
        p_raw = test$p.value,
        rank_biserial =
          (positive_sum - negative_sum) / (positive_sum + negative_sum)
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adjusted = p.adjust(p_raw, method = "holm"),
      significance = case_when(
        p_adjusted < 0.001 ~ "***",
        p_adjusted < 0.01 ~ "**",
        p_adjusted < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  if (!is.null(labels)) {
    result <- result %>% mutate(Item = unname(labels[Topic]))
  } else {
    result <- result %>% mutate(Item = Topic)
  }

  result %>%
    select(Item, Topic, year, n, V, p_raw, p_adjusted, rank_biserial, significance) %>%
    arrange(Item, year)
}

make_likert_plot_list <- function(
    paired_data,
    response_levels,
    response_labels,
    title_labels = NULL,
    wilcoxon_results = NULL) {
  long_data <- paired_data %>%
    pivot_longer(
      cols = c(pre, post),
      names_to = "test",
      values_to = "Score"
    ) %>%
    mutate(
      year = factor(year),
      test = factor(test, levels = c("post", "pre")),
      Score = factor(Score, levels = response_levels)
    )

  long_data %>%
    split(.$Topic) %>%
    imap(function(data_subset, topic_name) {
      clean_title <- if (is.null(title_labels)) {
        topic_name %>%
          str_remove("^[^_]+_") %>%
          str_replace_all("([a-z])([A-Z])", "\\1 \\2") %>%
          str_to_title()
      } else {
        title_labels[[topic_name]]
      }

      summary_data <- data_subset %>%
        count(year, test, Score, name = "Count") %>%
        group_by(year, test) %>%
        mutate(Percent = Count / sum(Count)) %>%
        ungroup()

      year_annotations <- tibble(
        year = levels(data_subset$year),
        year_label = levels(data_subset$year)
      )

      if (!is.null(wilcoxon_results)) {
        year_annotations <- wilcoxon_results %>%
          filter(Topic == topic_name) %>%
          transmute(
            year = as.character(year),
            year_label = str_trim(
              paste(year, if_else(significance == "ns", "", significance))
            )
          )
      }

      summary_data <- summary_data %>%
        mutate(year = as.character(year)) %>%
        left_join(year_annotations, by = "year") %>%
        mutate(
          year_label = factor(
            year_label,
            levels = year_annotations$year_label
          )
        )

      plot <- ggplot(summary_data, aes(x = Percent, y = test, fill = Score)) +
        geom_col(
          position = position_fill(reverse = TRUE),
          color = "black",
          width = 0.7
        ) +
        geom_text(
          aes(
            label = ifelse(
              Percent > 0.05,
              scales::percent(Percent, accuracy = 1),
              ""
            )
          ),
          position = position_fill(vjust = 0.5, reverse = TRUE),
          size = 3.2,
          fontface = "bold",
          color = "black"
        ) +
        facet_grid(year_label ~ ., switch = "y") +
        scale_fill_brewer(
          palette = "RdYlBu",
          drop = FALSE,
          labels = response_labels
        ) +
        scale_x_continuous(labels = scales::percent_format()) +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        labs(
          title = clean_title,
          x = "Percentage of cohort",
          y = "Cohort year",
          fill = "Response Score"
        ) +
        theme_efec +
        theme(
          axis.text.y = element_text(face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"),
          strip.placement = "outside"
        )

      plot
    })
}

combine_likert_plots <- function(plot_list, ncol) {
  common_legend <- cowplot::get_legend(
    plot_list[[1]] +
      theme(
        legend.position = "bottom",
        legend.box = "vertical"
      )
  )

  panel_grid <- wrap_plots(
    map(plot_list, ~ .x + theme(legend.position = "none")),
    ncol = ncol
  )

  wilcoxon_legend <- ggplot() +
    annotate(
      "text",
      x = 0.02,
      y = 0.5,
      label = "Wilcoxon Test",
      hjust = 0,
      size = 3.2
    ) +
    annotate(
      "text",
      x = 0.49,
      y = c(0.80, 0.50, 0.20),
      label = c("***", "**", "*"),
      hjust = 0,
      size = 3.2
    ) +
    annotate(
      "text",
      x = 0.67,
      y = c(0.80, 0.50, 0.20),
      label = c("p < 0.001", "p < 0.01", "p < 0.05"),
      hjust = 0,
      size = 3.2
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()

  legend_row <- (wrap_elements(full = common_legend) | wilcoxon_legend) +
    plot_layout(widths = c(0.76, 0.24))

  (panel_grid / legend_row) +
    plot_layout(heights = c(1, 0.12))
}

summarize_change_scores <- function(data) {
  data %>%
    group_by(Outcome, year, year_index) %>%
    summarise(
      n = n(),
      mean_change = mean(Change),
      sd_change = sd(Change),
      se_change = sd_change / sqrt(n),
      lower = mean_change - qt(0.975, df = n - 1) * se_change,
      upper = mean_change + qt(0.975, df = n - 1) * se_change,
      .groups = "drop"
    )
}

make_violin_plot <- function(data, summary_data, outcome, title, y_label = NULL) {
  raw_subset <- data %>% filter(Outcome == outcome)
  summary_subset <- summary_data %>% filter(Outcome == outcome)

  year_labels <- setNames(
    paste0(summary_subset$year, "\n(n = ", summary_subset$n, ")"),
    summary_subset$year_index
  )

  ggplot(
    raw_subset,
    aes(x = year_index, y = Change, group = year)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey35",
      linewidth = 0.6
    ) +
    geom_violin(
      width = 0.78,
      trim = FALSE,
      adjust = 0.9,
      fill = "grey95",
      color = "grey65",
      linewidth = 0.55
    ) +
    geom_point(
      position = position_jitter(width = 0.09, height = 0, seed = 2026),
      alpha = 0.5,
      size = 1.7,
      color = "grey45"
    ) +
    geom_errorbar(
      data = summary_subset,
      aes(
        x = year_index,
        ymin = lower,
        ymax = upper
      ),
      inherit.aes = FALSE,
      width = 0.12,
      linewidth = 0.75,
      color = "black"
    ) +
    geom_point(
      data = summary_subset,
      aes(x = year_index, y = mean_change),
      inherit.aes = FALSE,
      shape = 21,
      fill = "black",
      color = "black",
      stroke = 0.7,
      size = 3.6
    ) +
    scale_x_continuous(
      breaks = summary_subset$year_index,
      labels = year_labels,
      expand = expansion(add = c(0.4, 0.4))
    ) +
    scale_y_continuous(breaks = seq(-1, 5, by = 1)) +
    coord_cartesian(ylim = c(-2, 6)) +
    labs(
      title = title,
      x = "Cohort Year",
      y = y_label
    ) +
    theme_efec
}

################################################################################
# LOAD DATA
################################################################################

load("data_clean.RData")

################################################################################
# FACTOR ANALYSIS AND RELIABILITY
################################################################################

unique_pretests <- data_clean %>%
  filter(test == "pre") %>%
  arrange(ID, year) %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup()

pre_efficacy_all <- unique_pretests %>%
  select(all_of(efficacy_items)) %>%
  mutate(across(all_of(reverse_efficacy_items), ~ 7 - .x)) %>%
  drop_na()

# Retain item 10 in the initial model to document its weak loading
efficacy_initial <- run_factor_checks(
  pre_efficacy_all,
  "Initial efficacy scale",
  factors = 1,
  run_parallel = TRUE
)

pre_efficacy_final <- pre_efficacy_all %>% select(-Efficacy_10)
efficacy_final <- run_factor_checks(
  pre_efficacy_final,
  "Final efficacy scale without item 10",
  factors = 1,
  run_parallel = FALSE
)

pre_understanding <- unique_pretests %>%
  select(all_of(understanding_items)) %>%
  drop_na()

understanding_final <- run_factor_checks(
  pre_understanding,
  "Self-reported climate understanding scale",
  factors = 1,
  run_parallel = TRUE
)

################################################################################
# COMPOSITE SCORES AND CANONICAL HLM INPUT
################################################################################

data_scored <- data_clean %>%
  mutate(across(all_of(reverse_efficacy_items), ~ 7 - .x)) %>%
  mutate(
    Efficacy_Score = row_mean_or_na(pick(all_of(efficacy_items_scored))),
    ccUnder_Score = row_mean_or_na(pick(all_of(understanding_items))),
    test = factor(test, levels = c("pre", "post")),
    year = factor(year),
    ID = factor(ID)
  )

# Save canonical RData input and a reviewable CSV copy
save(data_scored, file = "data_scored.RData")
write.csv(data_scored, "data_scored_full.csv", row.names = FALSE)

################################################################################
# COMPOSITE CHANGE SCORES AND VIOLIN FIGURE
################################################################################

paired_changes <- data_scored %>%
  select(ID, year, test, all_of(names(outcome_labels))) %>%
  pivot_longer(
    cols = all_of(names(outcome_labels)),
    names_to = "Outcome_Code",
    values_to = "Score"
  ) %>%
  pivot_wider(names_from = test, values_from = Score) %>%
  drop_na(pre, post) %>%
  mutate(
    Change = post - pre,
    Outcome = factor(
      unname(outcome_labels[Outcome_Code]),
      levels = unname(outcome_labels)
    ),
    year = factor(year, levels = sort(unique(as.character(year)))),
    year_index = as.numeric(year)
  )

change_summary <- summarize_change_scores(paired_changes)
print(change_summary)

# Plot raw cohort means and descriptive 95% confidence intervals
efficacy_violin <- make_violin_plot(
  paired_changes,
  change_summary,
  "Teaching efficacy",
  "Change in Teaching Efficacy",
  "Score Difference (Post - Pre)"
)

understanding_violin <- make_violin_plot(
  paired_changes,
  change_summary,
  "Self-reported climate understanding",
  "Change in Self-Reported Understanding"
)

violin_figure <- efficacy_violin + understanding_violin
print(violin_figure)

ggsave(
  "Figure_Violins.pdf",
  violin_figure,
  width = 8.5,
  height = 4.5,
  units = "in"
)

################################################################################
# TOPIC-LEVEL FAMILIARITY
################################################################################

topics_to_analyze <- c(
  "topic_climateCauses",
  "topic_modeling",
  "topic_misconceptions",
  "topic_evidence",
  "topic_effects",
  "topic_geoengineering",
  "topic_stories",
  "topic_inequity"
)

topic_pairs <- make_paired_items(data_clean, topics_to_analyze)

topic_response_labels <- c(
  "1" = "I am unfamiliar with this topic",
  "2" = "I know a little about this topic",
  "3" = "I know a moderate amount about this topic",
  "4" = "I know a lot about this topic"
)

topic_wilcoxon <- summarize_wilcoxon(topic_pairs)

topic_plot_list <- make_likert_plot_list(
  topic_pairs,
  response_levels = 1:4,
  response_labels = topic_response_labels,
  wilcoxon_results = topic_wilcoxon
)

topic_figure <- combine_likert_plots(
  topic_plot_list,
  ncol = 2
)

print(topic_figure)

ggsave(
  "Figure_Topics.pdf",
  topic_figure,
  width = 10,
  height = 12,
  units = "in"
)

topic_increase_table <- summarize_increases(topic_pairs)
print(topic_increase_table)
print(topic_wilcoxon)
write.csv(topic_wilcoxon, "Table_Wilcoxon_Topics.csv", row.names = FALSE)
