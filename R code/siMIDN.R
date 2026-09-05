library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)

# ============================
# DATA
# ============================
df_E <- data.frame(
  repeat_id = rep(1:3, times = 3),
  group = rep(c("Mock", "siEMC3", "siMIDN"), each = 3),
  value = c(
    100.00000000, 100.00000000, 100.00000000,
    114.9741079, 113.4627910, 77.8695014,
    176.8363352, 249.1691470, 143.2354370
  )
)

df_E$group <- factor(df_E$group, levels = c("Mock", "siEMC3", "siMIDN"))

summary_E <- df_E %>%
  group_by(group) %>%
  summarise(mean_val = mean(value), sem_val = sd(value)/sqrt(n()), .groups = "drop")

# ============================
# STATS — one-sample t-test vs Mock, log-transformed (same method as panel C)
# ============================
# NOTE: as in the panel C script, the "control" vector for each treatment is
# just rep(100, n), so a paired t-test against it is mathematically identical
# to a one-sample t-test against mu = log2(100). Written explicitly here so
# the test no longer silently depends on row ordering / implicit pairing.
treatment_groups_E <- levels(df_E$group)[-1]

stat_results_E <- lapply(treatment_groups_E, function(g) {
  treat_vals <- df_E$value[df_E$group == g]
  test <- t.test(log2(treat_vals) - log2(100))
  data.frame(group1 = "Mock", group2 = g, p_value = test$p.value)
}) %>% bind_rows()

stat_results_E$p_adj <- p.adjust(stat_results_E$p_value, method = "BH")

stat_results_E <- stat_results_E %>%
  mutate(
    p.adj.signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    group1 = factor(group1, levels = levels(df_E$group)),
    group2 = factor(group2, levels = levels(df_E$group))
  ) %>%
  arrange(group2)

print(stat_results_E)

# ============================
# COLOR PALETTE — matched to panel C
# ============================
color_map_E <- c("Mock" = "#999999", "siEMC3" = "#4A90D9", "siMIDN" = "#8E5FBF")

# ============================
# ANNOTATION — built dynamically from stat_results_E so it can never drift
# out of sync with the actual computed p-values
# ============================
annot_lines_E <- sprintf("Mock vs. %s, p.adj = %.3f",
                         stat_results_E$group2, stat_results_E$p_adj)
annot_text_E <- paste0(paste(annot_lines_E, collapse = "; "),
                       "\npaired t-test, log-transformed, BH-corrected")

# ============================
# PLOT — same styling as panel C
# ============================
p_E <- ggplot(summary_E, aes(x = group, y = mean_val, fill = group)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", linewidth = 1.1) +
  geom_errorbar(aes(ymin = mean_val - sem_val, ymax = mean_val + sem_val),
                width = 0.15, linewidth = 0.5) +
  geom_jitter(data = df_E, aes(x = group, y = value),
              width = 0.06, size = 2.4, shape = 24, fill = "black",
              color = "black", stroke = 0.4) +
  scale_fill_manual(values = color_map_E) +
  scale_x_discrete(expand = expansion(add = 0.4)) +
  labs(x = "Treatment",
       y = "EGR1 levels relative to loading control (%)",
       title = "Densitometric analysis of EGR1 levels") +
  annotate("text", x = 2, y = max(df_E$value) * 1.1,
           label = annot_text_E,
           size = 3.3, fontface = "italic") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    plot.margin = margin(12, 12, 25, 12)
  ) +
  coord_cartesian(ylim = c(0, max(df_E$value) * 1.25), clip = "off")

p_E
ggsave("panel_E.pdf", plot = p_E, width = 3.2, height = 4.5, units = "in")
ggsave("panel_E.tiff", plot = p_E, width = 3.2, height = 4.5, units = "in", dpi = 600)

# ============================
# BY-REPEAT PLOTS (QC / supplementary)
# ============================
make_repeat_plot_E <- function(rep_num) {
  sub_df <- df_E %>% filter(repeat_id == rep_num)
  
  ggplot(sub_df, aes(x = group, y = value, fill = group)) +
    geom_bar(stat = "identity", width = 0.5, color = "black", linewidth = 0.4) +
    geom_point(size = 2, shape = 24, fill = "black", color = "black", stroke = 0.4) +
    scale_fill_manual(values = color_map_E) +
    labs(x = NULL, y = "EGR1 levels (%)",
         title = paste("Repeat", rep_num)) +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 8),
      plot.title = element_text(size = 9, hjust = 0.5, face = "bold")
    )
}

repeat_plots_E <- lapply(1:3, make_repeat_plot_E)

p_repeats_E <- wrap_plots(repeat_plots_E, nrow = 1) &
  ylim(0, max(df_E$value) * 1.05)

p_repeats_E
ggsave("panel_E_by_repeat.pdf", plot = p_repeats_E, width = 7, height = 3, units = "in")