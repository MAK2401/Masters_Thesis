library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)

# ============================
# 1. DATA
# ============================
df <- data.frame(
  repeat_id = rep(1:4, times = 5),
  group = rep(c("EGR1 + EV",
                "EGR1 + 250ng MIDN-S",
                "EGR1 + 125ng MIDN-S",
                "EGR1 + 62.5ng MIDN-S",
                "EGR1 + 31.25ng MIDN-S"), each = 4),
  value = c(
    100, 100, 100, 100,
    18.21322001, 93.67833511, 14.21513071, 30.50353382,
    12.65383781, 241.0926024, 68.48392857, 51.54543745,
    42.87661452, 212.3236262, 97.75958405, 71.76279957,
    53.80274795, 258.8763104, 224.2190311, 58.1252958
  )
)

df$group <- factor(df$group,
                   levels = c("EGR1 + EV",
                              "EGR1 + 250ng MIDN-S",
                              "EGR1 + 125ng MIDN-S",
                              "EGR1 + 62.5ng MIDN-S",
                              "EGR1 + 31.25ng MIDN-S"))

# ============================
# 2. SUMMARY STATS
# ============================
summary_df <- df %>%
  group_by(group) %>%
  summarise(mean_val = mean(value), sem_val = sd(value)/sqrt(n()), .groups = "drop")

# ============================
# 3. STATISTICS — one-sample t-test vs EV (log-transformed)
# ============================
# NOTE: the "control" vector for each treatment is just rep(100, n), so a
# paired t-test against it is mathematically identical to a one-sample
# t-test against mu = log2(100). Written explicitly here so the test no
# longer silently depends on row ordering / implicit pairing by position.
treatment_groups <- levels(df$group)[-1]

stat_results <- lapply(treatment_groups, function(g) {
  treat_vals <- df$value[df$group == g]
  test <- t.test(log2(treat_vals) - log2(100))
  data.frame(group1 = "EGR1 + EV", group2 = g, p_value = test$p.value)
}) %>% bind_rows()

stat_results$p_adj <- p.adjust(stat_results$p_value, method = "BH")

stat_results <- stat_results %>%
  mutate(
    p.adj.signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    group1 = factor(group1, levels = levels(df$group)),
    group2 = factor(group2, levels = levels(df$group))
  ) %>%
  arrange(group2)

print(stat_results)

max_y <- max(df$value)
stat_results$y.position <- max_y + seq(15, by = 20, length.out = nrow(stat_results))

# ============================
# 4. COLOR PALETTE (matched to reference style)
# ============================
color_map <- c("EGR1 + EV" = "#999999",
               "EGR1 + 250ng MIDN-S" = "#D64545",
               "EGR1 + 125ng MIDN-S" = "#E8934A",
               "EGR1 + 62.5ng MIDN-S" = "#E8C547",
               "EGR1 + 31.25ng MIDN-S" = "#5DC9AA")

# ============================
# 5. COMBINED PLOT
# ============================
p_combined <- ggplot(summary_df, aes(x = group, y = mean_val, fill = group)) +
  geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean_val - sem_val, ymax = mean_val + sem_val),
                width = 0.2, linewidth = 0.5) +
  geom_jitter(data = df, aes(x = group, y = value),
              width = 0.08, size = 2.2, shape = 24, fill = "black",
              color = "black", stroke = 0.4) +
  stat_pvalue_manual(stat_results, label = "p.adj.signif",
                     tip.length = 0.01, bracket.size = 0.4, hide.ns = FALSE) +
  scale_fill_manual(values = color_map) +
  labs(x = "Treatment",
       y = "Relative EGR1 protein expression (%)\n(normalised to loading control)",
       title = "Densitometric analysis of EGR1 levels") +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 8, angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    plot.margin = margin(10, 10, 20, 10)
  )

p_combined
ggsave("densitometry_combined.pdf", plot = p_combined, width = 4.5, height = 4, units = "in")
ggsave("densitometry_combined.tiff", plot = p_combined, width = 4.5, height = 4, units = "in", dpi = 600)

# ============================
# 6. BY-REPEAT PLOTS (QC / supplementary)
# ============================
make_repeat_plot <- function(rep_num) {
  sub_df <- df %>% filter(repeat_id == rep_num)
  
  ggplot(sub_df, aes(x = group, y = value, fill = group)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.4) +
    geom_point(size = 2, shape = 24, fill = "black", color = "black", stroke = 0.4) +
    scale_fill_manual(values = color_map) +
    labs(x = NULL, y = "Relative EGR1\nexpression (%)",
         title = paste("Repeat", rep_num)) +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
      plot.title = element_text(size = 9, hjust = 0.5, face = "bold")
    )
}

repeat_plots <- lapply(1:4, make_repeat_plot)

p_repeats_combined <- wrap_plots(repeat_plots, nrow = 1) &
  ylim(0, max(df$value) * 1.05)

p_repeats_combined
ggsave("densitometry_by_repeat.pdf", plot = p_repeats_combined,
       width = 10, height = 3, units = "in")

# ============================
# PANEL C — full-size figure with dynamic stats annotation
# ============================
# Recompute the full stats table first (needed for the annotation string
# below) so Panel C can never drift out of sync with the actual test results.
stat_results_check <- lapply(treatment_groups, function(g) {
  treat_vals <- df$value[df$group == g]
  test <- t.test(log2(treat_vals) - log2(100))
  data.frame(
    group = g,
    t_stat = round(test$statistic, 3),
    df_denom = test$parameter,
    p_value = test$p.value,
    mean_diff_log2 = test$estimate
  )
}) %>% bind_rows()

stat_results_check$p_adj <- p.adjust(stat_results_check$p_value, method = "BH")

print(stat_results_check)

# Build the annotation text directly from the computed (BH-corrected) p-values
# so the label on the plot always matches what was actually calculated.
annot_lines <- sprintf("%s: p.adj = %.3f",
                       gsub("EGR1 \\+ ", "", stat_results_check$group),
                       stat_results_check$p_adj)
annot_text <- paste0(paste(annot_lines, collapse = "; "),
                     "\npaired t-test vs. EV, BH-corrected")

p_C <- ggplot(summary_df, aes(x = group, y = mean_val, fill = group)) +
  geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 1.1) +
  geom_errorbar(aes(ymin = mean_val - sem_val, ymax = mean_val + sem_val),
                width = 0.2, linewidth = 0.5) +
  geom_jitter(data = df, aes(x = group, y = value),
              width = 0.08, size = 2.4, shape = 24, fill = "black",
              color = "black", stroke = 0.4) +
  scale_fill_manual(values = color_map) +
  labs(x = "Treatment",
       y = "EGR1 levels relative to loading control (%)",
       title = "Densitometric analysis of EGR1 levels") +
  annotate("text", x = 3, y = max(df$value) * 1.1,
           label = annot_text,
           size = 3.3, fontface = "italic") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, face = "bold", angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    plot.margin = margin(12, 12, 25, 12)
  ) +
  coord_cartesian(ylim = c(0, max(df$value) * 1.2), clip = "off")

p_C
ggsave("panel_C.pdf", plot = p_C, width = 5, height = 4.5, units = "in")
ggsave("panel_C.tiff", plot = p_C, width = 5, height = 4.5, units = "in", dpi = 600)