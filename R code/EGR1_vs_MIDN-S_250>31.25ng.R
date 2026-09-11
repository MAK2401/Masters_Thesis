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
# 2. COLOR PALETTE (matched to reference style)
# ============================
color_map <- c("EGR1 + EV" = "#999999",
               "EGR1 + 250ng MIDN-S" = "#D64545",
               "EGR1 + 125ng MIDN-S" = "#E8934A",
               "EGR1 + 62.5ng MIDN-S" = "#E8C547",
               "EGR1 + 31.25ng MIDN-S" = "#5DC9AA")


# ============================
# 3. BY-REPEAT PLOTS (QC / supplementary)
# ============================
make_repeat_plot <- function(rep_num) {
  sub_df <- df %>% filter(repeat_id == rep_num)
  
  ggplot(sub_df, aes(x = group, y = value, fill = group)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.4) +
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
