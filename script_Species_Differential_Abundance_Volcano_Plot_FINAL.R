#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# STANDALONE SCRIPT: SPECIES DIFFERENTIAL ABUNDANCE VOLCANO PLOT
# 
# Purpose: Visualize species significantly associated with LOW vs HIGH caries risk
#          using volcano plot with publication-quality formatting
#
# Output: Publication-ready volcano plot (PNG + PDF) showing:
#         - X-axis: log2(Fold Change)
#         - Y-axis: -log10(p-value)
#         - Light purple: Non-significant (p >= 0.05)
#         - Dark purple: Significant (p < 0.05)
#         - Species labels: OUTSIDE plot area with arrows
#
# Author: Analysis Pipeline
# Date: September 2026
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(ggrepel)
library(tidyverse)

cat("\n╔══════════════════════════════════════════════════════════════════╗\n")
cat("║  SPECIES DIFFERENTIAL ABUNDANCE BY CARIES RISK - VOLCANO PLOT   ║\n")
cat("║  Publication-Quality Visualization                             ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# Get the volcano plot data
volcano_data <- p_volcano_final$data

cat("Creating publication-quality volcano plot...\n\n")

# Create volcano plot
p_volcano_publication <- ggplot(volcano_data, aes(x = log2FC, y = neg_log10_pval)) +
  
  geom_point(
    aes(color = Significant, size = Significant),
    alpha = 0.8,
    stroke = 0.4,
    show.legend = TRUE
  ) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999", linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999", linewidth = 0.6) +
  
  scale_color_manual(
    values = c("Not significant" = "#d0c9e8", "Significant" = "#6b5b95"),
    name = ""
  ) +
  scale_size_manual(
    values = c("Not significant" = 2.5, "Significant" = 4),
    guide = "none"
  ) +
  
  geom_text_repel(
    data = volcano_data %>% filter(Significant == "Significant"),
    aes(label = Species),
    size = 4,
    fontface = "italic",
    box.padding = 1.0,
    point.padding = 0.5,
    arrow = arrow(length = unit(0.02, 'npc'), type = 'closed', ends = 'last'),
    force = 3,
    show.legend = FALSE,
    xlim = c(NA, NA),
    ylim = c(NA, NA)
  ) +
  
  labs(
    x = "log2(Fold Change)",
    y = "-log10(p-value)",
    title = "Species Differential Abundance by Caries Risk",
    caption = "← more abundant in LOW risk | more abundant in HIGH risk →"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 11, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    plot.caption = element_text(size = 11, face = "plain", hjust = 0.5, margin = margin(t = 10)),
    aspect.ratio = 0.65,
    plot.margin = margin(20, 20, 15, 15)
  )

print(p_volcano_publication)

# Save PNG
ggsave(
  "~/Desktop/16s may/Species_Differential_Abundance_Volcano_Plot_FINAL.png",
  plot = p_volcano_publication,
  width = 14,
  height = 9,
  dpi = 300,
  bg = "white"
)

# Save PDF
ggsave(
  "~/Desktop/16s may/Species_Differential_Abundance_Volcano_Plot_FINAL.pdf",
  plot = p_volcano_publication,
  width = 14,
  height = 9,
  dpi = 300,
  bg = "white"
)

cat("\n════════════════════════════════════════════════════════════════════\n")
cat("✅ VOLCANO PLOT SAVED\n")
cat("📊 PNG:  Species_Differential_Abundance_Volcano_Plot_FINAL.png\n")
cat("📄 PDF:  Species_Differential_Abundance_Volcano_Plot_FINAL.pdf\n")
cat("════════════════════════════════════════════════════════════════════\n")
