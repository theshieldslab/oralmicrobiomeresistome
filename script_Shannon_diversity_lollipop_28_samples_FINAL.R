#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# SHANNON DIVERSITY LOLLIPOP — 28 SAMPLES (CORRECTED GENUS DATA)
# Lollipop chart with Shannon values displayed on each point
# Colored by caries risk group
# Uses CORRECTED Genus_Abundance sheet (82 truly aggregated genera)
# Saves PNG + PDF plots + Excel interpretation with 5 sheets
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(readxl)
library(openxlsx2)
library(ggplot2)

# ── FILE PATHS ─────────────────────────────────────────────────────────────────
main_folder <- "/Users/macbookair/Desktop/16s may"
master_path <- file.path(main_folder, "16S_Master_Summary_UPDATED_28_samples_FINAL.xlsx")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Plot files (no timestamp)
output_lollipop_png <- file.path(main_folder, 
  "Shannon_diversity_lollipop_28_samples_FINAL.png")
output_lollipop_pdf <- file.path(main_folder, 
  "Shannon_diversity_lollipop_28_samples_FINAL.pdf")

# Interpretation Excel (with timestamp, starts with "interpretation_")
output_excel <- file.path(main_folder, 
  paste0("interpretation_Shannon_diversity_lollipop_28_samples_FINAL_", timestamp, ".xlsx"))

# ── COLOR PALETTE ──────────────────────────────────────────────────────────────
col_palette <- c(low = "#4E84C4", high = "#C9503A")

# ── LOAD DATA ──────────────────────────────────────────────────────────────────
# Load the CORRECTED Genus_Abundance sheet (82 truly aggregated genera)
genus_data <- read_excel(master_path, sheet = "Genus_Abundance")

# Load the caries risk metadata
caries_raw <- read_excel(master_path, sheet = "caries risk")

cat("════════════════════════════════════════════════════════════════════════════\n")
cat("SHANNON DIVERSITY ANALYSIS — 28 SAMPLES (CORRECTED DATA)\n")
cat("════════════════════════════════════════════════════════════════════════════\n")
cat("Genus_Abundance loaded:", nrow(genus_data), "unique genera\n\n")

# ── COMPUTE SHANNON DIVERSITY ──────────────────────────────────────────────────
# Transpose: patients as rows, genera as columns
genus_matrix <- as.matrix(genus_data[, -1])
rownames(genus_matrix) <- genus_data$Genus

# Compute Shannon index for each patient (column)
shannon_vals <- numeric(ncol(genus_matrix))
patient_names <- colnames(genus_matrix)

for (i in 1:ncol(genus_matrix)) {
  abund <- genus_matrix[, i]
  abund <- abund[abund > 0]
  shannon_vals[i] <- -sum(abund * log(abund), na.rm = TRUE)
}

# Prepare metadata with caries risk
# Use only as many caries values as we have patients
num_patients <- ncol(genus_data) - 1
caries_clean <- as.character(unlist(caries_raw[[1]][1:num_patients]))
caries_clean <- tolower(trimws(caries_clean))

shannon_df <- data.frame(
  Sample = patient_names,
  Shannon = shannon_vals,
  caries_risk = factor(caries_clean, levels = c("low", "high")),
  stringsAsFactors = FALSE
)

# ── WILCOXON TEST ──────────────────────────────────────────────────────────────
wt <- wilcox.test(Shannon ~ caries_risk, data = shannon_df, exact = FALSE)
p_fmt <- ifelse(wt$p.value < 0.001, "p < 0.001",
                paste0("p = ", format(round(wt$p.value, 3), nsmall = 3)))

# ── SUMMARY STATISTICS ─────────────────────────────────────────────────────────
shannon_summary <- shannon_df %>%
  group_by(caries_risk) %>%
  summarise(
    n = n(),
    Mean = round(mean(Shannon), 3),
    Median = round(median(Shannon), 3),
    SD = round(sd(Shannon), 3),
    Min = round(min(Shannon), 3),
    Max = round(max(Shannon), 3),
    .groups = "drop"
  )

# ── FIND MIN/MAX PATIENTS ──────────────────────────────────────────────────────
shannon_sorted <- shannon_df[order(shannon_df$Shannon), ]
min_patient <- shannon_sorted$Sample[1]
min_value <- shannon_sorted$Shannon[1]
max_patient <- shannon_sorted$Sample[nrow(shannon_sorted)]
max_value <- shannon_sorted$Shannon[nrow(shannon_sorted)]

cat("Wilcoxon test:", p_fmt, "\n\n")
print(shannon_summary)
cat("\n")

# ── PREPARE LOLLIPOP DATA (sorted by Shannon) ──────────────────────────────────
lollipop_df <- shannon_df[order(shannon_df$Shannon), ]
lollipop_df$Sample <- factor(lollipop_df$Sample, levels = lollipop_df$Sample)

# ── CREATE LOLLIPOP (with Shannon values on points) ───────────────────────────
p_lollipop <- ggplot(lollipop_df, aes(x = Sample, y = Shannon)) +
  geom_segment(
    aes(xend = Sample, y = 0, yend = Shannon, color = caries_risk),
    linewidth = 1.2, alpha = 0.85
  ) +
  geom_point(aes(color = caries_risk), size = 5, alpha = 0.95) +
  # Values displayed on top of points
  geom_text(
    aes(label = round(Shannon, 2), color = caries_risk),
    vjust = -1.0, size = 3.2, fontface = "bold", show.legend = FALSE
  ) +
  scale_color_manual(values = col_palette, labels = c("Low", "High"), name = "Caries risk") +
  scale_y_continuous(limits = c(0, max(lollipop_df$Shannon) + 0.6), expand = c(0, 0)) +
  theme_bw(base_size = 13) +
  labs(
    title = "Shannon Diversity Index - 28 Patients",
    subtitle = paste("Wilcoxon test:", p_fmt),
    x = "Patient (sorted by Shannon index)",
    y = "Shannon index"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "darkgrey"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# ──────────────────────────────────────────────────────────────────────────────
# SAVE PNG & PDF PLOTS
# ──────────────────────────────────────────────────────────────────────────────
ggsave(output_lollipop_png, plot = p_lollipop, width = 14, height = 8, 
       dpi = 600, bg = "white")
ggsave(output_lollipop_pdf, plot = p_lollipop, width = 14, height = 8, 
       bg = "white")
cat("Plots saved:\n")
cat("  PNG:", output_lollipop_png, "\n")
cat("  PDF:", output_lollipop_pdf, "\n")

# ──────────────────────────────────────────────────────────────────────────────
# PREPARE EXCEL DATA SHEETS
# ──────────────────────────────────────────────────────────────────────────────
# Sheet 1: All Shannon data (sorted by Shannon index)
shannon_export <- shannon_df[order(shannon_df$Shannon), ]
colnames(shannon_export) <- c("Patient", "Shannon_Index", "Caries_Risk")

# Sheet 2: Summary statistics (overall)
lollipop_summary <- tibble(
  Metric = c("Total samples", "Mean Shannon", "Median Shannon", "SD",
             "Minimum Shannon", "Maximum Shannon"),
  Value = c(
    as.character(nrow(shannon_df)),
    as.character(round(mean(shannon_df$Shannon), 3)),
    as.character(round(median(shannon_df$Shannon), 3)),
    as.character(round(sd(shannon_df$Shannon), 3)),
    paste0(round(min_value, 3), " (", min_patient, ")"),
    paste0(round(max_value, 3), " (", max_patient, ")")
  )
)

# Sheet 3: Statistical test details
stat_summary <- tibble(
  Metric = c("Test used", "Test statistic", "p-value", "Interpretation"),
  Value = c(
    "Wilcoxon rank-sum test (low vs high caries risk)",
    paste0("W = ", round(wt$statistic, 2)),
    format(round(wt$p.value, 4), nsmall = 4),
    ifelse(wt$p.value < 0.05, 
           "Significant difference between groups", 
           "No significant difference detected")
  )
)

# Sheet 4: Summary statistics by group
group_summary <- shannon_summary

# Sheet 5: Interpretation text (human-readable)
overall_text <- paste0(
  "OVERALL INTERPRETATION\n\n",
  "Test: The mean Shannon diversity index across the 28 samples was ",
  round(mean(shannon_df$Shannon), 2), ", with a median of ",
  round(median(shannon_df$Shannon), 2), ".\n\n",
  "Lowest diversity: The lowest Shannon diversity was observed in ",
  min_patient, " with a Shannon index of ", round(min_value, 2), ".\n\n",
  "Highest diversity: The highest Shannon diversity was observed in ",
  max_patient, " with a Shannon index of ", round(max_value, 2), ".\n\n",
  "Detailed interpretation:\n",
  "Shannon diversity analysis at the genus level revealed variability in alpha diversity among ",
  "supragingival biofilm samples. Across the 28-sample dataset, the Shannon diversity index ranged from ",
  round(min_value, 2), " in ", min_patient, " to ", round(max_value, 2), " in ", max_patient, 
  ", with a mean value of ", round(mean(shannon_df$Shannon), 2), ". These findings indicate ",
  "inter-individual differences in genus-level community richness and evenness.\n\n",
  "STATISTICAL TEST RESULTS\n\n",
  "Test: Wilcoxon rank-sum test\n",
  "Statistic: W = ", round(wt$statistic, 2), "\n",
  "p-value: ", format(round(wt$p.value, 4), nsmall = 4), "\n",
  "Result: ", ifelse(wt$p.value < 0.05, 
                     "There is a statistically significant difference in Shannon diversity between low and high caries risk groups.", 
                     "There is no statistically significant difference in Shannon diversity between low and high caries risk groups.")
)
interpretation_df <- data.frame(Text = overall_text)

# ──────────────────────────────────────────────────────────────────────────────
# SAVE EXCEL WITH INTERPRETATION SHEETS
# ──────────────────────────────────────────────────────────────────────────────
wb <- openxlsx2::wb_workbook()
# Sheet 1: All Shannon data (sorted)
wb$add_worksheet("Shannon_All_Data_Sorted")
wb$add_data(x = shannon_export)
# Sheet 2: Summary statistics (overall)
wb$add_worksheet("Summary_Statistics")
wb$add_data(x = lollipop_summary)
# Sheet 3: Statistical test details
wb$add_worksheet("Statistical_Test")
wb$add_data(x = stat_summary)
# Sheet 4: Summary by group (low vs high caries risk)
wb$add_worksheet("Summary_By_Group")
wb$add_data(x = group_summary)
# Sheet 5: Full interpretation text
wb$add_worksheet("Interpretation")
wb$add_data(x = interpretation_df)
wb$save(output_excel)

cat("\nInterpretation Excel file saved:\n")
cat("  ", output_excel, "\n\n")
cat("Sheets included:\n")
cat("  1. Shannon_All_Data_Sorted (sorted by Shannon index)\n")
cat("  2. Summary_Statistics (overall)\n")
cat("  3. Statistical_Test (Wilcoxon test results)\n")
cat("  4. Summary_By_Group (low vs high caries risk)\n")
cat("  5. Interpretation (full text interpretation)\n")
cat("\n════════════════════════════════════════════════════════════════════════════\n")
cat("SHANNON DIVERSITY ANALYSIS COMPLETE!\n")
cat("════════════════════════════════════════════════════════════════════════════\n")
