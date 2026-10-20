# ============================================================================
# qPCR RELATIVE ABUNDANCE ANALYSIS - ermF GENE
# ============================================================================
# ermF (erythromycin resistance) gene analysis
# Using RAW SQ (Standard Quantification) copy numbers with LOG SCALE
# Analysis of 28 patients with 3 replicates each
# Created: October 2026
# Version: Adapted from MASTER TEMPLATE
# ============================================================================

library(readxl)
library(tidyverse)
library(ggplot2)
library(openxlsx)

# ============================================================================
# CONFIGURATION - ERMF GENE
# ============================================================================
GENE_NAME <- "ermF"  # Erythromycin resistance gene
SHEET_NAME <- "emF"  # Sheet name in Excel file
FILE_PATH <- "/Users/macbookair/Desktop/qpcr/qPCR_normalized_by_16S_NEW_updated.xlsx"
OUTPUT_FOLDER <- "/Users/macbookair/Desktop/qpcr/"

# ============================================================================
# LOAD RAW DATA (SQ values - NOT normalized)
# ============================================================================
gene_raw <- read_excel(FILE_PATH, sheet = SHEET_NAME)

# ============================================================================
# PREPARE LONG FORMAT FOR PLOTTING
# ============================================================================
gene_long <- gene_raw |>
  pivot_longer(
    cols = c(Rep1, Rep2, Rep3),
    names_to = "Replicate",
    values_to = "copies"
  ) |>
  mutate(
    Patient = as.numeric(gsub("P", "", Sample))
  )

# ============================================================================
# SUMMARY STATISTICS WITH 95% CONFIDENCE INTERVALS
# ============================================================================
gene_summary <- gene_long |>
  group_by(Sample, Patient) |>
  summarise(
    mean_copies = mean(copies, na.rm = TRUE),
    sd_copies = sd(copies, na.rm = TRUE),
    n = n(),
    se_copies = sd_copies / sqrt(n),
    ci_lower = mean_copies - (qt(0.975, n-1) * se_copies),
    ci_upper = mean_copies + (qt(0.975, n-1) * se_copies),
    min_copies = min(copies, na.rm = TRUE),
    max_copies = max(copies, na.rm = TRUE),
    cv_pct = (sd_copies / mean_copies) * 100,
    .groups = "drop"
  ) |>
  arrange(mean_copies) |>
  mutate(Sample = fct_reorder(Sample, mean_copies))

# ============================================================================
# CREATE PLOT - LOG SCALE WITH SUPERSCRIPT NOTATION
# ============================================================================
p_gene <- ggplot(
  gene_long |> mutate(Sample = fct_reorder(Sample, copies, .fun = mean)), 
  aes(x = Sample, y = copies)
) +
  # 95% Confidence Interval band
  geom_ribbon(
    data = gene_summary |> mutate(SampleOrd = fct_reorder(Sample, mean_copies)),
    aes(x = SampleOrd, y = mean_copies, ymin = ci_lower, ymax = ci_upper),
    fill = "gray80",
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  # Individual replicate points
  geom_point(
    aes(color = Replicate), 
    size = 3.5, 
    alpha = 0.7,
    position = position_jitter(width = 0.12, height = 0)
  ) +
  # Min-Max error bars
  geom_errorbar(
    data = gene_summary |> mutate(SampleOrd = fct_reorder(Sample, mean_copies)),
    aes(x = SampleOrd, y = mean_copies, ymin = min_copies, ymax = max_copies),
    color = "gray40",
    width = 0.12,
    linewidth = 0.5,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  # Mean points - grey diamond
  geom_point(
    data = gene_summary |> mutate(SampleOrd = fct_reorder(Sample, mean_copies)),
    aes(x = SampleOrd, y = mean_copies),
    color = "gray50",
    size = 2.5,
    shape = 23,
    fill = "gray70",
    stroke = 0.6,
    inherit.aes = FALSE
  ) +
  # Log scale with superscript notation
  scale_y_log10(
    limits = c(1, 1e6),
    breaks = 10^(1:6),
    labels = function(x) {
      exponents <- log10(x)
      parse(text = paste0("10^", exponents))
    }
  ) +
  # Colors for replicates
  scale_color_manual(
    values = c("Rep1" = "#20B2AA", "Rep2" = "#9370DB", "Rep3" = "#20B2AA"),
    breaks = c("Rep1", "Rep2", "Rep3"),
    labels = rep("Measured replicate", 3)
  ) +
  ggtitle(paste("Relative abundance of", GENE_NAME)) +
  labs(
    subtitle = "qPCR quantification (by standard curve)",
    x = "Patient",
    y = paste(GENE_NAME, "abundance (% of 16S copies)"),
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
    axis.title.y = element_text(hjust = 0.5, size = 11, face = "bold"),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# ============================================================================
# STATISTICAL ANALYSIS & INTERPRETATION
# ============================================================================

# Overall statistics
overall_stats <- list(
  Gene = GENE_NAME,
  N_Patients = nrow(gene_summary),
  N_Replicates_Total = nrow(gene_long),
  Mean_Copies = round(mean(gene_summary$mean_copies), 2),
  Median_Copies = round(median(gene_summary$mean_copies), 2),
  SD_Copies = round(sd(gene_summary$mean_copies), 2),
  Min_Copies = round(min(gene_summary$mean_copies), 2),
  Max_Copies = round(max(gene_summary$mean_copies), 2),
  Range = round(max(gene_summary$mean_copies) - min(gene_summary$mean_copies), 2),
  CV_Overall = round((sd(gene_summary$mean_copies) / mean(gene_summary$mean_copies)) * 100, 2)
)

# High expressors (top 25%)
high_threshold <- quantile(gene_summary$mean_copies, 0.75)
high_expressors <- gene_summary |>
  filter(mean_copies >= high_threshold) |>
  select(Sample, mean_copies, ci_lower, ci_upper, cv_pct) |>
  arrange(desc(mean_copies))

# Low expressors (bottom 25%)
low_threshold <- quantile(gene_summary$mean_copies, 0.25)
low_expressors <- gene_summary |>
  filter(mean_copies <= low_threshold) |>
  select(Sample, mean_copies, ci_lower, ci_upper, cv_pct) |>
  arrange(mean_copies)

# ============================================================================
# CREATE EXCEL REPORT
# ============================================================================

wb <- createWorkbook()

# Sheet 1: Summary Statistics
addWorksheet(wb, "Summary_Statistics")
writeData(wb, "Summary_Statistics", 
          data.frame(
            Metric = names(overall_stats),
            Value = as.character(unlist(overall_stats))
          ))

# Sheet 2: Patient Details
addWorksheet(wb, "Patient_Details")
patient_data <- gene_summary |>
  select(Sample, Patient, mean_copies, sd_copies, se_copies, ci_lower, ci_upper, n, cv_pct) |>
  mutate(across(where(is.numeric), ~round(., 2)))
writeData(wb, "Patient_Details", patient_data)

# Sheet 3: High Expressors (Top 25%)
addWorksheet(wb, "High_Expressors_Top25")
high_exp_display <- high_expressors |>
  mutate(across(where(is.numeric), ~round(., 2)))
writeData(wb, "High_Expressors_Top25", high_exp_display)

# Sheet 4: Low Expressors (Bottom 25%)
addWorksheet(wb, "Low_Expressors_Bottom25")
low_exp_display <- low_expressors |>
  mutate(across(where(is.numeric), ~round(., 2)))
writeData(wb, "Low_Expressors_Bottom25", low_exp_display)

# Sheet 5: Interpretation & Findings
addWorksheet(wb, "Interpretation")
interpretation <- data.frame(
  Finding = c(
    "Gene Name",
    "Number of Patients",
    "Mean Copy Number",
    "Median Copy Number",
    "Range",
    "Coefficient of Variation",
    "High Expressors (>75th percentile)",
    "Low Expressors (<25th percentile)",
    "Variability Assessment",
    "Key Observation"
  ),
  Details = c(
    GENE_NAME,
    as.character(overall_stats$N_Patients),
    paste0(round(overall_stats$Mean_Copies, 2), " copies"),
    paste0(round(overall_stats$Median_Copies, 2), " copies"),
    paste0(round(overall_stats$Range, 2), " copies"),
    paste0(round(overall_stats$CV_Overall, 2), "%"),
    paste0(nrow(high_expressors), " patients (", 
           paste0(unique(high_expressors$Sample), collapse = ", "), ")"),
    paste0(nrow(low_expressors), " patients (", 
           paste0(unique(low_expressors$Sample), collapse = ", "), ")"),
    if_else(overall_stats$CV_Overall > 100, "Very high variability (CV > 100%)", 
            if_else(overall_stats$CV_Overall > 50, "High variability (50% < CV < 100%)", "Moderate variability")),
    paste0("The ", GENE_NAME, " gene shows ", 
           if_else(overall_stats$Range > 10000, "dramatic", "moderate"),
           " variation across patients with mean of ", 
           round(overall_stats$Mean_Copies, 0), " copies and median of ",
           round(overall_stats$Median_Copies, 0), " copies.")
  )
)
writeData(wb, "Interpretation", interpretation)

# Save Excel file
excel_filename <- paste0(OUTPUT_FOLDER, GENE_NAME, "_qPCR_analysis_", 
                         format(Sys.Date(), "%Y%m%d"), ".xlsx")
saveWorkbook(wb, excel_filename, overwrite = TRUE)

cat("✓ Excel report saved to:", excel_filename, "\n")

# ============================================================================
# SAVE PLOT
# ============================================================================

plot_filename <- paste0(OUTPUT_FOLDER, GENE_NAME, "_relative_abundance_plot_", 
                        format(Sys.Date(), "%Y%m%d"), ".png")
ggsave(plot_filename, p_gene, width = 14, height = 8, dpi = 300, bg = "white")

cat("✓ Plot saved to:", plot_filename, "\n")

# ============================================================================
# DISPLAY RESULTS
# ============================================================================

cat("\n")
cat(paste(rep("=", 75), collapse = ""), "\n")
cat(toupper(paste("ANALYSIS COMPLETE:", GENE_NAME)), "\n")
cat(paste(rep("=", 75), collapse = ""), "\n")
cat("\nSummary Statistics:\n")
print(data.frame(Metric = names(overall_stats), Value = unlist(overall_stats)))
cat("\nHigh Expressors (Top 25%):\n")
print(high_expressors)
cat("\nLow Expressors (Bottom 25%):\n")
print(low_expressors)
cat("\n✓ Files saved to:", OUTPUT_FOLDER, "\n\n")

# ============================================================================
# DISPLAY PLOT IN RSTUDIO PLOTS PANE
# ============================================================================
print(p_gene)
