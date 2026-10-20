# ============================================================================
# qPCR RELATIVE ABUNDANCE ANALYSIS - ermB
# ============================================================================
# Antimicrobial resistance gene: ermB (erythromycin resistance)
# Analysis of 28 patients with 3 replicates each
# Updated: August 2026
# Version: FINAL - UPDATED with Patient 22 correction
# ============================================================================

library(readxl)
library(tidyverse)
library(ggplot2)
library(openxlsx)

# ============================================================================
# CONFIGURATION
# ============================================================================
GENE_NAME <- "ermB"
FILE_PATH <- "/Users/macbookair/Desktop/qpcr/qPCR_first7_.xlsx"
OUTPUT_FOLDER <- "/Users/macbookair/Desktop/qpcr/"

# ============================================================================
# LOAD DATA
# ============================================================================
gene_raw <- read_excel(FILE_PATH, sheet = GENE_NAME)
s16_raw <- read_excel(FILE_PATH, sheet = "16s")

# Get 16S mean values
s16_mean_df <- s16_raw |> 
  select(Sample, starts_with("Mean")) |>
  set_names(c("Sample", "s16_mean"))

# ============================================================================
# CALCULATE PERCENTAGES (normalized by 16S mean)
# ============================================================================
gene_percentages <- gene_raw |>
  left_join(s16_mean_df, by = "Sample") |>
  mutate(
    Rep1_pct = (Rep1 / s16_mean) * 100,
    Rep2_pct = (Rep2 / s16_mean) * 100,
    Rep3_pct = (Rep3 / s16_mean) * 100
  ) |>
  select(Sample, Rep1_pct, Rep2_pct, Rep3_pct)

# ============================================================================
# PREPARE LONG FORMAT FOR PLOTTING
# ============================================================================
gene_pct_long <- gene_percentages |>
  pivot_longer(
    cols = c(Rep1_pct, Rep2_pct, Rep3_pct),
    names_to = "Replicate",
    values_to = "percentage"
  ) |>
  mutate(
    Replicate = str_replace(Replicate, "_pct", ""),
    Patient = as.numeric(str_extract(Sample, "\\d+"))
  )

# ============================================================================
# SUMMARY STATISTICS WITH 95% CONFIDENCE INTERVALS
# ============================================================================
gene_summary <- gene_pct_long |>
  group_by(Sample, Patient) |>
  summarise(
    mean_pct = mean(percentage, na.rm = TRUE),
    sd_pct = sd(percentage, na.rm = TRUE),
    n = n(),
    se_pct = sd_pct / sqrt(n),
    ci_lower = mean_pct - (qt(0.975, n-1) * se_pct),
    ci_upper = mean_pct + (qt(0.975, n-1) * se_pct),
    min_pct = min(percentage, na.rm = TRUE),
    max_pct = max(percentage, na.rm = TRUE),
    cv_pct = (sd_pct / mean_pct) * 100,  # Coefficient of variation
    .groups = "drop"
  ) |>
  arrange(mean_pct) |>
  mutate(Sample = fct_reorder(Sample, mean_pct))

# ============================================================================
# CREATE PLOT - FINAL TEMPLATE
# ============================================================================
p_gene <- ggplot(
  gene_pct_long |> mutate(Sample = fct_reorder(Sample, percentage, .fun = mean)), 
  aes(x = Sample, y = percentage)
) +
  # 95% Confidence Interval band
  geom_ribbon(
    data = gene_summary |> mutate(SampleOrd = fct_reorder(Sample, mean_pct)),
    aes(x = SampleOrd, y = mean_pct, ymin = ci_lower, ymax = ci_upper),
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
    data = gene_summary |> mutate(SampleOrd = fct_reorder(Sample, mean_pct)),
    aes(x = SampleOrd, y = mean_pct, ymin = min_pct, ymax = max_pct),
    color = "gray40",
    width = 0.12,
    linewidth = 0.5,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  # Mean points - grey diamond
  geom_point(
    data = gene_summary |> mutate(SampleOrd = fct_reorder(Sample, mean_pct)),
    aes(x = SampleOrd, y = mean_pct),
    color = "gray50",
    size = 2.5,
    shape = 23,
    fill = "gray70",
    stroke = 0.6,
    inherit.aes = FALSE
  ) +
  # Log scale with superscript notation - showing -4 to 8
  scale_y_log10(
    limits = c(1e-4, 1e8),
    breaks = c(1e-4, 1e-1, 1e2, 1e5, 1e8),  # Only these 5 positions
    minor_breaks = NULL,  # No minor breaks
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
    # Closed frame with all borders at thickness 0.8
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()  # We'll use custom approach
  )

# ============================================================================
# STATISTICAL ANALYSIS & INTERPRETATION
# ============================================================================

# Overall statistics
overall_stats <- list(
  Gene = GENE_NAME,
  N_Patients = nrow(gene_summary),
  N_Replicates_Total = nrow(gene_pct_long),
  Mean_Abundance = round(mean(gene_summary$mean_pct), 4),
  Median_Abundance = round(median(gene_summary$mean_pct), 4),
  SD_Abundance = round(sd(gene_summary$mean_pct), 4),
  Min_Abundance = round(min(gene_summary$mean_pct), 4),
  Max_Abundance = round(max(gene_summary$mean_pct), 4),
  Range = round(max(gene_summary$mean_pct) - min(gene_summary$mean_pct), 4),
  CV_Overall = round((sd(gene_summary$mean_pct) / mean(gene_summary$mean_pct)) * 100, 2)
)

# High expressors (top 25%)
high_threshold <- quantile(gene_summary$mean_pct, 0.75)
high_expressors <- gene_summary |>
  filter(mean_pct >= high_threshold) |>
  select(Sample, mean_pct, ci_lower, ci_upper, cv_pct) |>
  arrange(desc(mean_pct))

# Low expressors (bottom 25%)
low_threshold <- quantile(gene_summary$mean_pct, 0.25)
low_expressors <- gene_summary |>
  filter(mean_pct <= low_threshold) |>
  select(Sample, mean_pct, ci_lower, ci_upper, cv_pct) |>
  arrange(mean_pct)

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
  select(Sample, Patient, mean_pct, sd_pct, se_pct, ci_lower, ci_upper, n, cv_pct) |>
  mutate(across(where(is.numeric), ~round(., 6)))
writeData(wb, "Patient_Details", patient_data)

# Sheet 3: High Expressors (Top 25%)
addWorksheet(wb, "High_Expressors_Top25")
high_exp_display <- high_expressors |>
  mutate(across(where(is.numeric), ~round(., 6)))
writeData(wb, "High_Expressors_Top25", high_exp_display)

# Sheet 4: Low Expressors (Bottom 25%)
addWorksheet(wb, "Low_Expressors_Bottom25")
low_exp_display <- low_expressors |>
  mutate(across(where(is.numeric), ~round(., 6)))
writeData(wb, "Low_Expressors_Bottom25", low_exp_display)

# Sheet 5: Interpretation & Findings
addWorksheet(wb, "Interpretation")
interpretation <- data.frame(
  Finding = c(
    "Gene Name",
    "Number of Patients",
    "Mean Relative Abundance",
    "Median Relative Abundance",
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
    paste0(round(overall_stats$Mean_Abundance, 4), "%"),
    paste0(round(overall_stats$Median_Abundance, 4), "%"),
    paste0(round(overall_stats$Range, 4), "%"),
    paste0(round(overall_stats$CV_Overall, 2), "%"),
    paste0(nrow(high_expressors), " patients (", 
           paste0(unique(high_expressors$Sample), collapse = ", "), ")"),
    paste0(nrow(low_expressors), " patients (", 
           paste0(unique(low_expressors$Sample), collapse = ", "), ")"),
    if_else(overall_stats$CV_Overall > 100, "High variability", 
            if_else(overall_stats$CV_Overall > 50, "Moderate variability", "Low variability")),
    paste0("The ", GENE_NAME, " gene shows ", 
           if_else(overall_stats$Range > 10, "substantial", "modest"),
           " variation across patients with a ", 
           if_else(nrow(high_expressors) > 5, "broad", "narrow"),
           " distribution of relative abundance.")
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
cat("\n✓ Files saved to:", OUTPUT_FOLDER, "\n\n")

# ============================================================================
# DISPLAY PLOT IN RSTUDIO PLOTS PANE
# ============================================================================
print(p_gene)
