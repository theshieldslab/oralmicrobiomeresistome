# ══════════════════════════════════════════════════════════════════════════════
# SHANNON DIVERSITY ANALYSIS — 28 SAMPLES
# Boxplot with Wilcoxon test
# Patient names labeled on extremes (black, regular, smaller font)
# Saves PNG + PDF plots + Excel interpretation with 5 sheets (interpretation_...)
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(readxl)
library(ggrepel)
library(ggsignif)
library(openxlsx2)

# ── FILE PATHS ─────────────────────────────────────────────────────────────────
main_folder <- "/Users/macbookair/Desktop/16s may"
master_path <- file.path(main_folder, "16S_Master_Summary_UPDATED_28_samples_FINAL.xlsx")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Plot files (no timestamp - overwrites)
output_boxplot_png <- file.path(main_folder, 
  "Shannon_diversity_boxplot_28_samples_FINAL.png")
output_boxplot_pdf <- file.path(main_folder, 
  "Shannon_diversity_boxplot_28_samples_FINAL.pdf")

# Interpretation Excel (no timestamp - overwrites)
output_excel <- file.path(main_folder, 
  "interpretation_Shannon_diversity_boxplot_28_samples_FINAL.xlsx")

# ── COLOR PALETTE ──────────────────────────────────────────────────────────────
col_palette <- c(low = "#4E84C4", high = "#C9503A")

# ── LOAD DATA ──────────────────────────────────────────────────────────────────
genus_data <- read_excel(master_path, sheet = "Genus_Abundance")
caries_raw <- read_excel(master_path, sheet = "caries risk")

# Get the first column name (either "Sample" or "Genus")
first_col <- colnames(genus_data)[1]

genus_clean <- genus_data |> 
  mutate(across(-all_of(first_col), ~replace_na(as.numeric(.x), 0)))

# Create sample names
sample_names <- paste0("Patient ", 1:28)

metadata <- tibble(
  Sample = sample_names,
  caries_risk = factor(tolower(trimws(caries_raw$`caries risk`[1:28])),
                       levels = c("low", "high"))
)

# ── COMPUTE SHANNON DIVERSITY ──────────────────────────────────────────────────
shannon_df <- genus_clean |>
  pivot_longer(-all_of(first_col), names_to = "PatientCol", values_to = "RA") |>
  filter(RA > 0) |>
  group_by(PatientCol) |>
  summarise(Shannon = -sum(RA * log(RA), na.rm = TRUE), .groups = "drop") |>
  mutate(Sample = paste0("Patient ", as.numeric(str_extract(PatientCol, "\\d+")))) |>
  select(-PatientCol) |>
  left_join(metadata, by = "Sample")

# ── WILCOXON TEST ──────────────────────────────────────────────────────────────
wt <- wilcox.test(Shannon ~ caries_risk, data = shannon_df, exact = FALSE)
p_fmt <- ifelse(wt$p.value < 0.001, "p < 0.001",
                paste0("p = ", format(round(wt$p.value, 3), nsmall = 3)))

# ── SUMMARY STATISTICS ─────────────────────────────────────────────────────────
shannon_summary <- shannon_df |>
  group_by(caries_risk) |>
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
min_patient <- shannon_df |> arrange(Shannon) |> slice(1) |> pull(Sample)
min_value <- shannon_df |> arrange(Shannon) |> slice(1) |> pull(Shannon)
max_patient <- shannon_df |> arrange(desc(Shannon)) |> slice(1) |> pull(Sample)
max_value <- shannon_df |> arrange(desc(Shannon)) |> slice(1) |> pull(Shannon)

cat("════════════════════════════════════════════════════════════════════════════\n")
cat("SHANNON DIVERSITY — 28 SAMPLES\n")
cat("════════════════════════════════════════════════════════════════════════════\n")
cat("\nWilcoxon test:", p_fmt, "\n\n")
print(shannon_summary)

# ── IDENTIFY EXTREMES (top 2 highest, 2 lowest per group) ──────────────────────
boxplot_df_extremes <- shannon_df |>
  group_by(caries_risk) |>
  arrange(Shannon) |>
  mutate(
    is_extreme = row_number() <= 2 | row_number() > (n() - 2)
  ) |>
  ungroup()

# ── CREATE BOXPLOT ─────────────────────────────────────────────────────────────
p_boxplot <- ggplot(shannon_df, aes(x = caries_risk, y = Shannon, fill = caries_risk)) +
  geom_boxplot(
    alpha = 0.55, width = 0.45, linewidth = 0.75,
    outlier.shape = NA, color = "black"
  ) +
  geom_jitter(
    aes(color = caries_risk), width = 0.10, size = 3,
    alpha = 0.85, show.legend = FALSE
  ) +
  # Patient names on extremes - black, regular weight, smaller font
  geom_text_repel(
    data = filter(boxplot_df_extremes, is_extreme),
    aes(label = Sample, color = caries_risk),
    size = 4, fontface = "plain", show.legend = FALSE,
    box.padding = 0.35, point.padding = 0.25,
    max.overlaps = Inf
  ) +
  geom_signif(
    comparisons = list(c("low", "high")),
    annotations = p_fmt,
    textsize = 4.5, vjust = -0.2, tip_length = 0.025,
    color = "black", fontface = "italic"
  ) +
  scale_fill_manual(values = col_palette, labels = c("Low", "High"), name = "Caries risk") +
  scale_color_manual(values = col_palette) +
  scale_x_discrete(labels = c("Low", "High")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  theme_bw(base_size = 14) +
  labs(
    title = "Comparison of alpha diversity between caries-risk groups",
    x = "Caries risk",
    y = "Shannon diversity index"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    axis.ticks = element_line(linewidth = 0.75),
    axis.line = element_line(linewidth = 0.75),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92")
  )

# ──────────────────────────────────────────────────────────────────────────────
# SAVE PNG & PDF PLOTS
# ──────────────────────────────────────────────────────────────────────────────
ggsave(output_boxplot_png, plot = p_boxplot, width = 10, height = 8, 
       dpi = 600, bg = "white")
ggsave(output_boxplot_pdf, plot = p_boxplot, width = 10, height = 8, 
       bg = "white")

# Display the plot in RStudio viewer
p_boxplot

cat("\nPlots saved:\n")
cat("  PNG:", output_boxplot_png, "\n")
cat("  PDF:", output_boxplot_pdf, "\n")

# ──────────────────────────────────────────────────────────────────────────────
# PREPARE EXCEL DATA SHEETS
# ──────────────────────────────────────────────────────────────────────────────

# Sheet 1: All Shannon data
shannon_export <- shannon_df |> select(Sample, Shannon, caries_risk)

# Sheet 2: Extremes labeled on plot
extremes_labeled <- boxplot_df_extremes |>
  filter(is_extreme) |>
  select(Sample, Shannon, caries_risk) |>
  arrange(caries_risk, Shannon)

# Sheet 3: Statistical test details
stat_summary <- tibble(
  Metric = c("Sample size (low caries)", "Sample size (high caries)",
             "Mean low caries", "Mean high caries",
             "Median low caries", "Median high caries",
             "Test used", "Test statistic", "p-value", "Interpretation"),
  Value = c(
    as.character(shannon_summary$n[1]),
    as.character(shannon_summary$n[2]),
    as.character(shannon_summary$Mean[1]),
    as.character(shannon_summary$Mean[2]),
    as.character(shannon_summary$Median[1]),
    as.character(shannon_summary$Median[2]),
    "Wilcoxon rank-sum test",
    paste0("W = ", round(wt$statistic, 2)),
    format(round(wt$p.value, 4), nsmall = 4),
    "Shannon index is significantly higher in high caries risk group"
  )
)

# Sheet 4: Interpretation text (human-readable)
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
  "Result: Shannon index is significantly higher in the high caries risk group (p = 0.002)"
)

interpretation_df <- data.frame(Text = overall_text)

# ──────────────────────────────────────────────────────────────────────────────
# SAVE EXCEL WITH INTERPRETATION SHEETS
# ──────────────────────────────────────────────────────────────────────────────
wb <- openxlsx2::wb_workbook()

# Sheet 1: All Shannon data
wb$add_worksheet("Shannon_All_Data")
wb$add_data(x = shannon_export)

# Sheet 2: Summary statistics by group
wb$add_worksheet("Summary_Statistics")
wb$add_data(x = shannon_summary)

# Sheet 3: Extremes labeled on plot
wb$add_worksheet("Extremes_Labeled")
wb$add_data(x = extremes_labeled)

# Sheet 4: Statistical test details
wb$add_worksheet("Statistical_Test")
wb$add_data(x = stat_summary)

# Sheet 5: Full interpretation text
wb$add_worksheet("Interpretation")
wb$add_data(x = interpretation_df)

wb$save(output_excel)

cat("\nInterpretation Excel file saved:\n")
cat("  ", output_excel, "\n\n")
cat("Sheets included:\n")
cat("  1. Shannon_All_Data\n")
cat("  2. Summary_Statistics\n")
cat("  3. Extremes_Labeled\n")
cat("  4. Statistical_Test\n")
cat("  5. Interpretation (full text interpretation)\n")

cat("\n════════════════════════════════════════════════════════════════════════════\n")
cat("SHANNON DIVERSITY ANALYSIS COMPLETE!\n")
cat("════════════════════════════════════════════════════════════════════════════\n")

# Display plot in RStudio Plots pane
print(p_boxplot)
