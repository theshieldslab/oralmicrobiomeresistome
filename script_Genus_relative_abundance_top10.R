#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# GENUS-LEVEL RELATIVE ABUNDANCE - TOP 10 GENERA (28 SAMPLES)
# Simplified stacked bar chart with proper scale + Statistical Analysis
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(readxl)
library(ggplot2)
library(openxlsx2)

main_folder <- "/Users/macbookair/Desktop/16s may"
master_path <- file.path(main_folder, "16S_Master_Summary_UPDATED_28_samples_FINAL.xlsx")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_png <- file.path(main_folder, "Genus_relative_abundance_top10_28samples.png")
output_jpeg <- file.path(main_folder, "Genus_relative_abundance_top10_28samples.jpg")
output_pdf <- file.path(main_folder, "Genus_relative_abundance_top10_28samples.pdf")
output_excel <- file.path(main_folder, paste0("interpretation_Genus_abundance_top10_", timestamp, ".xlsx"))

cat("════════════════════════════════════════════════════════════════════════════\n")
cat("GENUS RELATIVE ABUNDANCE - TOP 10 GENERA\n")
cat("════════════════════════════════════════════════════════════════════════════\n\n")

# ── LOAD DATA ──────────────────────────────────────────────────────────────────
genus_data <- read_excel(master_path, sheet = "Genus_Abundance")
cat("Loaded Genus_Abundance:", nrow(genus_data), "rows\n\n")

# ── EXTRACT GENUS NAME (first word only) ──────────────────────────────────────
# Get the first column name (either "Sample" or "Genus")
first_col <- colnames(genus_data)[1]

genus_data <- genus_data %>%
  mutate(Genus_Simple = word(.data[[first_col]], 1))

cat("Example genus names extracted:\n")
print(head(genus_data %>% select(all_of(first_col), Genus_Simple), 10))
cat("\nUnique simple genus names:", n_distinct(genus_data$Genus_Simple), "\n\n")

# ── GET TOP 10 GENERA (by mean abundance) ──────────────────────────────────────
# Aggregate by simple genus name (sum all species within each genus)
genus_matrix <- genus_data %>%
  select(-all_of(first_col), -Genus_Simple) %>%
  as.matrix()

genus_means <- rowMeans(genus_matrix, na.rm = TRUE)
names(genus_means) <- genus_data$Genus_Simple

# Get unique top 10
top10_unique <- names(sort(genus_means, decreasing = TRUE))
top10_unique <- unique(top10_unique[1:15])[1:10]  # Get unique, ensuring 10

cat("Top 10 unique genera by mean abundance:\n")
for (i in seq_along(top10_unique)) {
  mean_val <- mean(genus_means[names(genus_means) == top10_unique[i]], na.rm = TRUE)
  cat(sprintf("%2d. %s (%.2f%%)\n", i, top10_unique[i], mean_val*100))
}
cat("\n")

# ── AGGREGATE TO GENUS LEVEL (simple) ──────────────────────────────────────────
# Sum all species within each genus
genus_agg <- genus_data %>%
  select(-all_of(first_col)) %>%
  group_by(Genus_Simple) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")

# Convert to long format
plot_data <- genus_agg %>%
  filter(Genus_Simple %in% top10_unique) %>%
  pivot_longer(-Genus_Simple, names_to = "Patient", values_to = "Abundance") %>%
  mutate(
    Genus_Simple = factor(Genus_Simple, levels = top10_unique),
    Patient = factor(Patient, levels = paste0("Patient ", 1:28)),
    Abundance_Pct = Abundance * 100
  ) %>%
  arrange(Patient, desc(Genus_Simple))

cat("Plot data prepared:\n")
cat("  Rows:", nrow(plot_data), "\n")
cat("  Patients:", n_distinct(plot_data$Patient), "\n")
cat("  Genera:", n_distinct(plot_data$Genus_Simple), "\n\n")

# Verify stacking
sample_check <- plot_data %>%
  group_by(Patient) %>%
  summarise(Total_Pct = sum(Abundance_Pct, na.rm = TRUE), .groups = "drop")

cat("Sample totals (% of top 10 genera only):\n")
print(summary(sample_check$Total_Pct))
cat("\n")

# ── COLOR PALETTE ──────────────────────────────────────────────────────────────
# Vibrant distinct palette - extracted from reference image
colors_10 <- c(
  "Streptococcus" = "#1F77B4",      # Dark blue
  "Selenomonas" = "#17BECF",        # Cyan/teal
  "Neisseria" = "#FF7F0E",          # Orange
  "Campylobacter" = "#2CA02C",      # Green
  "Veillonella" = "#9467BD",        # Purple
  "Aggregatibacter" = "#D62728",    # Red
  "Cardiobacterium" = "#8C564B",    # Brown
  "Fusobacterium" = "#E377C2",      # Pink/magenta
  "Gemella" = "#BCBD22",            # Yellow-green
  "Actinomyces" = "#7F7F7F"         # Gray
)

# ── CREATE STACKED BAR CHART ───────────────────────────────────────────────────
p_genus <- ggplot(plot_data, aes(x = Patient, y = Abundance_Pct, fill = Genus_Simple)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(
    values = colors_10,
    name = "Genus",
    guide = guide_legend(nrow = 5, title.position = "top")
  ) +
  scale_y_continuous(
    name = "Relative Abundance (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    expand = c(0, 0)
  ) +
  labs(
    title = "Genus-level Relative Abundance (Top 10 genera - 28 Samples)",
    x = "Patient Sample"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 10)),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 0.8)
  )

# ── SAVE OUTPUTS ───────────────────────────────────────────────────────────────
cat("Saving plots...\n")
ggsave(output_png, plot = p_genus, width = 13, height = 6, dpi = 600, bg = "white")
ggsave(output_jpeg, plot = p_genus, width = 13, height = 6, dpi = 300, bg = "white")
ggsave(output_pdf, plot = p_genus, width = 13, height = 6, bg = "white")

cat("\n✓ PNG saved:", output_png, "\n")
cat("✓ JPEG saved:", output_jpeg, "\n")
cat("✓ PDF saved:", output_pdf, "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# STATISTICAL ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

# Load caries risk metadata
caries_risk_sheet <- read_excel(master_path, sheet = "caries risk")
num_patients <- 28
caries_clean <- as.character(unlist(caries_risk_sheet[[1]][1:num_patients]))
caries_clean <- tolower(trimws(caries_clean))

# Calculate Shannon diversity for each patient
shannon_vals <- numeric(28)
patient_cols <- paste0("Patient ", 1:28)

for (i in 1:28) {
  # Get abundances for this patient
  patient_data <- genus_agg %>%
    select(Genus_Simple, all_of(paste0("Patient ", i)))
  
  abund <- patient_data[[paste0("Patient ", i)]]
  abund <- abund[abund > 0]
  shannon_vals[i] <- -sum(abund * log(abund), na.rm = TRUE)
}

# Create diversity dataframe
diversity_df <- data.frame(
  Patient = patient_cols,
  Shannon = shannon_vals,
  Caries_Risk = factor(caries_clean, levels = c("low", "high")),
  stringsAsFactors = FALSE
)

# Wilcoxon test
wt <- wilcox.test(Shannon ~ Caries_Risk, data = diversity_df, exact = FALSE)
p_fmt <- ifelse(wt$p.value < 0.001, "p < 0.001",
                paste0("p = ", format(round(wt$p.value, 4), nsmall = 4)))

# Summary statistics
diversity_summary <- diversity_df %>%
  group_by(Caries_Risk) %>%
  summarise(
    n = n(),
    Mean = round(mean(Shannon), 3),
    Median = round(median(Shannon), 3),
    SD = round(sd(Shannon), 3),
    Min = round(min(Shannon), 3),
    Max = round(max(Shannon), 3),
    .groups = "drop"
  )

cat("════════════════════════════════════════════════════════════════════════════\n")
cat("ALPHA DIVERSITY ANALYSIS (SHANNON INDEX)\n")
cat("════════════════════════════════════════════════════════════════════════════\n\n")
cat("Wilcoxon test:", p_fmt, "\n\n")
print(diversity_summary)
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# PREPARE INTERPRETATION EXCEL FILE
# ══════════════════════════════════════════════════════════════════════════════

# Sheet 1: Top 10 genera summary
top10_summary <- tibble(
  Rank = 1:10,
  Genus = top10_unique,
  Mean_Abundance_Pct = sapply(top10_unique, function(g) {
    mean(genus_agg$Genus_Simple[genus_agg$Genus_Simple == g], na.rm = TRUE) * 100
  })
)

# Recalculate with proper aggregation
top10_summary <- tibble(
  Rank = 1:10,
  Genus = top10_unique,
  Mean_Abundance_Pct = sapply(top10_unique, function(g) {
    genus_data_subset <- genus_data %>% filter(Genus_Simple == g) %>% select(-Genus, -Genus_Simple)
    mean(as.matrix(genus_data_subset), na.rm = TRUE) * 100
  })
)

# Sheet 2: Diversity statistics
diversity_export <- diversity_df
colnames(diversity_export) <- c("Patient", "Shannon_Index", "Caries_Risk")

# Sheet 3: Diversity summary by group
diversity_summary_export <- diversity_summary

# Sheet 4: Statistical test details
stat_test <- tibble(
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

# Sheet 5: Full interpretation text
overall_interpretation <- paste0(
  "GENUS-LEVEL RELATIVE ABUNDANCE ANALYSIS - TOP 10 GENERA\n",
  "Dataset: 28 oral biofilm samples\n\n",
  
  "RELATIVE ABUNDANCE FINDINGS:\n",
  "The microbial community composition was analyzed at the genus level across all samples. ",
  "The top 10 dominant genera were identified by mean abundance:\n\n"
)

for (i in 1:nrow(top10_summary)) {
  overall_interpretation <- paste0(
    overall_interpretation,
    sprintf(
      "%2d. %s: %.2f%% mean abundance\n",
      i, top10_summary$Genus[i],
      top10_summary$Mean_Abundance_Pct[i]
    )
  )
}

overall_interpretation <- paste0(
  overall_interpretation,
  "\nALPHA DIVERSITY ANALYSIS:\n",
  "Shannon diversity index was calculated for each patient sample to assess microbial richness and evenness.\n\n",
  
  "Summary statistics:\n",
  sprintf("  Overall mean Shannon index: %.3f\n", mean(diversity_df$Shannon)),
  sprintf("  Overall median Shannon index: %.3f\n", median(diversity_df$Shannon)),
  sprintf("  Range: %.3f - %.3f\n\n", min(diversity_df$Shannon), max(diversity_df$Shannon)),
  
  "By Caries Risk Group:\n",
  sprintf("  Low caries risk (n=%d): Mean Shannon = %.3f ± %.3f\n",
          diversity_summary$n[diversity_summary$Caries_Risk == "low"],
          diversity_summary$Mean[diversity_summary$Caries_Risk == "low"],
          diversity_summary$SD[diversity_summary$Caries_Risk == "low"]),
  sprintf("  High caries risk (n=%d): Mean Shannon = %.3f ± %.3f\n\n",
          diversity_summary$n[diversity_summary$Caries_Risk == "high"],
          diversity_summary$Mean[diversity_summary$Caries_Risk == "high"],
          diversity_summary$SD[diversity_summary$Caries_Risk == "high"]),
  
  "STATISTICAL TEST RESULTS:\n",
  "Wilcoxon rank-sum test was used to compare Shannon diversity between caries risk groups.\n",
  sprintf("Test statistic (W): %.2f\n", wt$statistic),
  sprintf("p-value: %.4f\n", wt$p.value),
  sprintf("Result: %s\n\n",
          ifelse(wt$p.value < 0.05,
                 "There is a statistically significant difference in Shannon diversity between low and high caries risk groups.",
                 "There is no statistically significant difference in Shannon diversity between low and high caries risk groups.")),
  
  "CONCLUSION:\n",
  "This analysis reveals the composition and diversity of the oral microbiome focusing on the top 10 dominant genera. ",
  "The identification of dominant genera and assessment of alpha diversity provide insights into ",
  "microbial community structure and its relationship with caries risk status."
)

interpretation_df <- data.frame(Text = overall_interpretation)

# Create and save workbook
cat("Creating interpretation Excel file...\n\n")
wb <- openxlsx2::wb_workbook()

wb$add_worksheet("Top10_Genera_Summary")
wb$add_data(x = top10_summary)

wb$add_worksheet("Diversity_All_Patients")
wb$add_data(x = diversity_export)

wb$add_worksheet("Diversity_By_Group")
wb$add_data(x = diversity_summary_export)

wb$add_worksheet("Statistical_Tests")
wb$add_data(x = stat_test)

wb$add_worksheet("Interpretation")
wb$add_data(x = interpretation_df)

wb$save(output_excel)

cat("✓ Interpretation Excel:", output_excel, "\n\n")

cat("════════════════════════════════════════════════════════════════════════════\n")
cat("GENUS ABUNDANCE ANALYSIS COMPLETE - TOP 10 GENERA!\n")
cat("════════════════════════════════════════════════════════════════════════════\n")
cat("Outputs generated:\n")
cat("  ✓ PNG (600 dpi publication quality)\n")
cat("  ✓ JPEG (300 dpi)\n")
cat("  ✓ PDF (vector format)\n")
cat("  ✓ Excel with 5 sheets (top 10 genera, diversity, statistics, interpretation)\n")
cat("════════════════════════════════════════════════════════════════════════════\n\n")

# Display
print(p_genus)

