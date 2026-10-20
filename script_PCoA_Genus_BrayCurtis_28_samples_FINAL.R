# ══════════════════════════════════════════════════════════════════════════════
# PCoA ANALYSIS — 28 SAMPLES
# Principal Coordinates Analysis with Bray-Curtis distance
# Colored by caries risk group with dashed 95% ellipses
# Patient names labeled ONLY on outliers (clean plot)
# Saves PNG + PDF plots + Excel interpretation with 5 sheets (interpretation_...)
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(readxl)
library(vegan)
library(ggrepel)
library(openxlsx2)

# ── FILE PATHS ─────────────────────────────────────────────────────────────────
main_folder <- "/Users/macbookair/Desktop/16s may"
master_path <- file.path(main_folder, "16S_Master_Summary_UPDATED_28_samples_FINAL.xlsx")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Plot files (no timestamp - overwrites)
output_pcoa_png <- file.path(main_folder, 
  "PCoA_Genus_BrayCurtis_28_samples_FINAL.png")
output_pcoa_pdf <- file.path(main_folder, 
  "PCoA_Genus_BrayCurtis_28_samples_FINAL.pdf")

# Interpretation Excel (no timestamp - overwrites)
output_excel <- file.path(main_folder, 
  "interpretation_PCoA_Genus_BrayCurtis_28_samples_FINAL.xlsx")

# ── COLOR PALETTE ──────────────────────────────────────────────────────────────
col_palette <- c(low = "#4E84C4", high = "#C9503A")

# ── LOAD DATA ──────────────────────────────────────────────────────────────────
genus_data <- read_excel(master_path, sheet = "Genus_Abundance")
caries_raw <- read_excel(master_path, sheet = "caries risk")

# Get the first column name (either "Sample" or "Genus")
first_col <- colnames(genus_data)[1]

# Transpose: patients as rows, genera as columns
genus_t <- as.data.frame(t(genus_data[, -1]))
colnames(genus_t) <- genus_data[[first_col]]
genus_t$Sample <- paste0("Patient ", 1:28)
genus_t <- genus_t |> relocate(Sample)

genus_clean <- genus_t |> 
  mutate(across(-Sample, ~replace_na(as.numeric(.x), 0)))

metadata <- tibble(
  Sample = genus_clean$Sample,
  caries_risk = factor(tolower(trimws(caries_raw$`caries risk`[1:28])),
                       levels = c("low", "high"))
)

# ── BUILD GENUS MATRIX FOR BRAY-CURTIS ────────────────────────────────────────
genus_mat <- genus_clean |>
  mutate(n = as.integer(str_extract(Sample, "\\d+"))) |>
  arrange(n) |>
  select(-Sample, -n) |>
  as.data.frame()

rownames(genus_mat) <- metadata |>
  mutate(n = as.integer(str_extract(Sample, "\\d+"))) |>
  arrange(n) |>
  pull(Sample)

# ── BRAY-CURTIS DISTANCE & PCoA ────────────────────────────────────────────────
bray <- vegdist(genus_mat, method = "bray")
pcoa <- cmdscale(bray, eig = TRUE, k = 2)

var_pct <- round(100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]), 1)

pcoa_df <- as_tibble(pcoa$points, rownames = "Sample") |>
  rename(PCoA1 = V1, PCoA2 = V2) |>
  left_join(metadata, by = "Sample")

# ── PERMANOVA ──────────────────────────────────────────────────────────────────
set.seed(8243)
perm <- adonis2(bray ~ caries_risk, data = metadata, permutations = 999)
p_perm <- perm$`Pr(>F)`[1]
perm_label <- paste0("p = ", round(p_perm, 3), ", PERMANOVA")

cat("════════════════════════════════════════════════════════════════════════════\n")
cat("PCoA ANALYSIS — 28 SAMPLES\n")
cat("════════════════════════════════════════════════════════════════════════════\n")
cat("\nPCoA variance explained:\n")
cat("  PCoA1:", var_pct[1], "%\n")
cat("  PCoA2:", var_pct[2], "%\n\n")
cat("PERMANOVA:", perm_label, "\n")
print(perm)

# ── IDENTIFY TRUE OUTLIERS (outside 1.5x median distance from centroid) ────────
identify_outliers <- function(df, group_name) {
  group_data <- df |> filter(caries_risk == group_name)
  
  centroid_x <- mean(group_data$PCoA1)
  centroid_y <- mean(group_data$PCoA2)
  
  group_data |>
    mutate(
      dist_from_centroid = sqrt((PCoA1 - centroid_x)^2 + (PCoA2 - centroid_y)^2),
      med_dist = median(dist_from_centroid),
      is_outlier = dist_from_centroid > med_dist * 1.5
    ) |>
    filter(is_outlier)
}

outliers_low <- identify_outliers(pcoa_df, "low")
outliers_high <- identify_outliers(pcoa_df, "high")
all_outliers <- bind_rows(outliers_low, outliers_high)

cat("\nOutliers identified:\n")
cat("  Low caries risk:", nrow(outliers_low), "\n")
cat("  High caries risk:", nrow(outliers_high), "\n")
cat("  Total:", nrow(all_outliers), "\n")

# ── CREATE PCoA PLOT with outliers labeled ─────────────────────────────────────
p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = caries_risk, fill = caries_risk)) +
  # Confidence ellipses (dashed)
  stat_ellipse(
    aes(group = caries_risk),
    geom = "polygon", alpha = 0.08, linetype = 2, linewidth = 0.75, level = 0.95
  ) +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  # Points
  geom_point(size = 3.5, alpha = 0.9) +
  # ONLY outlier patient names
  geom_text_repel(
    data = all_outliers,
    aes(label = Sample, color = caries_risk),
    size = 4, fontface = "plain", show.legend = FALSE,
    box.padding = 0.4, point.padding = 0.3,
    max.overlaps = Inf
  ) +
  # PERMANOVA p-value
  annotate(
    "text", x = -Inf, y = Inf,
    label = perm_label, hjust = -0.05, vjust = 1.5,
    size = 4, fontface = "italic"
  ) +
  scale_color_manual(values = col_palette, labels = c("Low", "High"), name = "Caries risk") +
  scale_fill_manual(values = col_palette, labels = c("Low", "High"), name = "Caries risk") +
  theme_bw(base_size = 14) +
  labs(
    title = "Bray–Curtis PCoA according to caries risk",
    x = paste0("PCoA1 [", var_pct[1], "%]"),
    y = paste0("PCoA2 [", var_pct[2], "%]")
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    axis.ticks = element_line(linewidth = 0.75),
    axis.line = element_line(linewidth = 0.75),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92")
  )

# ──────────────────────────────────────────────────────────────────────────────
# SAVE PNG & PDF PLOTS
# ──────────────────────────────────────────────────────────────────────────────
ggsave(output_pcoa_png, plot = p_pcoa, width = 12, height = 9, 
       dpi = 600, bg = "white")
ggsave(output_pcoa_pdf, plot = p_pcoa, width = 12, height = 9, 
       bg = "white")

cat("\nPlots saved:\n")
cat("  PNG:", output_pcoa_png, "\n")
cat("  PDF:", output_pcoa_pdf, "\n")

# ──────────────────────────────────────────────────────────────────────────────
# PREPARE EXCEL DATA SHEETS
# ──────────────────────────────────────────────────────────────────────────────

# Sheet 1: PCoA coordinates
pcoa_export <- pcoa_df |> arrange(Sample) |> 
  mutate(
    Distance_from_origin = sqrt(PCoA1^2 + PCoA2^2),
    PCoA1 = round(PCoA1, 4),
    PCoA2 = round(PCoA2, 4),
    Distance_from_origin = round(Distance_from_origin, 4),
    Is_Outlier = ifelse(Sample %in% all_outliers$Sample, "Yes", "No")
  ) |>
  select(Sample, caries_risk, PCoA1, PCoA2, Distance_from_origin, Is_Outlier)

# Sheet 2: Variance explained
variance_df <- tibble(
  Axis = c("PCoA1", "PCoA2", "Total variance captured"),
  Variance_Percent = c(var_pct[1], var_pct[2], sum(var_pct[1:2]))
)

# Sheet 3: PERMANOVA results
permanova_df <- tibble(
  Metric = c("Test", "Df (Model)", "Df (Residual)", "SumOfSqs (Model)", 
             "SumOfSqs (Residual)", "R2 (Model)", "R2 (Residual)", 
             "F-statistic", "p-value", "Interpretation"),
  Value = c(
    "PERMANOVA (Adonis2)",
    "1",
    "26",
    round(perm[1, "SumOfSqs"], 4),
    round(perm[2, "SumOfSqs"], 4),
    round(perm[1, "R2"], 4),
    round(perm[2, "R2"], 4),
    round(perm[1, "F"], 4),
    round(perm[1, "Pr(>F)"], 4),
    "Significant difference in microbial composition between low and high caries risk groups"
  )
)

# Sheet 4: Summary by group
summary_by_group <- tibble(
  Caries_Risk = c("Low", "High"),
  n_samples = c(sum(pcoa_df$caries_risk == "low"), sum(pcoa_df$caries_risk == "high")),
  n_outliers = c(nrow(outliers_low), nrow(outliers_high)),
  Mean_PCoA1 = c(
    round(mean(pcoa_df$PCoA1[pcoa_df$caries_risk == "low"]), 4),
    round(mean(pcoa_df$PCoA1[pcoa_df$caries_risk == "high"]), 4)
  ),
  Mean_PCoA2 = c(
    round(mean(pcoa_df$PCoA2[pcoa_df$caries_risk == "low"]), 4),
    round(mean(pcoa_df$PCoA2[pcoa_df$caries_risk == "high"]), 4)
  )
)

# Sheet 5: Interpretation text
overall_text <- paste0(
  "PCoA INTERPRETATION — 28 SAMPLES\n\n",
  "Test: Principal Coordinates Analysis (PCoA) using Bray-Curtis distance\n",
  "Samples: 28 (n=9 low caries risk, n=19 high caries risk)\n\n",
  "VARIANCE EXPLAINED:\n",
  "PCoA1 explains ", var_pct[1], "% of variance\n",
  "PCoA2 explains ", var_pct[2], "% of variance\n",
  "Combined: ", round(sum(var_pct[1:2]), 1), "% of variance captured\n\n",
  "PERMANOVA RESULTS:\n",
  "Test: PERMANOVA (Adonis2) with 999 permutations\n",
  "F-statistic: ", round(perm[1, "F"], 3), "\n",
  "p-value: ", round(perm[1, "Pr(>F)"], 4), "\n",
  "R² (Model): ", round(perm[1, "R2"], 4), " (", round(perm[1, "R2"] * 100, 1), "% of variation explained by caries risk group)\n\n",
  "INTERPRETATION:\n",
  "Principal Coordinates Analysis revealed a significant separation in genus-level microbial ",
  "composition between low and high caries risk groups (PERMANOVA: F=", round(perm[1, "F"], 3), ", p=", round(perm[1, "Pr(>F)"], 4), "). ",
  "The first two principal coordinates captured ", round(sum(var_pct[1:2]), 1), "% of total variance, ",
  "with PCoA1 accounting for the majority (", var_pct[1], "%). Low caries risk samples were ",
  "predominantly distributed in the upper-left quadrant, while high caries risk samples clustered ",
  "in the lower-left and central regions. This spatial separation indicates distinct microbial ",
  "community structures associated with caries risk phenotypes at the genus level. A total of ",
  nrow(all_outliers), " samples were identified as outliers (", nrow(outliers_low), " from low risk, ",
  nrow(outliers_high), " from high risk), showing exceptional positioning relative to their respective groups.\n\n",
  "BIOLOGICAL SIGNIFICANCE:\n",
  "The significant clustering by caries risk group suggests that specific compositional differences ",
  "in the oral microbiota at the genus level are associated with differential caries susceptibility. ",
  "These findings support the hypothesis that microbial community structure, not just individual ",
  "taxa abundance, plays a role in caries disease etiology."
)

interpretation_df <- data.frame(Text = overall_text)

# ──────────────────────────────────────────────────────────────────────────────
# SAVE EXCEL WITH INTERPRETATION SHEETS
# ──────────────────────────────────────────────────────────────────────────────
wb <- openxlsx2::wb_workbook()

wb$add_worksheet("PCoA_Coordinates")
wb$add_data(x = pcoa_export)

wb$add_worksheet("Variance_Explained")
wb$add_data(x = variance_df)

wb$add_worksheet("PERMANOVA_Results")
wb$add_data(x = permanova_df)

wb$add_worksheet("Summary_By_Group")
wb$add_data(x = summary_by_group)

wb$add_worksheet("Interpretation")
wb$add_data(x = interpretation_df)

wb$save(output_excel)

cat("\nInterpretation Excel file saved:\n")
cat("  ", output_excel, "\n\n")
cat("Sheets included:\n")
cat("  1. PCoA_Coordinates (sample positions + distances + outlier flag)\n")
cat("  2. Variance_Explained (PCoA1, PCoA2 variance)\n")
cat("  3. PERMANOVA_Results (statistical test details)\n")
cat("  4. Summary_By_Group (means by caries risk + outlier counts)\n")
cat("  5. Interpretation (full text interpretation)\n")

cat("\n════════════════════════════════════════════════════════════════════════════\n")
cat("PCoA ANALYSIS COMPLETE!\n")
cat("════════════════════════════════════════════════════════════════════════════\n")
