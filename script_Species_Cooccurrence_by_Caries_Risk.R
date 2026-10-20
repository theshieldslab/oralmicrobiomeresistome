#!/usr/bin/env Rscript
# ============================================================
# SPECIES CO-OCCURRENCE/CO-PRESENCE BY CARIES RISK
# Which species appear together in LOW vs HIGH risk groups
# ComplexHeatmap correlation analysis
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(stats)

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("SPECIES CO-OCCURRENCE ANALYSIS BY CARIES RISK\n")
cat("Which species appear together in LOW vs HIGH risk groups\n")
cat("════════════════════════════════════════════════════════════════\n\n")

# ============================================================
# FILE PATHS
# ============================================================

main_folder <- "/Users/macbookair/Desktop/16s may"
master_path <- file.path(main_folder, "16S_Master_Summary_UPDATED_28_samples_FINAL.xlsx")

output_png_low <- file.path(main_folder, "Species_Cooccurrence_Low_Risk_FINAL.png")
output_pdf_low <- file.path(main_folder, "Species_Cooccurrence_Low_Risk_FINAL.pdf")

output_png_high <- file.path(main_folder, "Species_Cooccurrence_High_Risk_FINAL.png")
output_pdf_high <- file.path(main_folder, "Species_Cooccurrence_High_Risk_FINAL.pdf")

output_xlsx <- file.path(main_folder, "Species_Cooccurrence_by_Caries_Risk_Analysis_FINAL.xlsx")

# ============================================================
# 1. LOAD DATA
# ============================================================

cat("Loading data from master file...\n\n")

species_data <- read_excel(master_path, sheet = "Species_Abundance")
caries_data <- read_excel(master_path, sheet = "caries risk")

first_col <- colnames(species_data)[1]
patient_cols <- colnames(species_data)[-1]
n_samples <- length(patient_cols)

# ============================================================
# 2. PREPARE DATA: CREATE MATRICES BY CARIES RISK
# ============================================================

cat("Preparing data matrices by caries risk group...\n\n")

# Create metadata
metadata <- tibble(
  Patient = patient_cols,
  Patient_Num = as.numeric(str_extract(patient_cols, "\\d+")),
  Caries_Risk = tolower(trimws(caries_data$`caries risk`[1:n_samples]))
) |>
  arrange(Caries_Risk, Patient_Num)

# Normalize species data
species_normalized <- species_data |>
  mutate(across(-all_of(first_col), ~as.numeric(.) * 100))

colnames(species_normalized)[1] <- "Species"

# Select top species overall
top_species <- species_normalized |>
  rowwise() |>
  mutate(Mean = mean(c_across(-Species), na.rm = TRUE)) |>
  ungroup() |>
  arrange(desc(Mean)) |>
  slice_head(n = 15) |>
  pull(Species)

cat(paste0("✓ Selected top 15 species for co-occurrence analysis\n"))
cat(paste0("✓ Low caries risk samples: ", sum(metadata$Caries_Risk == "low"), "\n"))
cat(paste0("✓ High caries risk samples: ", sum(metadata$Caries_Risk == "high"), "\n\n"))

# ============================================================
# 3. CREATE CORRELATION MATRICES FOR EACH GROUP
# ============================================================

cat("Computing co-occurrence correlations...\n\n")

# Low risk group
low_risk_patients <- metadata$Patient[metadata$Caries_Risk == "low"]
low_risk_data <- species_normalized |>
  filter(Species %in% top_species) |>
  select(Species, all_of(low_risk_patients)) |>
  column_to_rownames("Species") |>
  as.matrix()

cor_low <- cor(t(low_risk_data), method = "pearson", use = "complete.obs")

# High risk group
high_risk_patients <- metadata$Patient[metadata$Caries_Risk == "high"]
high_risk_data <- species_normalized |>
  filter(Species %in% top_species) |>
  select(Species, all_of(high_risk_patients)) |>
  column_to_rownames("Species") |>
  as.matrix()

cor_high <- cor(t(high_risk_data), method = "pearson", use = "complete.obs")

cat("✓ Low risk correlation matrix computed\n")
cat("✓ High risk correlation matrix computed\n\n")

# ============================================================
# 4. COLOR SCALE FOR CORRELATION
# ============================================================

col_fun_corr <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#1a237e", "#4fc3f7", "#f5f5f5", "#ff7043", "#bf360c")
)

# ============================================================
# 5. CREATE HEATMAP FOR LOW RISK
# ============================================================

cat("Creating co-occurrence heatmap (LOW caries risk)...\n\n")

ht_low <- Heatmap(
  
  cor_low,
  
  name = "Pearson\nCorrelation",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  
  col = col_fun_corr,
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  row_names_side = "left",
  row_names_gp = gpar(
    fontsize = 9,
    fontface = "italic",
    col = "black"
  ),
  
  column_names_gp = gpar(
    fontsize = 9,
    fontface = "italic",
    col = "black"
  ),
  column_names_rot = 45,
  
  width = unit(ncol(cor_low) * 8, "mm"),
  height = unit(nrow(cor_low) * 8, "mm"),
  
  heatmap_legend_param = list(
    title = "Pearson\nCorrelation",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    at = c(-1, -0.5, 0, 0.5, 1),
    legend_height = unit(8, "cm"),
    legend_width = unit(1.2, "cm")
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- cor_low[i, j]
    if (!is.na(value) && abs(value) >= 0.3) {
      text_color <- ifelse(abs(value) >= 0.7, "white", "black")
      grid.text(
        sprintf("%.2f", value),
        x = x, y = y,
        gp = gpar(fontsize = 6.5, col = text_color)
      )
    }
  }
)

# Draw and save low risk
draw(
  ht_low,
  heatmap_legend_side = "right",
  padding = unit(c(2, 15, 12, 12), "mm")
)

grid.text(
  "Species Co-occurrence Network: LOW Caries Risk",
  x = unit(0.5, "npc"),
  y = unit(0.97, "npc"),
  gp = gpar(fontsize = 15, fontface = "bold", col = "black")
)

# Save PNG
png(output_png_low, width = 3600, height = 1700, res = 300)
draw(ht_low, heatmap_legend_side = "right", padding = unit(c(2, 15, 12, 12), "mm"))
grid.text("Species Co-occurrence Network: LOW Caries Risk",
          x = unit(0.5, "npc"), y = unit(0.97, "npc"),
          gp = gpar(fontsize = 15, fontface = "bold", col = "black"))
dev.off()

# Save PDF
pdf(output_pdf_low, width = 12, height = 5.5)
draw(ht_low, heatmap_legend_side = "right", padding = unit(c(2, 15, 12, 12), "mm"))
grid.text("Species Co-occurrence Network: LOW Caries Risk",
          x = unit(0.5, "npc"), y = unit(0.97, "npc"),
          gp = gpar(fontsize = 15, fontface = "bold", col = "black"))
dev.off()

cat("✓ Low risk heatmap saved (PNG + PDF)\n\n")

# ============================================================
# 6. CREATE HEATMAP FOR HIGH RISK
# ============================================================

cat("Creating co-occurrence heatmap (HIGH caries risk)...\n\n")

ht_high <- Heatmap(
  
  cor_high,
  
  name = "Pearson\nCorrelation",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  
  col = col_fun_corr,
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  row_names_side = "left",
  row_names_gp = gpar(
    fontsize = 9,
    fontface = "italic",
    col = "black"
  ),
  
  column_names_gp = gpar(
    fontsize = 9,
    fontface = "italic",
    col = "black"
  ),
  column_names_rot = 45,
  
  width = unit(ncol(cor_high) * 8, "mm"),
  height = unit(nrow(cor_high) * 8, "mm"),
  
  heatmap_legend_param = list(
    title = "Pearson\nCorrelation",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    at = c(-1, -0.5, 0, 0.5, 1),
    legend_height = unit(8, "cm"),
    legend_width = unit(1.2, "cm")
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- cor_high[i, j]
    if (!is.na(value) && abs(value) >= 0.3) {
      text_color <- ifelse(abs(value) >= 0.7, "white", "black")
      grid.text(
        sprintf("%.2f", value),
        x = x, y = y,
        gp = gpar(fontsize = 6.5, col = text_color)
      )
    }
  }
)

# Draw and save high risk
draw(
  ht_high,
  heatmap_legend_side = "right",
  padding = unit(c(2, 15, 12, 12), "mm")
)

grid.text(
  "Species Co-occurrence Network: HIGH Caries Risk",
  x = unit(0.5, "npc"),
  y = unit(0.97, "npc"),
  gp = gpar(fontsize = 15, fontface = "bold", col = "black")
)

# Save PNG
png(output_png_high, width = 3600, height = 1700, res = 300)
draw(ht_high, heatmap_legend_side = "right", padding = unit(c(2, 15, 12, 12), "mm"))
grid.text("Species Co-occurrence Network: HIGH Caries Risk",
          x = unit(0.5, "npc"), y = unit(0.97, "npc"),
          gp = gpar(fontsize = 15, fontface = "bold", col = "black"))
dev.off()

# Save PDF
pdf(output_pdf_high, width = 12, height = 5.5)
draw(ht_high, heatmap_legend_side = "right", padding = unit(c(2, 15, 12, 12), "mm"))
grid.text("Species Co-occurrence Network: HIGH Caries Risk",
          x = unit(0.5, "npc"), y = unit(0.97, "npc"),
          gp = gpar(fontsize = 15, fontface = "bold", col = "black"))
dev.off()

cat("✓ High risk heatmap saved (PNG + PDF)\n\n")

# ============================================================
# 7. EXTRACT STRONG CORRELATIONS
# ============================================================

cat("Extracting strong co-occurrence pairs...\n\n")

# Function to extract correlations
extract_corr_pairs <- function(cor_mat, threshold = 0.5, exclude_diag = TRUE) {
  pairs <- data.frame()
  for (i in 1:nrow(cor_mat)) {
    for (j in 1:ncol(cor_mat)) {
      if (exclude_diag && i >= j) next
      if (abs(cor_mat[i, j]) >= threshold) {
        pairs <- rbind(pairs, data.frame(
          Species1 = rownames(cor_mat)[i],
          Species2 = colnames(cor_mat)[j],
          Correlation = round(cor_mat[i, j], 3)
        ))
      }
    }
  }
  return(pairs |> arrange(desc(abs(Correlation))))
}

low_pairs <- extract_corr_pairs(cor_low, threshold = 0.5)
high_pairs <- extract_corr_pairs(cor_high, threshold = 0.5)

cat(paste0("✓ Strong positive correlations (r >= 0.5):\n"))
cat(paste0("  - Low risk: ", nrow(low_pairs), " species pairs\n"))
cat(paste0("  - High risk: ", nrow(high_pairs), " species pairs\n\n"))

# ============================================================
# 8. EXPORT TO EXCEL
# ============================================================

cat("Exporting results to Excel...\n\n")

sheets_export <- list(
  
  "Low_Risk_Correlations" = low_pairs,
  
  "High_Risk_Correlations" = high_pairs,
  
  "Low_Risk_Matrix" = as.data.frame(cor_low) |>
    rownames_to_column("Species"),
  
  "High_Risk_Matrix" = as.data.frame(cor_high) |>
    rownames_to_column("Species"),
  
  "Summary" = data.frame(
    Metric = c(
      "Analysis Type",
      "Correlation Method",
      "Species Analyzed",
      "Low Risk Samples",
      "High Risk Samples",
      "Strong Pairs (Low Risk, r >= 0.5)",
      "Strong Pairs (High Risk, r >= 0.5)"
    ),
    Value = c(
      "Species co-occurrence by caries risk",
      "Pearson correlation",
      length(top_species),
      sum(metadata$Caries_Risk == "low"),
      sum(metadata$Caries_Risk == "high"),
      nrow(low_pairs),
      nrow(high_pairs)
    ),
    stringsAsFactors = FALSE
  )
)

write_xlsx(sheets_export, output_xlsx)

cat("✓ Excel file saved\n\n")

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("════════════════════════════════════════════════════════════════\n")
cat("SPECIES CO-OCCURRENCE ANALYSIS COMPLETE\n")
cat("════════════════════════════════════════════════════════════════\n\n")

cat("OUTPUT FILES:\n")
cat(paste0("  ✓ ", basename(output_png_low), "\n"))
cat(paste0("  ✓ ", basename(output_pdf_low), "\n"))
cat(paste0("  ✓ ", basename(output_png_high), "\n"))
cat(paste0("  ✓ ", basename(output_pdf_high), "\n"))
cat(paste0("  ✓ ", basename(output_xlsx), "\n\n"))

cat("KEY INSIGHTS:\n")
cat(paste0("  • Low risk group: ", nrow(low_pairs), " strong co-occurrence pairs (r >= 0.5)\n"))
cat(paste0("  • High risk group: ", nrow(high_pairs), " strong co-occurrence pairs (r >= 0.5)\n\n"))

if (nrow(low_pairs) > 0) {
  cat("Top 3 Low Risk Pairs:\n")
  print(head(low_pairs, 3))
  cat("\n")
}

if (nrow(high_pairs) > 0) {
  cat("Top 3 High Risk Pairs:\n")
  print(head(high_pairs, 3))
  cat("\n")
}

cat("════════════════════════════════════════════════════════════════\n")
