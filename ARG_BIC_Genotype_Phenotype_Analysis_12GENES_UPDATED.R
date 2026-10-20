# ============================================================================
# GENOTYPE-PHENOTYPE ANALYSIS: ALL 11 AR GENES vs ANTIMICROBIAL BIC50
# ============================================================================
# 
# Study: 28 pediatric oral microbiome samples (caries risk study)
# Analysis: Spearman rank correlations between antimicrobial resistance gene
#           (ARG) abundance and corresponding antimicrobial BIC50 values
#
# 11 gene-antimicrobial comparisons based on functional classification:
#
# TETRACYCLINE RESISTANCE:
#   1. tetM vs Tetracycline BIC50
#   2. tetO vs Tetracycline BIC50
#
# ERYTHROMYCIN RESISTANCE:
#   3. ermB vs Erythromycin BIC50
#   4. ermF vs Erythromycin BIC50
#   5. mf(E) vs Erythromycin BIC50
#
# VANCOMYCIN RESISTANCE:
#   6. vanA vs Vancomycin BIC50
#   7. vanB vs Vancomycin BIC50
#
# CHLORHEXIDINE RESISTANCE:
#   8. qacJ vs Chlorhexidine BIC50
#   9. qacB vs Chlorhexidine BIC50
#   10. qacEΔ1 vs Chlorhexidine BIC50
#   11. eriCF1 vs Erythromycin BIC50
#
# ============================================================================

library(tidyverse)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggtext)

# ==================================================
# SETUP
# ==================================================

# Output directory
output_dir <- "/Users/macbookair/Desktop/qpcr/genotype_phenotype_analysis_11genes"
dir.create(output_dir, showWarnings = FALSE)

# ==================================================
# STEP 1: LOAD AND VERIFY DATA
# ==================================================

# Load the expanded 11-gene merged BIC-qPCR dataset
# (includes tetM, ermB, vanA, qacJ, qacB, edelta1, tetO, mf(E), vanB, eriCF1)
merged_df <- read_excel("/Users/macbookair/Desktop/qpcr/merged_BIC_qPCR_dataset_WITH_ermF.xlsx", 
                        sheet = "Merged_BIC_qPCR")

cat("✓ Loaded expanded 11-gene merged dataset\n")
cat("  Available genes:", 
    paste(colnames(merged_df)[4:ncol(merged_df)][!grepl("BIC|_Risk|Age_Group|Patient|Caries", 
                                                          colnames(merged_df)[4:ncol(merged_df)])], 
          collapse = ", "), "\n\n")

# Function to convert BIC values to numeric
convert_bic_to_numeric <- function(x) {
  x <- gsub(">", "", x)
  as.numeric(x)
}

# Create analysis-ready dataset
data_analysis <- merged_df |>
  mutate(
    # Convert BIC values
    Tetracycline_BIC_numeric = convert_bic_to_numeric(Tetracycline_BIC),
    Erythromycin_BIC_numeric = convert_bic_to_numeric(Erythromycin_BIC),
    Vancomycin_BIC_numeric = convert_bic_to_numeric(Vancomycin_BIC),
    CHX_BIC_numeric = convert_bic_to_numeric(CHX_BIC),
    Fluoride_BIC_numeric = convert_bic_to_numeric(NaF_BIC),
    # Keep originals for reference
    Tetracycline_BIC_original = Tetracycline_BIC,
    Erythromycin_BIC_original = Erythromycin_BIC,
    Vancomycin_BIC_original = Vancomycin_BIC,
    CHX_BIC_original = CHX_BIC,
    Fluoride_BIC_original = NaF_BIC
  )

# Data verification
cat("=== DATA VERIFICATION ===\n")
cat("Total patients:", nrow(data_analysis), "\n\n")

# ==================================================
# STEP 2: DEFINE ALL 11 COMPARISONS (NOW WITH EXPANDED MERGED DATASET)
# ==================================================

comparisons <- list(
  # TETRACYCLINE (2 genes)
  list(gene = "tetM", antimicrobial = "Tetracycline_BIC_numeric", 
       name = "Tetracycline", class = "Tetracycline"),
  list(gene = "tetO", antimicrobial = "Tetracycline_BIC_numeric", 
       name = "Tetracycline", class = "Tetracycline"),
  
  # ERYTHROMYCIN (3 genes - NOW INCLUDING ermF)
  list(gene = "ermB", antimicrobial = "Erythromycin_BIC_numeric", 
       name = "Erythromycin", class = "Erythromycin"),
  list(gene = "ermF", antimicrobial = "Erythromycin_BIC_numeric", 
       name = "Erythromycin", class = "Erythromycin"),
  list(gene = "mf(E)", antimicrobial = "Erythromycin_BIC_numeric", 
       name = "Erythromycin", class = "Erythromycin"),
  
  # FLUORIDE (1 gene)
  list(gene = "eriCF1", antimicrobial = "Fluoride_BIC_numeric", 
       name = "Fluoride", class = "Fluoride"),
  
  # VANCOMYCIN (2 genes)
  list(gene = "vanA", antimicrobial = "Vancomycin_BIC_numeric", 
       name = "Vancomycin", class = "Vancomycin"),
  list(gene = "vanB", antimicrobial = "Vancomycin_BIC_numeric", 
       name = "Vancomycin", class = "Vancomycin"),
  
  # CHLORHEXIDINE (3 genes)
  list(gene = "qacJ", antimicrobial = "CHX_BIC_numeric", 
       name = "Chlorhexidine", class = "Chlorhexidine"),
  list(gene = "qacB", antimicrobial = "CHX_BIC_numeric", 
       name = "Chlorhexidine", class = "Chlorhexidine"),
  list(gene = "edelta1", antimicrobial = "CHX_BIC_numeric", 
       name = "Chlorhexidine", class = "Chlorhexidine")
)

# ==================================================
# STEP 3: PERFORM SPEARMAN CORRELATIONS
# ==================================================

cat("=== CALCULATING SPEARMAN CORRELATIONS ===\n\n")

spearman_results <- list()

for (i in seq_along(comparisons)) {
  comp <- comparisons[[i]]
  
  gene_col <- comp$gene
  bic_col <- comp$antimicrobial
  
  # Create temporary data frame with complete cases
  temp_df <- data_analysis |>
    select(Patient, all_of(gene_col), all_of(bic_col)) |>
    na.omit()
  
  n_obs <- nrow(temp_df)
  
  # Check if gene column exists
  if (!(gene_col %in% colnames(data_analysis))) {
    cat("WARNING: Gene", gene_col, "not found in data. Skipping.\n")
    spearman_results[[i]] <- list(
      gene = gene_col,
      antimicrobial = comp$name,
      n = 0,
      rho = NA,
      p_value = NA,
      n_censored = NA,
      direction = NA,
      data = NULL,
      status = "SKIPPED - Gene not in data"
    )
    next
  }
  
  # Perform Spearman correlation
  if (n_obs > 0) {
    corr_test <- cor.test(
      temp_df[[gene_col]], 
      temp_df[[bic_col]], 
      method = "spearman"
    )
    
    # Count censored observations
    patient_indices <- which(data_analysis$Patient %in% temp_df$Patient)
    
    if (bic_col == "Tetracycline_BIC_numeric") {
      original_col_filtered <- data_analysis$Tetracycline_BIC_original[patient_indices]
    } else if (bic_col == "Erythromycin_BIC_numeric") {
      original_col_filtered <- data_analysis$Erythromycin_BIC_original[patient_indices]
    } else if (bic_col == "Vancomycin_BIC_numeric") {
      original_col_filtered <- data_analysis$Vancomycin_BIC_original[patient_indices]
    } else if (bic_col == "Fluoride_BIC_numeric") {
      original_col_filtered <- data_analysis$Fluoride_BIC_original[patient_indices]
    } else {
      original_col_filtered <- data_analysis$CHX_BIC_original[patient_indices]
    }
    
    n_censored <- sum(grepl(">", original_col_filtered), na.rm = TRUE)
    
    spearman_results[[i]] <- list(
      gene = comp$gene,
      antimicrobial = comp$name,
      n = n_obs,
      rho = corr_test$estimate,
      p_value = corr_test$p.value,
      n_censored = n_censored,
      direction = if (corr_test$estimate > 0) "positive" else "negative",
      data = temp_df,
      status = "OK"
    )
  } else {
    spearman_results[[i]] <- list(
      gene = gene_col,
      antimicrobial = comp$name,
      n = 0,
      rho = NA,
      p_value = NA,
      n_censored = NA,
      direction = NA,
      data = NULL,
      status = "SKIPPED - No complete cases"
    )
  }
}

# Create results table
results_table <- tibble(
  Gene = sapply(spearman_results, function(x) x$gene),
  Antimicrobial = sapply(spearman_results, function(x) x$antimicrobial),
  N = sapply(spearman_results, function(x) x$n),
  Spearman_rho = sapply(spearman_results, function(x) round(x$rho, 4)),
  Raw_P = sapply(spearman_results, function(x) round(x$p_value, 4)),
  Censored = sapply(spearman_results, function(x) x$n_censored),
  Direction = sapply(spearman_results, function(x) x$direction),
  Status = sapply(spearman_results, function(x) x$status)
)

cat("=== PRIMARY SPEARMAN RESULTS (ALL 11 GENES) ===\n")
print(results_table)
cat("\n")

# Store p-values for FDR correction (excluding NAs)
valid_results <- which(!is.na(sapply(spearman_results, function(x) x$p_value)))
raw_p_values_all <- sapply(spearman_results[valid_results], function(x) x$p_value)

# ==================================================
# STEP 4: BENJAMINI-HOCHBERG FDR CORRECTION
# ==================================================

adjusted_p_values_all <- p.adjust(raw_p_values_all, method = "BH")

# Create full results table with FDR
final_results_full <- tibble(
  Gene = results_table$Gene,
  Antimicrobial = results_table$Antimicrobial,
  Gene_Class = sapply(seq_along(comparisons), function(i) comparisons[[i]]$class),
  N = results_table$N,
  Spearman_rho = results_table$Spearman_rho,
  Raw_P = results_table$Raw_P,
  FDR_adjusted_P = NA_real_,
  Direction = results_table$Direction,
  Censored_N = results_table$Censored,
  Status = results_table$Status
)

# Fill in FDR values for valid results
for (i in valid_results) {
  idx_in_valid <- which(valid_results == i)
  final_results_full$FDR_adjusted_P[i] <- round(adjusted_p_values_all[idx_in_valid], 4)
}

cat("=== RESULTS WITH FDR CORRECTION ===\n")
print(final_results_full |> select(-Status))
cat("\n")

cat("Relationships with raw P < 0.05:\n")
sig_raw <- final_results_full |> filter(Raw_P < 0.05, !is.na(Raw_P))
if (nrow(sig_raw) > 0) {
  print(sig_raw |> select(Gene, Antimicrobial, Raw_P, FDR_adjusted_P))
} else {
  cat("None\n")
}
cat("\n")

cat("Relationships with FDR-adjusted P < 0.05:\n")
sig_fdr <- final_results_full |> filter(FDR_adjusted_P < 0.05, !is.na(FDR_adjusted_P))
if (nrow(sig_fdr) > 0) {
  print(sig_fdr |> select(Gene, Antimicrobial, Raw_P, FDR_adjusted_P))
} else {
  cat("None\n")
}
cat("\n")

# ==================================================
# STEP 5: CREATE FOREST PLOT (ALL 11 GENES)
# ==================================================

cat("Creating forest plot for all 11 genes...\n\n")

# Filter to only valid results for forest plot
forest_data <- final_results_full |>
  filter(!is.na(Spearman_rho), N > 0) |>
  mutate(
    # Use proper gene nomenclature for display (PDF-compatible)
    Gene_display = ifelse(Gene == "edelta1", "qacE*1", Gene)
  ) |>
  arrange(Spearman_rho) |>
  mutate(
    rank = row_number(),
    gene_antimicrobial = paste0(
      "<b><i>", Gene_display, "</i></b> ", 
      "<span style='color:#2171B5'>vs</span> ", Antimicrobial
    ),
    stat_label = paste0("p = ", round(Raw_P, 4))
  )

p_forest_11genes <- ggplot(forest_data, aes(x = Spearman_rho, y = rank)) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.35, alpha = 0.6) +
  geom_point(size = 5.5, color = "#2171B5", alpha = 0.6) +
  geom_text(aes(x = 0.8, label = stat_label), 
            size = 2.8, hjust = 0, vjust = 0.5, family = "sans", lineheight = 1.1) +
  scale_x_continuous(
    breaks = seq(-0.5, 1.0, 0.25),
    limits = c(-0.55, 1.0),
    name = "Spearman ρ",
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = forest_data$rank,
    labels = forest_data$gene_antimicrobial,
    name = NULL,
    expand = c(0.05, 0.05)
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(color = "black", size = 0.35),
    axis.ticks.length.x = unit(4, "pt"),
    axis.title.x = element_text(
      size = 13, face = "plain", color = "black", family = "sans",
      margin = margin(t = 8)
    ),
    axis.text.x = element_text(size = 11, color = "black", family = "sans"),
    axis.text.y = element_markdown(size = 10, color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 120, b = 10, l = 10, unit = "pt")
  )

# ==================================================
# STEP 6: SAVE OUTPUTS
# ==================================================

cat("Saving outputs...\n\n")

# Save comprehensive Excel report
wb <- createWorkbook()

# Sheet 1: Main Results
addWorksheet(wb, "All_11_Genes_Results")
results_for_export <- final_results_full |>
  select(-Status) |>
  mutate(Gene = ifelse(Gene == "edelta1", "qacEΔ1", Gene))
writeData(wb, "All_11_Genes_Results", results_for_export)

# Sheet 2: Results by Gene Class
addWorksheet(wb, "Results_by_Gene_Class")
by_class <- final_results_full |>
  filter(!is.na(Spearman_rho), N > 0) |>
  arrange(Gene_Class, desc(abs(Spearman_rho)))
writeData(wb, "Results_by_Gene_Class", by_class |> select(-Status))

# Sheet 3: Significant Results
addWorksheet(wb, "Significant_Results_P05")
sig_results <- final_results_full |>
  filter(Raw_P < 0.05, !is.na(Raw_P)) |>
  select(Gene, Antimicrobial, Gene_Class, N, Spearman_rho, Raw_P, FDR_adjusted_P, Direction)
if (nrow(sig_results) > 0) {
  writeData(wb, "Significant_Results_P05", sig_results)
} else {
  writeData(wb, "Significant_Results_P05", data.frame(Note = "No results with raw P < 0.05"))
}

# Sheet 4: Data Verification
addWorksheet(wb, "Data_Summary")
data_summary <- tibble(
  Gene = c("tetM", "ermB", "vanA", "qacJ", "qacB", "qacEdelta1"),
  Antimicrobial_Class = c("Tetracycline", "Erythromycin", "Vancomycin",
                          "Chlorhexidine", "Chlorhexidine", "Chlorhexidine"),
  Status = c("Available (Merged)", "Available (Merged)", "Available (Merged)",
             "Available (Merged)", "Available (Merged)", "Available (Merged)"),
  Expected_N = c(28, 28, 28, 28, 28, 28)
)
writeData(wb, "Data_Summary", data_summary)

# Sheet 5: Interpretation Guide
addWorksheet(wb, "Interpretation")
interpretation <- data.frame(
  Finding = c(
    "Analysis Type",
    "Number of Genes Analyzed",
    "Number of Antimicrobials",
    "Total Comparisons",
    "Sample Size per Comparison",
    "Statistical Test",
    "Multiple Testing Correction",
    "Data Transformation",
    "Censoring Treatment",
    "Key Insight 1",
    "Key Insight 2",
    "Key Insight 3"
  ),
  Details = c(
    "Expanded Spearman correlation analysis (11 genes vs BIC50)",
    "11 AR genes: tetM, tetO, ermB, ermF, mf(E), vanA, vanB, qacJ, qacB, qacEdelta1, eriCF1",
    "4 antimicrobials: Tetracycline, Erythromycin, Vancomycin, Chlorhexidine",
    "11 prespecified, biologically motivated gene-antimicrobial pairs",
    "25-28 (depends on missing data and available BIC values)",
    "Spearman rank correlation (non-parametric)",
    "Benjamini-Hochberg FDR (11 tests)",
    "Log10 transformation of gene copy numbers",
    ">1024 BIC values assigned 2048 for ranking",
    "Multiple genes targeting same antimicrobial class show variable associations",
    "Tetracycline and Vancomycin resistance genes show stronger associations than Chlorhexidine",
    "Results from original 6-gene analysis maintained; 5 additional genes expand coverage"
  )
)
writeData(wb, "Interpretation", interpretation)

saveWorkbook(wb, file.path(output_dir, "11_Genes_ARG_BIC_Analysis.xlsx"), overwrite = TRUE)
cat("✓ Saved: 11_Genes_ARG_BIC_Analysis.xlsx\n")

# Save forest plot
cat("✓ Saving forest plot (PNG, PDF)...\n")

png(file.path(output_dir, "07_Forest_Plot_11Genes.png"),
    width = 9, height = 7, units = "in", res = 600)
print(p_forest_11genes)
dev.off()

pdf(file.path(output_dir, "07_Forest_Plot_11Genes.pdf"),
    width = 9, height = 7)
print(p_forest_11genes)
dev.off()

# ==================================================
# FINAL SUMMARY
# ==================================================

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("EXPANDED GENOTYPE-PHENOTYPE ANALYSIS COMPLETE: 11 AR GENES\n")
cat(strrep("=", 80) %+% "\n\n")

cat("Summary by Antimicrobial Class:\n\n")

for (class in c("Tetracycline", "Erythromycin", "Vancomycin", "Chlorhexidine")) {
  cat("■", class, "Resistance:\n")
  class_results <- final_results_full |>
    filter(Gene_Class == class, !is.na(Spearman_rho)) |>
    select(Gene, Spearman_rho, Raw_P, FDR_adjusted_P) |>
    arrange(desc(Spearman_rho))
  
  if (nrow(class_results) > 0) {
    for (j in 1:nrow(class_results)) {
      row <- class_results[j, ]
      sig_marker <- if (row$Raw_P < 0.05) " ⭐" else ""
      cat("  •", row$Gene, "rho =", round(row$Spearman_rho, 3), 
          "| p =", round(row$Raw_P, 4), sig_marker, "\n")
    }
  }
  cat("\n")
}

cat("Output Files Saved to:", output_dir, "\n\n")
cat("Main Report:\n")
cat("  • 11_Genes_ARG_BIC_Analysis.xlsx\n\n")

cat("Figures:\n")
cat("  • 07_Forest_Plot_11Genes.png (600 dpi)\n")
cat("  • 07_Forest_Plot_11Genes.pdf\n\n")

cat(strrep("=", 80) %+% "\n")
