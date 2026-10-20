# ==========================================================
# SPEARMAN CORRELATION: ERYTHROMYCIN BIOMASS vs CFU
# 28 PATIENTS - STANDARDIZED STYLING
# Matching AUC Response Plot Parameters
# ==========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(stringr)

# =============================================
# 1. SETUP
# =============================================

main_folder <- "/Users/macbookair/Desktop"
output_dir <- file.path(main_folder, "Biomass_CFU_results_28_patients")
cfu_file <- file.path(output_dir, "CFU_results_all_antimicrobials.xlsx")

# =============================================
# 2. LOAD DATA
# =============================================

# Load merged data from existing Excel
merged_all <- read_excel(file.path(output_dir, "Spearman_Biomass_vs_CFU_28patients.xlsx"), sheet = "Merged_Data")

# Filter for Erythromycin only
data_ery <- merged_all %>%
  filter(Antimicrobial == "Erythromycin") %>%
  mutate(Patient = as.character(Patient))

# =============================================
# 3. CALCULATE SPEARMAN
# =============================================

spearman_ery <- cor.test(data_ery$Biomass, data_ery$CFU, method = "spearman", exact = FALSE)
rho_ery <- round(as.numeric(spearman_ery$estimate), 2)
p_ery <- signif(spearman_ery$p.value, 2)

# =============================================
# 4. IDENTIFY EXTREME POINTS
# =============================================

model <- lm(CFU ~ Biomass, data = data_ery)

data_ery <- data_ery %>%
  mutate(
    residuals = resid(model),
    extreme_score = abs(residuals)
  )

# Get top 5 extreme patients
outlier_patients <- data_ery %>%
  group_by(Patient) %>%
  summarise(max_extreme_score = max(extreme_score, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_extreme_score)) %>%
  slice_head(n = 5) %>%
  pull(Patient)

label_data <- data_ery %>%
  filter(Patient %in% outlier_patients) %>%
  group_by(Patient) %>%
  slice_max(order_by = extreme_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# =============================================
# 5. CREATE PLOT - STANDARDIZED STYLING
# =============================================

p_ery_spearman <- ggplot(data_ery, aes(x = Biomass, y = CFU)) +
  geom_point(
    aes(fill = Group),
    size = 3.8,
    alpha = 0.6,
    shape = 21,
    stroke = 0.4,
    color = "black"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.7,
    color = "black"
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = Patient),
    size = 4.5,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.color = "grey70",
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = paste0("Spearman ρ = ", rho_ery, "\np = ", p_ery),
    hjust = -0.05,
    vjust = 1.2,
    size = 4.5,
    fontface = "plain"
  ) +
  scale_fill_manual(
    values = c(
      "CT" = "#4C78A8",
      "BIC" = "#54A24B",
      "2x BIC" = "#E45756"
    ),
    name = "Group"
  ) +
  coord_cartesian(
    xlim = c(0.3, 1.05),
    ylim = c(0, 7)
  ) +
  labs(
    title = "Erythromycin Biomass-CFU correlation",
    x = "Biofilm biomass",
    y = expression(bold(Log[10]~"(CFU/mL)"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.margin = margin(t = 60, r = 10, b = 10, l = 10),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold", angle = 90),
    axis.text = element_text(size = 16),
    axis.line = element_line(linewidth = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),
    legend.position = c(0.98, 0.05),
    legend.justification = c(1, 0),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.key = element_rect(fill = NA, color = NA),
    legend.key.size = unit(0.5, "cm")
  )

p_ery_spearman

# =============================================
# 6. SAVE
# =============================================

# Save as PDF
ggsave(
  filename = file.path(output_dir, "Erythromycin_Biomass_vs_CFU_STANDARDIZED.pdf"),
  plot = p_ery_spearman,
  width = 8,
  height = 6,
  dpi = 600
)

# Save as PNG
ggsave(
  filename = file.path(output_dir, "Erythromycin_Biomass_vs_CFU_STANDARDIZED.png"),
  plot = p_ery_spearman,
  width = 8,
  height = 6,
  dpi = 600
)

cat("✓ Erythromycin Spearman plot saved (PDF + PNG)!\n")
