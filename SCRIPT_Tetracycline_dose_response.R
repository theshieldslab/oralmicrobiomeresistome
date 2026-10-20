
# =====================================================
# Tetracycline DOSE-RESPONSE (FINAL VERSION - with caries risk shapes)
# =====================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(writexl)

main_folder    <- "~/Desktop"
output_folder  <- "~/Desktop/CV may 28 patients"
metadata_file  <- file.path(output_folder, "16S_Master_Summary_all_18_samples_FINAL.xlsx")
input_file     <- file.path(main_folder, "pos doc results _all the tests.xlsx")

# ── Caries risk metadata (rows 1–28 = Patient 1–28) ──────────────────────────
caries_risk_meta <- read_excel(metadata_file, sheet = "caries risk") |>
  slice(1:28) |>
  mutate(
    Patient     = factor(paste("Patient", 1:28), levels = paste("Patient", 1:28)),
    caries_risk = `caries risk`
  ) |>
  select(Patient, caries_risk)

# ── Raw data ──────────────────────────────────────────────────────────────────
raw_data <- read_excel(input_file, sheet = "BIC_for R", col_names = FALSE)

tet_raw <- raw_data |>
  select(
    CP = 40, `2` = 41, `4` = 42, `8` = 43, `16` = 44, `32` = 45,
    `64` = 46, `128` = 47, `256` = 48, `512` = 49, `1024` = 50,
    Patient = 51
  ) |>
  filter(row_number() > 1) |>
  fill(Patient, .direction = "down") |>
  filter(!is.na(Patient)) |>
  filter(!Patient %in% c("Tetracycline", "Chx", "tet"))

tet_long <- tet_raw |>
  mutate(across(CP:`1024`, as.character)) |>
  pivot_longer(cols = CP:`1024`, names_to = "Concentration_label", values_to = "OD_percent") |>
  mutate(
    OD_percent          = as.numeric(str_remove(OD_percent, "%")),
    OD_percent          = ifelse(OD_percent < 2, OD_percent * 100, OD_percent),
    OD_percent          = ifelse(Concentration_label == "CP", 100, OD_percent),
    Patient             = str_to_title(Patient),
    Patient             = factor(Patient, levels = paste("Patient", 1:28)),
    Concentration_label = factor(
      Concentration_label,
      levels = c("CP", "2", "4", "8", "16", "32", "64", "128", "256", "512", "1024")
    )
  )

tet_summary <- tet_long |>
  group_by(Patient, Concentration_label) |>
  summarise(
    mean_OD = mean(OD_percent, na.rm = TRUE),
    se_OD   = sd(OD_percent, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    mean_OD = ifelse(Concentration_label == "CP", 100, mean_OD),
    se_OD   = ifelse(Concentration_label == "CP", 0,   se_OD)
  ) |>
  left_join(caries_risk_meta, by = "Patient")

# ── Plot ──────────────────────────────────────────────────────────────────────
p_tet <- ggplot(
  tet_summary,
  aes(x = Concentration_label, y = mean_OD,
      group = Patient, color = Patient, shape = caries_risk)
) +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  geom_point(size = 1.8) +
  geom_errorbar(
    aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD),
    width = 0.08, linewidth = 0.25, alpha = 0.8
  ) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.5, color = "black") +
  scale_shape_manual(
    values = c("high" = 16, "low" = 15),
    labels = c("high" = "High risk", "low" = "Low risk")
  ) +
  coord_cartesian(ylim = c(0, 105)) +
  labs(
    title = "Tetracycline dose-response",
    x     = "Concentration (µg/mL)",
    y     = "Biomass (%)",
    color = "Patient",
    shape = "Caries risk"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x  = element_text(size = 12, face = "bold"),
    axis.title.y  = element_text(size = 12, face = "bold"),
    axis.text.x   = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y   = element_text(size = 10),
    legend.title  = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.text   = element_text(size = 9)
  ) +
  guides(
    color = guide_legend(ncol = 2),
    shape = guide_legend(override.aes = list(size = 3, color = "black"))
  )

print(p_tet)

ggsave(
  file.path(output_folder, "Tetracycline_dose_response_FINAL.png"),
  plot = p_tet, width = 10, height = 8, dpi = 600
)

# ── BIC interpretation ────────────────────────────────────────────────────────
tet_interpretation <- tet_summary |>
  mutate(
    conc_num = case_when(
      Concentration_label == "CP" ~ 0,
      TRUE ~ as.numeric(as.character(Concentration_label))
    )
  ) |>
  group_by(Patient) |>
  summarise(
    CP_biomass          = mean_OD[Concentration_label == "CP"],
    threshold_50        = CP_biomass / 2,
    first_conc_below_50 = {
      below <- conc_num[mean_OD <= threshold_50 & Concentration_label != "CP"]
      if (length(below) == 0) NA else min(below)
    },
    BIC = case_when(
      is.na(first_conc_below_50) ~ ">1024",
      TRUE ~ as.character(first_conc_below_50)
    ),
    .groups = "drop"
  )

print(tet_interpretation, n = 28)

write_xlsx(tet_interpretation, file.path(output_folder, "Tetracycline_BIC_interpretation.xlsx"))

cat("\nSaved Tetracycline figure and interpretation.\n")
