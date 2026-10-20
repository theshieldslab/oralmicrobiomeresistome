
# =====================================================
# NaF DOSE-RESPONSE (FINAL VERSION - with caries risk colors & formatted legend)
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

naf_raw <- raw_data |>
  select(
    CP = 1, `2` = 2, `4` = 3, `8` = 4, `16` = 5, `32` = 6,
    `64` = 7, `128` = 8, `256` = 9, `512` = 10, `1024` = 11,
    Patient = 12
  ) |>
  filter(row_number() > 1) |>
  fill(Patient, .direction = "down") |>
  filter(!is.na(Patient)) |>
  filter(!Patient %in% c("NaF", "Chx", "naf"))

naf_long <- naf_raw |>
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

naf_summary <- naf_long |>
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

# ── Prepare Data for Plotting ─────────────────────────────────────────────────
# We prepare the summary data for plotting, ensuring patients are ordered correctly.
# Note: We no longer use dummy entries for the legend.

high_patients <- caries_risk_meta |> filter(caries_risk == "high") |> pull(Patient)
low_patients  <- caries_risk_meta |> filter(caries_risk == "low")  |> pull(Patient)

# Define colors for patients based on risk
high_colors <- colorRampPalette(c("#8B0000", "#FF9999"))(length(high_patients))
low_colors  <- colorRampPalette(c("#1a5c1a", "#90EE90"))(length(low_patients))

all_colors <- c(setNames(high_colors, high_patients), 
                setNames(low_colors, low_patients))

# Ensure plot_data follows the patient order
plot_data <- naf_summary |>
  mutate(Patient = factor(Patient, levels = c(high_patients, low_patients)))

# ── Plot ──────────────────────────────────────────────────────────────────────
p_naf <- ggplot(
  plot_data,
  aes(x = Concentration_label, y = mean_OD,
      group = Patient, color = Patient)
) +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  geom_point(size = 1.5) +
  geom_errorbar(
    aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD),
    width = 0.08, linewidth = 0.25, alpha = 0.8
  ) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.5, color = "black") +
  
  scale_color_manual(values = all_colors) +
  coord_cartesian(ylim = c(0, 110)) +
  labs(
    title = "NaF dose-response",
    x     = "Concentration (µg/mL)",
    y     = "Biomass (%)",
    color = "Patient"
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
  guides(color = guide_legend(ncol = 2))

print(p_naf)

ggsave(
  file.path(output_folder, "NaF_dose_response_FINAL.png"),
  plot = p_naf, width = 10, height = 8, dpi = 600
)


# ── BIC interpretation ────────────────────────────────────────────────────────
naf_interpretation <- naf_summary |>
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

print(naf_interpretation, n = 28)

write_xlsx(naf_interpretation, file.path(output_folder, "NaF_BIC_interpretation.xlsx"))

cat("\nSaved NaF figure and interpretation.\n")
