# =============================================================================
# Bottle Bioassay - response boxplots by site
# -----------------------------------------------------------------------------
# Produces one figure of boxplots:
#   Rows    = endpoint (60-min knockdown, 24-h mortality) x insecticide
#   Columns = site
#   Boxes   = replicate variation per colony; points = individual replicates
#
# Input : Bottle_Bioassay_Data.xlsx  (single sheet, analysis-ready)
#         Columns: Site, Strain, Insecticide, Concentration (ug/bottle),
#                  Replicate, 60 min kd, 24 hr mort, total,
#                  60 min kd %, 24 hr mort %
# Output: bottle_bioassay_boxplots.png
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

input_file  <- "Bottle_Bioassay_Data.xlsx"
output_file <- "bottle_bioassay_boxplots.png"

# Colony order and colours
strain_order <- c("Covè", "Cove", "Tiassalé", "Tiassale", "KDR", "Siaya", "Kisumu")
colony_cols  <- c(
  "Kisumu"   = "#6EC9FF",
  "KDR"      = "#9EE37D",
  "Tiassale" = "#FFB347", "Tiassalé" = "#FFB347",
  "Cove"     = "#B266FF", "Covè"     = "#B266FF",
  "Siaya"    = "#FF6B6B"
)

# =============================================================================
# READ AND CLEAN DATA
# =============================================================================

df <- read_excel(input_file)
colnames(df) <- str_squish(colnames(df))

df <- df %>%
  rename(
    Concentration = `Concentration (ug/bottle)`,
    KD_pct        = `60 min kd %`,
    Mort_pct      = `24 hr mort %`
  ) %>%
  mutate(
    Site          = str_squish(as.character(Site)),
    Strain        = str_squish(as.character(Strain)),
    Insecticide   = str_squish(as.character(Insecticide)),
    Concentration = as.numeric(Concentration),
    KD_pct        = as.numeric(KD_pct),
    Mort_pct      = as.numeric(Mort_pct),
    # Percentages stored as proportions are rescaled to 0-100
    KD_pct   = ifelse(!is.na(KD_pct)   & KD_pct   <= 1, KD_pct   * 100, KD_pct),
    Mort_pct = ifelse(!is.na(Mort_pct) & Mort_pct <= 1, Mort_pct * 100, Mort_pct),
    # Treatment axis: Control vs labelled dose
    Concentration_f = factor(
      ifelse(Concentration == 0, "Control", paste0(Concentration, " \u00b5g")),
      levels = c("Control", paste0(sort(unique(Concentration[Concentration > 0])), " \u00b5g"))
    ),
    Strain = factor(Strain, levels = strain_order)
  )

# =============================================================================
# LONG FORMAT FOR COMBINED PLOTTING
# =============================================================================

plot_df <- df %>%
  pivot_longer(cols = c(KD_pct, Mort_pct), names_to = "Endpoint", values_to = "Value") %>%
  mutate(
    Endpoint = recode(Endpoint,
                      KD_pct   = "60 min knockdown (%)",
                      Mort_pct = "24 h mortality (%)"),
    Endpoint = factor(Endpoint, levels = c("60 min knockdown (%)", "24 h mortality (%)"))
  ) %>%
  drop_na(Value, Concentration_f, Site, Strain)

present_strains <- levels(droplevels(plot_df$Strain))
named_present   <- colony_cols[names(colony_cols) %in% present_strains]

# =============================================================================
# PLOT
# =============================================================================

p <- ggplot(plot_df, aes(x = Concentration_f, y = Value, fill = Strain, colour = Strain)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.65,
               position = position_dodge(width = 0.8)) +
  geom_point(aes(group = Strain),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 1.3, alpha = 0.55, colour = "black") +
  facet_grid(Endpoint + Insecticide ~ Site, scales = "free_x") +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  scale_fill_manual(values = named_present, breaks = strain_order, drop = TRUE) +
  scale_colour_manual(values = named_present, breaks = strain_order, drop = TRUE) +
  labs(x = "Treatment", y = "Response (%)", fill = "Colony", colour = "Colony") +
  theme_bw(base_size = 11) +
  theme(
    strip.text         = element_text(face = "bold", size = 9),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing      = unit(0.8, "lines"),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

# =============================================================================
# SAVE
# =============================================================================

n_sites <- n_distinct(plot_df$Site)
n_rows  <- n_distinct(plot_df$Endpoint) * n_distinct(plot_df$Insecticide)
ggsave(output_file, p,
       width = max(9, 3.5 * n_sites), height = max(7, 2.6 * n_rows), dpi = 300)
message("Saved: ", output_file)
