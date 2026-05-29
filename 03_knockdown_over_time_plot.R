# =============================================================================
# Two-Tube Volatile Bioassay - knockdown over time figure
# -----------------------------------------------------------------------------
# Produces one figure: knockdown (%) vs time post exposure.
#   Rows    = sites,  Columns = doses
#   Thin lines = individual replicate tubes
#   Bold lines = strain mean per site
#
# Input : Two_Tube_Time_Course_Data.xlsx  (single sheet, analysis-ready)
#         Columns: Date, Site, Strain, Treatment, Dose (mg), Replicate, Tube,
#                  KD, Mort 24hrs, N tested, Time (mins)
# Output: knockdown_over_time.png
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

input_file  <- "Two_Tube_Time_Course_Data.xlsx"
output_file <- "knockdown_over_time.png"

# Strain plotting order and colours
strain_order <- c("Covè", "Cove", "Tiassalé", "Tiassale", "KDR", "Siaya", "Kisumu")
strain_cols  <- c(
  "Kisumu"   = "#6EC9FF",
  "KDR"      = "#9EE37D",
  "Tiassale" = "#FFB347", "Tiassalé" = "#FFB347",
  "Cove"     = "#B266FF", "Covè"     = "#B266FF",
  "Siaya"    = "#FF6B6B"
)

dose_order_mg   <- c(0, 0.005, 0.1, 2)
time_breaks     <- c(0, 5, 10, 15, 30, 45, 60)
replicate_alpha <- 0.65

# =============================================================================
# READ AND CLEAN DATA
# =============================================================================

raw <- read_excel(input_file)
names(raw) <- str_squish(names(raw))

dat <- raw %>%
  rename(
    Dose_mg  = `Dose (mg)`,
    KD_raw   = `KD`,
    N_tested = `N tested`,
    Time_min = `Time (mins)`
  ) %>%
  mutate(
    Site      = str_squish(as.character(Site)),
    Strain    = str_squish(as.character(Strain)),
    Replicate = str_squish(as.character(Replicate)),
    Tube      = str_squish(as.character(Tube)),
    Dose_mg   = suppressWarnings(as.numeric(Dose_mg)),
    N_tested  = suppressWarnings(as.numeric(N_tested)),
    Time_min  = suppressWarnings(as.numeric(Time_min)),
    KD_raw    = suppressWarnings(as.numeric(KD_raw)),
    KD_pct    = ifelse(!is.na(KD_raw) & N_tested > 0, KD_raw / N_tested * 100, NA_real_)
  ) %>%
  filter(!is.na(Site), !is.na(Strain), !is.na(Dose_mg),
         Site != "", Strain != "")

# =============================================================================
# FACTORS AND COLOURS
# =============================================================================

site_levels   <- unique(dat$Site)
strain_levels <- unique(dat$Strain)
dose_levels   <- c(dose_order_mg, setdiff(sort(unique(dat$Dose_mg)), dose_order_mg))

dat <- dat %>%
  mutate(
    Site       = factor(Site, levels = site_levels),
    Strain     = factor(Strain, levels = strain_levels),
    Dose_panel = factor(Dose_mg, levels = dose_levels, labels = as.character(dose_levels))
  )

plot_colours <- strain_cols[strain_levels]

# =============================================================================
# REPLICATE AND MEAN DATA FRAMES
# =============================================================================

kd_dat <- dat %>%
  filter(!is.na(Time_min), !is.na(KD_pct), !is.na(N_tested), N_tested > 0)

kd_rep_df <- kd_dat %>%
  group_by(Site, Strain, Dose_mg, Dose_panel, Replicate, Tube, Time_min) %>%
  summarise(rep_kd = mean(KD_pct, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site, Strain, Dose_mg, Replicate, Tube) %>%
  mutate(rep_id = cur_group_id()) %>%
  ungroup() %>%
  mutate(SiteStrain = paste(Site, Strain, sep = " - "))

kd_mean_df <- kd_rep_df %>%
  group_by(Site, Strain, Dose_mg, Dose_panel, Time_min) %>%
  summarise(mean_kd = mean(rep_kd, na.rm = TRUE), .groups = "drop") %>%
  mutate(SiteStrain = paste(Site, Strain, sep = " - "))

# =============================================================================
# PLOT
# =============================================================================

base_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    strip.background  = element_rect(fill = "grey92", colour = "grey60"),
    strip.text        = element_text(face = "plain"),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    legend.text       = element_text(size = 9),
    legend.key.height = unit(0.55, "lines")
  )

pct_y <- scale_y_continuous(
  limits = c(0, 100), breaks = seq(0, 100, 25),
  labels = function(x) paste0(x, "%"),
  expand = expansion(mult = c(0.01, 0.03))
)

p <- ggplot(mapping = aes()) +
  geom_line(
    data = kd_rep_df,
    aes(x = Time_min, y = rep_kd, colour = Strain,
        group = interaction(rep_id, SiteStrain)),
    linewidth = 0.35, alpha = replicate_alpha
  ) +
  geom_line(
    data = kd_mean_df,
    aes(x = Time_min, y = mean_kd, colour = Strain,
        group = interaction(Dose_panel, SiteStrain)),
    linewidth = 1.4
  ) +
  geom_point(
    data = kd_mean_df,
    aes(x = Time_min, y = mean_kd, colour = Strain),
    size = 2.5, shape = 19
  ) +
  facet_grid(Site ~ Dose_panel, drop = FALSE) +
  scale_colour_manual(values = plot_colours, drop = FALSE) +
  scale_x_continuous(breaks = time_breaks) +
  pct_y +
  labs(x = "Time post exposure (min)", y = "Knockdown (%)", colour = "Strain") +
  base_theme

n_doses <- length(dose_levels)
n_sites <- length(site_levels)
ggsave(output_file, p,
       width = max(10, n_doses * 3), height = max(6, n_sites * 3),
       units = "in", dpi = 300)
message("Saved: ", output_file)
