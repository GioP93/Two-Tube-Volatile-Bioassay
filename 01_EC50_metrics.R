# =============================================================================
# Two-Tube Volatile Bioassay - EC50 / Resistance Ratio metrics
# -----------------------------------------------------------------------------
# Produces a single results table:
#   Endpoint | Site | Strain | EC50 (mg) | SE | Lower 95CI | Upper 95CI | RR
#
# Endpoints:
#   - Knockdown 60 mins   (dose-response of 60-min knockdown)
#   - Mortality 24 hours  (dose-response of 24-h mortality)
#
# Method:
#   - Per Site, a single joint log-logistic LL.2 model (binomial) is fitted to
#     all strains at once (curveid = Strain), excluding the 0 mg control
#     (LL.2 is log-based and cannot take a zero dose).
#   - EC50, SE and 95% CI come from drc::ED (delta method).
#   - Resistance Ratio (RR) vs Kisumu comes from drc::EDcomp (delta method).
#   - A strain is only reported if every replicate reaches the 50% effect
#     threshold at the top dose; otherwise it is excluded.
#
# Input : Two_Tube_Time_Course_Data.xlsx  (single sheet, already analysis-ready)
#         Columns: Date, Site, Strain, Treatment, Dose (mg), Replicate, Tube,
#                  KD, Mort 24hrs, N tested, Time (mins)
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(drc)
})

# drc loads MASS, whose select()/filter() mask dplyr's. Pin them back to dplyr.
select <- dplyr::select
filter <- dplyr::filter

input_file     <- "Two_Tube_Time_Course_Data.xlsx"
ref_strain     <- "Kisumu"          # reference strain for resistance ratios
ed50_threshold <- 50                # % effect required at top dose to report

# =============================================================================
# READ AND CLEAN DATA
# =============================================================================

dat <- read_excel(input_file)
names(dat) <- str_squish(names(dat))

dat <- dat %>%
  rename(
    Dose_mg  = `Dose (mg)`,
    KD_count = `KD`,
    Mort_count = `Mort 24hrs`,
    N_tested = `N tested`,
    Time_min = `Time (mins)`
  ) %>%
  mutate(
    Site       = str_squish(as.character(Site)),
    Strain     = str_squish(as.character(Strain)),
    Dose_mg    = suppressWarnings(as.numeric(Dose_mg)),
    KD_count   = suppressWarnings(as.numeric(KD_count)),
    Mort_count = suppressWarnings(as.numeric(Mort_count)),
    N_tested   = suppressWarnings(as.numeric(N_tested)),
    Time_min   = suppressWarnings(as.numeric(Time_min))
  ) %>%
  filter(!is.na(Site), !is.na(Strain), Site != "", Strain != "")

# Knockdown is read at the 60-minute time point; mortality is one value per tube
kd60_dat <- dat %>%
  filter(Time_min == 60, !is.na(KD_count), !is.na(N_tested), N_tested > 0)

mort_dat <- dat %>%
  filter(!is.na(Mort_count), !is.na(N_tested), N_tested > 0) %>%
  distinct(Date, Site, Strain, Dose_mg, Replicate, Tube,
           Mort_count, N_tested)

# =============================================================================
# MODEL HELPERS
# =============================================================================

group_to_df <- function(g) {
  g <- as.data.frame(g, stringsAsFactors = FALSE)
  g[] <- lapply(g, function(x) if (is.factor(x)) as.character(x) else x)
  g
}

get_curve_names <- function(fit) {
  cid <- fit$dataList$curveid
  if (!is.null(cid) && length(cid) > 0) {
    lvls <- if (is.factor(cid)) levels(droplevels(cid)) else unique(as.character(cid))
    if (length(lvls) > 0) return(lvls)
  }
  unique(sub("^[a-zA-Z]+:", "", names(coef(fit))))
}

# Per-strain check: did every replicate reach the threshold at the top dose?
response_threshold_status <- function(df, x_var, resp_count, resp_total,
                                      threshold = ed50_threshold) {
  valid <- df %>%
    filter(
      !is.na(.data[[resp_count]]), !is.na(.data[[resp_total]]), !is.na(.data[[x_var]]),
      .data[[resp_total]] > 0, .data[[x_var]] > 0,
      .data[[resp_count]] >= 0, .data[[resp_count]] <= .data[[resp_total]]
    )
  if (nrow(valid) == 0) {
    return(data.frame(Strain = character(), threshold_reached = logical(),
                      min_x = numeric(), max_x = numeric()))
  }
  id_cols <- intersect(c("Strain", "Replicate", "Tube"), names(valid))
  
  by_unit <- valid %>%
    group_by(across(all_of(id_cols))) %>%
    filter(.data[[x_var]] == max(.data[[x_var]], na.rm = TRUE)) %>%
    summarise(
      response_at_max_pct = sum(.data[[resp_count]], na.rm = TRUE) /
        sum(.data[[resp_total]], na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  ranges <- valid %>%
    group_by(Strain) %>%
    summarise(min_x = min(.data[[x_var]], na.rm = TRUE),
              max_x = max(.data[[x_var]], na.rm = TRUE), .groups = "drop")
  
  by_unit %>%
    group_by(Strain) %>%
    summarise(threshold_reached = all(response_at_max_pct >= threshold, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(ranges, by = "Strain") %>%
    mutate(Strain = as.character(Strain))
}

# Joint LL.2 (binomial) fit across all strains at one site
fit_joint <- function(df, resp_count, resp_total, x_var, min_levels = 3) {
  df <- df %>%
    filter(
      !is.na(.data[[resp_count]]), !is.na(.data[[resp_total]]), !is.na(.data[[x_var]]),
      .data[[resp_total]] > 0, .data[[x_var]] > 0,
      .data[[resp_count]] >= 0, .data[[resp_count]] <= .data[[resp_total]]
    ) %>%
    mutate(Strain = droplevels(factor(Strain)))
  
  keep <- df %>%
    group_by(Strain) %>%
    summarise(n_lev = n_distinct(.data[[x_var]]), .groups = "drop") %>%
    filter(n_lev >= min_levels) %>%
    pull(Strain)
  
  df <- df %>% filter(Strain %in% keep) %>% mutate(Strain = droplevels(factor(Strain)))
  if (n_distinct(df$Strain) == 0) return(NULL)
  
  tryCatch(
    suppressWarnings(drm(
      as.formula(paste0(resp_count, " / ", resp_total, " ~ ", x_var)),
      curveid = Strain, weights = df[[resp_total]], data = df,
      fct = LL.2(), type = "binomial"
    )),
    error = function(e) NULL
  )
}

clean_names <- function(rn, curve_names, pct = 50) {
  out <- gsub("\\s+", "", rn)
  out <- gsub(paste0(":", pct, "$"), "", out)
  out <- gsub("^e:", "", out)
  out <- gsub("^ED", "", out)
  if (all(out %in% curve_names)) return(out)
  matched <- vapply(out, function(x) {
    hits <- curve_names[vapply(curve_names, function(cn) grepl(cn, x, fixed = TRUE), logical(1))]
    if (length(hits) == 1) hits else NA_character_
  }, character(1))
  if (all(!is.na(matched))) return(matched)
  curve_names
}

extract_ed50 <- function(fit, pct = 50) {
  if (is.null(fit)) return(NULL)
  ed <- tryCatch(ED(fit, pct, interval = "delta", type = "relative", display = FALSE),
                 error = function(e) NULL)
  if (is.null(ed) || length(ed) == 0) return(NULL)
  if (!is.matrix(ed)) ed <- matrix(ed, nrow = 1)
  cn <- get_curve_names(fit)
  nc <- ncol(ed)
  data.frame(
    Strain           = clean_names(rownames(ed), cn, pct),
    Model_Estimate   = ed[, 1],
    Model_SE         = ed[, 2],
    Model_Lower_95CI = if (nc >= 4) ed[, nc - 1] else NA_real_,
    Model_Upper_95CI = if (nc >= 4) ed[, nc] else NA_real_,
    row.names = NULL, stringsAsFactors = FALSE
  )
}

extract_rr <- function(fit, ref = ref_strain, pct = 50) {
  if (is.null(fit)) return(NULL)
  cn <- get_curve_names(fit)
  if (!(ref %in% cn) || length(setdiff(cn, ref)) == 0) return(NULL)
  
  rr <- tryCatch(EDcomp(fit, c(pct, pct), interval = "delta", display = FALSE),
                 error = function(e) NULL)
  if (is.null(rr) || length(rr) == 0) return(NULL)
  if (!is.matrix(rr)) rr <- matrix(rr, nrow = 1, dimnames = list(names(rr), NULL))
  
  comp  <- rownames(rr)
  clean <- gsub("\\s+", "", comp)
  clean <- gsub("e:", "", clean, fixed = TRUE)
  clean <- gsub(paste0(":", pct), "", clean, fixed = TRUE)
  clean <- gsub("\\(|\\)", "", clean)
  parts <- strsplit(clean, "/", fixed = TRUE)
  num   <- vapply(parts, function(x) x[1], character(1))
  den   <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  nc    <- ncol(rr)
  
  data.frame(
    Numerator = num, Denominator = den,
    Raw_RR = rr[, 1], Raw_Lower = if (nc >= 4) rr[, nc - 1] else NA_real_,
    Raw_Upper = if (nc >= 4) rr[, nc] else NA_real_, row.names = NULL,
    stringsAsFactors = FALSE
  ) %>%
    filter(Numerator == ref | Denominator == ref) %>%
    mutate(
      inverted    = Numerator == ref,
      Test_Strain = ifelse(inverted, Denominator, Numerator),
      RR          = ifelse(inverted, 1 / Raw_RR, Raw_RR)
    ) %>%
    select(Test_Strain, RR)
}

# =============================================================================
# RUN ANALYSIS PER ENDPOINT AND SITE
# =============================================================================

run_endpoint <- function(data, endpoint_label, resp_count) {
  sites <- unique(data$Site)
  out <- list()
  
  for (st in sites) {
    sub <- data %>% filter(Site == st, Dose_mg > 0)
    status <- response_threshold_status(sub, "Dose_mg", resp_count, "N_tested")
    fit    <- fit_joint(sub, resp_count, "N_tested", "Dose_mg")
    
    ed <- extract_ed50(fit)
    rr <- extract_rr(fit)
    if (is.null(ed)) next
    
    ed <- ed %>%
      left_join(status, by = "Strain") %>%
      mutate(
        threshold_reached = ifelse(is.na(threshold_reached), FALSE, threshold_reached),
        Reportable = threshold_reached & is.finite(Model_Estimate) &
          (is.na(max_x) | Model_Estimate <= max_x)
      )
    
    if (!is.null(rr)) ed <- ed %>% left_join(rr, by = c("Strain" = "Test_Strain"))
    if (!("RR" %in% names(ed))) ed$RR <- NA_real_
    
    ed <- ed %>%
      mutate(
        Endpoint = endpoint_label,
        Site     = st,
        EC50     = ifelse(Reportable, round(Model_Estimate, 3), NA_real_),
        SE       = ifelse(Reportable, round(Model_SE, 3), NA_real_),
        Lower_95CI = ifelse(Reportable, round(Model_Lower_95CI, 3), NA_real_),
        Upper_95CI = ifelse(Reportable, round(Model_Upper_95CI, 3), NA_real_),
        RR       = ifelse(Reportable & Strain != ref_strain, round(RR, 2), NA_real_),
        Note     = ifelse(Reportable, "",
                          "Excluded from analysis. 50% effect threshold not reached")
      ) %>%
      select(Endpoint, Site, Strain, EC50, SE, Lower_95CI, Upper_95CI, RR, Note)
    
    out[[st]] <- ed
  }
  bind_rows(out)
}

kd_table   <- run_endpoint(kd60_dat, "Knockdown 60 mins",  "KD_count")
mort_table <- run_endpoint(mort_dat, "Mortality 24 hours", "Mort_count")

results <- bind_rows(kd_table, mort_table)

# =============================================================================
# OUTPUT
# =============================================================================

write.csv(results, "EC50_RR_results_table.csv", row.names = FALSE)

print(results, row.names = FALSE)
message("\nSaved: EC50_RR_results_table.csv")