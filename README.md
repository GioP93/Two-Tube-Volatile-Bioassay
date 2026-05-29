# Two-Tube Volatile Bioassay — data and analysis code

This repository contains the underlying data and the analysis scripts for the
two-tube volatile (vapour) bioassay study across three sites (AIRID, KEMRI, LSTM)
and a benchtop bottle bioassay.

## Data files

- **Two_Tube_Time_Course_Data.xlsx** — one row per tube per time point.
  Columns: `Date`, `Site`, `Strain`, `Treatment`, `Dose (mg)`, `Replicate`,
  `Tube`, `KD`, `Mort 24hrs`, `N tested`, `Time (mins)`.
- **Bottle_Bioassay_Data.xlsx** — one row per replicate bottle.
  Columns: `Site`, `Strain`, `Insecticide`, `Concentration (ug/bottle)`,
  `Replicate`, `60 min kd`, `24 hr mort`, `total`, `60 min kd %`, `24 hr mort %`.

Both files are analysis-ready: column names match the scripts exactly and no
reshaping is required before running them.

## Scripts (R)

Each script reads the relevant Excel file from the working directory and writes
its output there. Set the working directory to the folder containing the data
and scripts (e.g. `setwd(...)`), then run.

| Script | Produces |
| --- | --- |
| `01_EC50_metrics.R` | `EC50_RR_results_table.csv` — EC50, SE, 95% CI and resistance ratio (vs Kisumu) for 60-min knockdown and 24-h mortality, per site and strain |
| `02_mortality_dose_response_plot.R` | `24h_mortality_dose_response.png` |
| `03_knockdown_over_time_plot.R` | `knockdown_over_time.png` |
| `04_bottle_bioassay_boxplots.R` | `bottle_bioassay_boxplots.png` |

### Required R packages

`readxl`, `dplyr`, `stringr`, `ggplot2`, `tidyr`, and `drc` (metrics only).

```r
install.packages(c("readxl", "dplyr", "stringr", "ggplot2", "tidyr", "drc"))
```

## Method notes

- EC50 estimates use a log-logistic two-parameter model (`drc::LL.2`, binomial),
  fitted jointly across strains within a site (`curveid = Strain`), excluding the
  0 mg control (the log-based model cannot use a zero dose).
- EC50, SE and 95% CI are from `drc::ED` (delta method); resistance ratios from
  `drc::EDcomp` (delta method).
- A strain is reported only if every replicate reaches the 50% effect threshold
  at the top dose; otherwise it is excluded from analysis.
