# Economic Shocks and Civil Conflict in the Maghreb
### A Regional Replication and Adaptation of Miguel, Satyanath & Sergenti (2004)

---

## Overview

This repository contains the replication code for a panel study examining the relationship between economic shocks and civil conflict in North Africa. The analysis follows the instrumental variables (IV) framework of Miguel, Satyanath & Sergenti's "Economic Shocks and Civil Conflict: An Instrumental Variables Approach," Journal of Political Economy, 112(4): 725–753. This study is adapted for a regional sample of five North African countries over the period 1981–2008.

Rather than using rainfall as an instrument for GDP growth (as in the original paper, which is appropriate for rain-fed agricultural economies in sub-Saharan Africa), this paper instruments GDP growth with international oil price shocks, reflecting the central role of hydrocarbon rents in North Africa.

---

## Countries and Sample Period

| Country | ISO3 Code | Sample Years |
|---------|-----------|--------------|
| Algeria | DZA | 1981–2008 |
| Egypt | EGY | 1981–2008 |
| Libya | LBY | 1981–2008 |
| Morocco | MAR | 1981–2008 |
| Tunisia | TUN | 1981–2008 |

- **Panel coverage:** 1980–2008 (restricted to PRIO/UCDP conflict data availability)
- **Estimation sample:** 1981–2008 (after one-year lag for oil price growth)
- **Observations:** 5 countries × 28 years = **140 country-year observations**

---

## Data Sources

| Variable | Source |
|----------|--------|
| GDP growth, population growth | World Bank World Development Indicators (WDI), accessed via `WDI` package |
| Civil conflict (25-death threshold, Types 3 & 4) | PRIO/UCDP Armed Conflict Dataset |
| International oil prices | IMF Primary Commodity Prices (`oil_imf.csv`) |
| Oil rents as % of GDP | World Bank WDI (local file `gdp_per.csv`) |

**The code is set up to automatically download the World Bank WDI data and should run from top to bottom with everything fully automated**.

---

## Repository Structure

```
.
├── Script.R                    # Main analysis script
├── README.md                       # This file
├── Data/
│   ├── Source Data/
│   │   ├── oil_imf.csv             # IMF oil price series
│   │   └── Main Conflict Table.xls # PRIO/UCDP conflict data
│   └── Cleaned Data/
│       └── gdp_per.csv             # Oil rents as % of GDP
└── Tables/                         # Output figures and tables (generated on run)
```
---

## Identification Strategy

The original Miguel et al. (2004) paper uses **rainfall shocks** as an instrument for GDP growth in sub-Saharan Africa, where agriculture is largely rain-fed. This instrument is not appropriate for the Maghreb, where irrigation is more common and economies are heavily tied to hydrocarbon export revenues.

This adaptation uses **international oil price growth** (current and lagged) as an instrument for GDP growth, motivated by the following:

1. Oil and gas exports constitute a substantial share of GDP and government revenue across the region (especially Algeria and Libya).
2. International oil prices are determined on world markets and are plausibly exogenous to any single country's political conditions.
3. The first-stage relationship between oil price shocks and GDP growth is testable and reported.

---

## Econometric Models

The script estimates the following model sequence:

### OLS Specifications
- **Pooled OLS:** `conflict_25 ~ gdp_growth`
- **Country FE:** `conflict_25 ~ gdp_growth + country`
- **Country FE + time trend:** `conflict_25 ~ gdp_growth + country + year`

### IV-2SLS Specifications (instrument: oil price growth)
- **IV with country FE:** instruments `gdp_growth` with `oil_growth`
- **IV with country FE + time trend**
- **Lagged IV specification**
- **Miguel-style current + lagged growth IV** (with and without FE/trend)
- **IV with population growth control**

### Robustness / Appendix
- Country-specific linear time trends (OLS and IV)
- Probit baseline
- Wooldridge test for serial correlation in panels

### Standard Errors
All reported standard errors are **clustered at the country level** to account for within-country serial correlation, motivated by the Wooldridge test result reported in Appendix Table B1.

---

## Output Tables and Figures

The script produces the following outputs in the `Tables/` directory:

| File | Description |
|------|-------------|
| `conflict_table_clean.png` | Table 1: Civil conflict incidence by country |
| `oil_rents_as_share_gdp.png` | Table 2: Oil rents as % of GDP, descriptive statistics |
| `first_stage_oil_growth_table.png` | Table 3: First-stage estimates (oil shocks → GDP growth) |
| `reduced_form_oil_conflict_table.png` | Table 4: Reduced-form estimates (oil shocks → conflict) |
| `economic_growth_conflict_table.png` | Table 5: OLS and IV-2SLS main results |
| `conflict_incidence_plot.png` | Figure 1: Conflict incidence heatmap by country-year |
| `oil_plot.png` / `oil_plot_paper.png` | Figure 2: Oil rents as % of GDP over time |
| `appendix_A1_first_stage_cstrend.png` | Appendix A1: First stage with country-specific trends |
| `appendix_A2_reduced_form_cstrend.png` | Appendix A2: Reduced form with country-specific trends |
| `appendix_A3_ols_iv_cstrend.png` | Appendix A3: OLS and IV with country-specific trends |
| `appendix_B_wooldridge_test.png` | Appendix B1: Wooldridge serial correlation test |

---

## Replication Instructions

### Prerequisites

Install the following R packages before running the script:

```r
install.packages(c(
  "WDI", "AER", "dplyr", "tidyr", "readr", "readxl",
  "knitr", "sandwich", "lmtest", "ggplot2", "gt",
  "broom", "car", "plm"
))
```

### Running the Script

1. Clone or download this repository.
2. Open `Script.R` and update the `setwd()` call at the top of the script to your local working directory.
3. Run the full script. Pulling World Bank data via `WDI()` may take several minutes.

> **Note on `setwd()`:** The script currently contains a hardcoded working directory path. Update this before running.

---End Readme---
