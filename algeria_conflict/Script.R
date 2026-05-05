# ============================================================
# Miguel et al. Replication/Adaptation
# North Africa Panel: Algeria, Libya, Egypt, Tunisia, Morocco
# Restricted to PRIO coverage: 1980–2008
# Estimation sample after oil-growth lag: 1981–2008
# ============================================================

# ---- 0. Load libraries ----
#If you don't have them you'll need to install them.
#Don't forget to setwd()!

library(WDI)
library(AER)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(knitr)
library(sandwich)
library(lmtest)
library(dplyr)
library(ggplot2)
library(gt)
library(broom)
library(car)

# Don't forget to setwd()!

# ---- 0.5. Load local data ----
oil_imf <- read_csv("Data/Source Data/oil_imf.csv")
conflict_data <- read_excel("Data/Source Data/Main Conflict Table.xls")
gdp_per <- read_csv("Data/Cleaned Data/gdp_per.csv")

gdp_per <- gdp_per %>%
  mutate(country = recode(country,
                          "Egypt, Arab Rep." = "Egypt"))

# ---- 1. Countries and years ----
countries <- c("DZA", "LBY", "EGY", "TUN", "MAR")

start_year <- 1980
end_year   <- 2008

# ---- 2. Pull World Bank data ----
# May take 4-5 minutes

wb_panel <- WDI(
  country = countries,
  indicator = c(
    gdp_growth = "NY.GDP.MKTP.KD.ZG",
    pop_growth = "SP.POP.GROW"
  ),
  start = start_year,
  end = end_year,
  extra = FALSE
) %>%
  transmute(
    country = iso3c,
    year = as.integer(year),
    gdp_growth = gdp_growth,
    pop_growth = pop_growth
  ) %>%
  arrange(country, year)

# ---- 3. Clean IMF oil data ----
oil_imf_clean <- oil_imf %>%
  transmute(
    year = as.integer(TIME_PERIOD),
    oil_price = OBS_VALUE
  ) %>%
  filter(year >= start_year, year <= end_year) %>%
  arrange(year) %>%
  mutate(
    oil_growth = c(NA, diff(log(oil_price)))
  )

# ---- 4. Expand oil data across all countries ----
oil_panel <- expand.grid(
  country = countries,
  year = oil_imf_clean$year
) %>%
  left_join(oil_imf_clean, by = "year") %>%
  arrange(country, year)

# ---- 5. Build PRIO/UCDP conflict dummy ----

# Handles either Year or YEAR column name
year_col <- if ("Year" %in% names(conflict_data)) "Year" else "YEAR"

# Quick check of conflict data coverage
summary(conflict_data[[year_col]])

conflict_panel <- conflict_data %>%
  filter(
    Location %in% c("Algeria", "Libya", "Egypt", "Tunisia", "Morocco"),
    Type %in% c(3, 4),
    .data[[year_col]] >= start_year,
    .data[[year_col]] <= end_year
  ) %>%
  transmute(
    country = case_when(
      Location == "Algeria" ~ "DZA",
      Location == "Libya"   ~ "LBY",
      Location == "Egypt"   ~ "EGY",
      Location == "Tunisia" ~ "TUN",
      Location == "Morocco" ~ "MAR"
    ),
    year = as.integer(.data[[year_col]]),
    conflict_25 = 1
  ) %>%
  distinct(country, year, .keep_all = TRUE)

# ---- 6. Create full country-year conflict panel ----
conflict_full_panel <- expand.grid(
  country = countries,
  year = start_year:end_year
) %>%
  left_join(conflict_panel, by = c("country", "year")) %>%
  mutate(
    conflict_25 = ifelse(is.na(conflict_25), 0, conflict_25)
  ) %>%
  arrange(country, year)

# ---- 7. Merge full IV dataset ----
panel_iv <- wb_panel %>%
  left_join(oil_panel, by = c("country", "year")) %>%
  left_join(conflict_full_panel, by = c("country", "year")) %>%
  filter(!is.na(oil_growth)) %>%
  arrange(country, year) %>%
  mutate(
    country = as.factor(country)
  )

# ---- 8. Create lagged variables ----
panel_iv <- panel_iv %>%
  group_by(country) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    gdp_growth_lag = lag(gdp_growth),
    oil_growth_lag = lag(oil_growth),
    conflict_lag = lag(conflict_25)
  ) %>%
  ungroup()

#Change country codes to country names for consistency
panel_iv <- panel_iv %>%
  mutate(
    country_label = recode(as.character(country),
                           "DZA" = "Algeria",
                           "EGY" = "Egypt",
                           "LBY" = "Libya",
                           "MAR" = "Morocco",
                           "TUN" = "Tunisia"
    )
  )
# ---- 9. Check dataset ----
head(panel_iv, 20)
table(panel_iv$country)
table(panel_iv$conflict_25)
summary(panel_iv)
colSums(is.na(panel_iv))

# Expected final sample:
# 1981–2008 = 28 years
# 5 countries x 28 years = 140 observations

# Optional: save clean dataset
write.csv(panel_iv, "Data/Cleaned Data/panel_iv_clean_1981_2008.csv", row.names = FALSE)

# ============================================================
# MODELS
# !!!NOTE!!! Clustered SEs are used. 
# ============================================================

# ---- 10. Pooled OLS: raw relationship ----
# Model
ols_pooled <- lm(
  conflict_25 ~ gdp_growth,
  data = panel_iv)

# Model fit info (R2, etc.)
summary(ols_pooled)

# ---- 11. Country fixed effects OLS ----
ols_country_fe <- lm(
  conflict_25 ~ gdp_growth + country,
  data = panel_iv
)

summary(ols_country_fe)

# ---- 12. Country FE + linear time trend OLS ----
ols_country_trend <- lm(
  conflict_25 ~ gdp_growth + country + year,
  data = panel_iv
)

summary(ols_country_trend)

# ============================================================
# IV MODELS
# ============================================================

# ---- 14. IV with country fixed effects ----
iv_panel <- ivreg(
  conflict_25 ~ gdp_growth + country |
    oil_growth + country,
  data = panel_iv
)

summary(iv_panel, diagnostics = TRUE)

# ---- 15. IV with country fixed effects + time trend ----
iv_panel_trend <- ivreg(
  conflict_25 ~ gdp_growth + country + year |
    oil_growth + country + year,
  data = panel_iv
)

summary(iv_panel_trend, diagnostics = TRUE)

# ---- 16. Lagged OLS specification ----
ols_lag <- lm(
  conflict_25 ~ gdp_growth_lag + country + year,
  data = panel_iv
)

summary(ols_lag)

# ---- 17. Lagged IV specification ----
iv_panel_lag <- ivreg(
  conflict_25 ~ gdp_growth_lag + country + year |
    oil_growth_lag + country + year,
  data = panel_iv
)

summary(iv_panel_lag, diagnostics = TRUE)

# ---- 18. Miguel-style current + lagged growth IV, no FE/trend ----
iv_growth_lag <- ivreg(
  conflict_25 ~ gdp_growth + gdp_growth_lag |
    oil_growth + oil_growth_lag,
  data = panel_iv
)

summary(iv_growth_lag, diagnostics = TRUE)

iv_growth_lag_results <- coeftest(
  iv_growth_lag,
  vcov = vcovHC(iv_growth_lag, type = "HC1")
)

print(iv_growth_lag_results)

# ---- 19. Miguel-style current + lagged growth IV ----
iv_panel_current_lag <- ivreg(
  conflict_25 ~ gdp_growth + gdp_growth_lag + country + year |
    oil_growth + oil_growth_lag + country + year,
  data = panel_iv
)

summary(iv_panel_current_lag, diagnostics = TRUE)

iv_panel_current_lag_results <- coeftest(
  iv_panel_current_lag,
  vcov = vcovHC(iv_panel_current_lag, type = "HC1")
)

print(iv_panel_current_lag_results)

# ---- 20. Optional control: population growth ----
iv_panel_pop <- ivreg(
  conflict_25 ~ gdp_growth + pop_growth + country + year |
    oil_growth + pop_growth + country + year,
  data = panel_iv
)

summary(iv_panel_pop, diagnostics = TRUE)

# ============================================================
# ---- Descriptive Statistics ----
# ============================================================
# Oil as a Percent of GDP
oil_summary <- gdp_per %>%
  group_by(country) %>%
  summarise(
    mean_oil = mean(oil_rents, na.rm = TRUE),
    sd_oil   = sd(oil_rents, na.rm = TRUE),
    min_oil  = min(oil_rents, na.rm = TRUE),
    max_oil  = max(oil_rents, na.rm = TRUE),
    n_obs    = sum(!is.na(oil_rents)),
    .groups = "drop"
  )

region_country_mean <- oil_summary %>%
  summarise(
    country = "Region (Country Means)",
    mean_oil = mean(mean_oil, na.rm = TRUE),
    sd_oil   = NA_real_,
    min_oil  = NA_real_,
    max_oil  = NA_real_,
    n_obs = n()
  )

# Country-level conflict summary:
conflict_summary <- panel_iv %>%
  group_by(country_label) %>%
  summarise(
    conflict_years = sum(conflict_25 == 1, na.rm = TRUE),
    total_years    = sum(!is.na(conflict_25)),
    conflict_rate  = mean(conflict_25, na.rm = TRUE) * 100,
    .groups = "drop"
  )

print(conflict_summary)

# Regional summary: 
region_conflict <- panel_iv %>%
  summarise(
    country_label = "Region",
    conflict_years = sum(conflict_25 == 1, na.rm = TRUE),
    total_years    = sum(!is.na(conflict_25)),
    conflict_rate  = mean(conflict_25, na.rm = TRUE) * 100
  )

final_conflict_table <- bind_rows(conflict_summary, region_conflict)

print(final_conflict_table)

# ============================================================
# ---- Create Tables for Paper ----
# Images Will Save in Tables Folder
# ============================================================
# ----Table 1: Civil Conflict Incidence ----
#Styled table:
conflict_table_plot <- final_conflict_table %>%
  gt() %>%
  tab_header(
    title = "Civil Conflict Incidence in the Maghreb",
    subtitle = "Conflict years by country, 1981–2008"
  ) %>%
  fmt_number(
    columns = conflict_rate,
    decimals = 1
  ) %>%
  fmt_number(
    columns = c(conflict_years, total_years),
    decimals = 0
  ) %>%
  cols_label(
    country_label = "Country",
    conflict_years = "Conflict Years",
    total_years = "Total Years",
    conflict_rate = "Conflict Rate (%)"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = country_label == "Region")
  ) %>%
  opt_row_striping() %>%
  tab_options(
    table.font.size = "small",
    heading.title.font.size = "medium",
    data_row.padding = px(6)
  )
#Save to Tables folder
gtsave(conflict_table_plot, "conflict_table_clean.png", path = "Tables")

# ============================================================
# ----Table 2: Oil Rents as a Share of GDP (Maghreb)----
# ============================================================
final_table <- bind_rows(oil_summary, region_country_mean)

table_plot <- final_table %>%
  gt() %>%
  
  # Title
  tab_header(
    title = "Oil Rents as a Share of GDP (Maghreb)",
    subtitle = "Descriptive Statistics by Country, 1981–2008"
  ) %>%
  
  # Format numbers
  fmt_number(
    columns = c(mean_oil, sd_oil, min_oil, max_oil),
    decimals = 2
  ) %>%
  fmt_number(
    columns = n_obs,
    decimals = 0
  ) %>%
  
  # Rename columns
  cols_label(
    country  = "Country",
    mean_oil = "Mean (%)",
    sd_oil   = "SD",
    min_oil  = "Min",
    max_oil  = "Max",
    n_obs    = "N"
  ) %>%
  
  # Align numbers nicely
  cols_align(
    align = "center",
    -country
  ) %>%
  
  # Bold region row
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = country == "Region (Country Means)"
    )
  ) %>%
  
  # Light zebra striping
  opt_row_striping() %>%
  
  # Clean theme
  tab_options(
    table.font.size = "small",
    heading.title.font.size = "medium",
    data_row.padding = px(6)
  )

gtsave(table_plot, "oil_rents_as_share_gdp.png", path = "Tables")

# ============================================================
# ---- TABLE 3: FIRST STAGE — OIL SHOCKS AND GDP GROWTH -----
# ============================================================
# ---- First-stage models ----
fs_1 <- lm(
  gdp_growth ~ oil_growth,
  data = panel_iv
)

fs_2 <- lm(
  gdp_growth ~ oil_growth + oil_growth_lag,
  data = panel_iv
)

fs_3 <- lm(
  gdp_growth ~ oil_growth + oil_growth_lag + country,
  data = panel_iv
)

fs_4 <- lm(
  gdp_growth ~ oil_growth + oil_growth_lag + country + year,
  data = panel_iv
)

# ---- Function: clustered coefficient + SE ----
get_robust_coef <- function(model, term_name, model_name) {
  
  robust <- coeftest(
    model,
    vcov = vcovCL(model, cluster = ~country)
  )
  
  if (!(term_name %in% rownames(robust))) {
    return(tibble(
      model = model_name,
      term = term_name,
      coef_display = "",
      se_display = ""
    ))
  }
  
  coef <- robust[term_name, "Estimate"]
  se   <- robust[term_name, "Std. Error"]
  p    <- robust[term_name, "Pr(>|t|)"]
  
  stars <- case_when(
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
  
  tibble(
    model = model_name,
    term = term_name,
    coef_display = paste0(sprintf("%.3f", coef), stars),
    se_display = paste0("(", sprintf("%.3f", se), ")")
  )
}

# ---- Function: clustered F-test for oil instruments ----
get_clustered_f <- function(model, terms_to_test) {
  
  available_terms <- terms_to_test[terms_to_test %in% names(coef(model))]
  
  if (length(available_terms) == 0) {
    return("---")
  }
  
  restrictions <- paste0(available_terms, " = 0")
  
  test <- car::linearHypothesis(
    model,
    restrictions,
    vcov. = vcovCL(model, cluster = model.frame(model)$country),
    test = "F"
  )
  
  sprintf("%.2f", test$F[2])
}

# ---- Extract coefficient rows ----
fs_results <- bind_rows(
  get_robust_coef(fs_1, "oil_growth",     "(1)"),
  get_robust_coef(fs_1, "oil_growth_lag", "(1)"),
  
  get_robust_coef(fs_2, "oil_growth",     "(2)"),
  get_robust_coef(fs_2, "oil_growth_lag", "(2)"),
  
  get_robust_coef(fs_3, "oil_growth",     "(3)"),
  get_robust_coef(fs_3, "oil_growth_lag", "(3)"),
  
  get_robust_coef(fs_4, "oil_growth",     "(4)"),
  get_robust_coef(fs_4, "oil_growth_lag", "(4)")
)

coef_rows <- fs_results %>%
  mutate(
    variable = case_when(
      term == "oil_growth" ~ "Oil Growth",
      term == "oil_growth_lag" ~ "Oil Growth (t-1)"
    )
  ) %>%
  select(variable, model, coef_display) %>%
  pivot_wider(names_from = model, values_from = coef_display)

se_rows <- fs_results %>%
  mutate(
    variable = case_when(
      term == "oil_growth" ~ "",
      term == "oil_growth_lag" ~ " "
    )
  ) %>%
  select(variable, model, se_display) %>%
  pivot_wider(names_from = model, values_from = se_display)

first_stage_coef_table <- bind_rows(
  coef_rows[1, ],
  se_rows[1, ],
  coef_rows[2, ],
  se_rows[2, ]
)

# ---- Model info rows ----
info_rows <- tibble(
  variable = c(
    "Country Fixed Effects",
    "Time Trend",
    "Observations",
    "R-squared",
    "Clustered F-stat"
  ),
  `(1)` = c(
    "No",
    "No",
    nobs(fs_1),
    sprintf("%.3f", summary(fs_1)$r.squared),
    get_clustered_f(fs_1, c("oil_growth"))
  ),
  `(2)` = c(
    "No",
    "No",
    nobs(fs_2),
    sprintf("%.3f", summary(fs_2)$r.squared),
    get_clustered_f(fs_2, c("oil_growth", "oil_growth_lag"))
  ),
  `(3)` = c(
    "Yes",
    "No",
    nobs(fs_3),
    sprintf("%.3f", summary(fs_3)$r.squared),
    get_clustered_f(fs_3, c("oil_growth", "oil_growth_lag"))
  ),
  `(4)` = c(
    "Yes",
    "Yes",
    nobs(fs_4),
    sprintf("%.3f", summary(fs_4)$r.squared),
    get_clustered_f(fs_4, c("oil_growth", "oil_growth_lag"))
  )
)

# ---- Final table ----
first_stage_table <- bind_rows(
  first_stage_coef_table,
  info_rows
)

# ---- Style and save ----
first_stage_table_plot <- first_stage_table %>%
  gt() %>%
  tab_header(
    title = "Oil Shocks and Economic Growth",
    subtitle = "First-stage estimates, 1981–2008"
  ) %>%
  cols_label(variable = "") %>%
  cols_align(align = "center", columns = -variable) %>%
  cols_align(align = "left", columns = variable) %>%
  cols_width(
    variable ~ px(145),
    everything() ~ px(80)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = variable %in% c("Oil Growth", "Oil Growth (t-1)"))
  ) %>%
  tab_style(
    style = cell_text(color = "gray40"),
    locations = cells_body(rows = variable %in% c("", " "))
  ) %>%
  tab_source_note(
    source_note = "Notes: Dependent variable is annual GDP growth. Standard errors clustered at the country level are reported in parentheses to account for heteroskedasticity and within-country serial correlation. *** p < 0.01, ** p < 0.05, * p < 0.10. F-statistics test the joint significance of the oil-growth instruments using country-clustered standard errors."
  ) %>%
  opt_row_striping() %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = "medium",
    data_row.padding = px(3)
  )

gtsave(
  first_stage_table_plot,
  filename = "first_stage_oil_growth_table.png",
  path = "Tables"
)

# ============================================================
# ----TABLE 4: REDUCED FORM — OIL SHOCKS AND CIVIL CONFLICT ----
# ============================================================

library(dplyr)
library(tidyr)
library(gt)
library(lmtest)
library(sandwich)

dir.create("Tables", showWarnings = FALSE)

# ---- Reduced-form models ----
rf_1 <- lm(
  conflict_25 ~ oil_growth,
  data = panel_iv
)

rf_2 <- lm(
  conflict_25 ~ oil_growth + oil_growth_lag,
  data = panel_iv
)

rf_3 <- lm(
  conflict_25 ~ oil_growth + oil_growth_lag + country,
  data = panel_iv
)

rf_4 <- lm(
  conflict_25 ~ oil_growth + oil_growth_lag + country + year,
  data = panel_iv
)

# ---- Function: clustered coefficient + SE ----
get_clustered_coef <- function(model, term_name, model_name) {
  
  clustered <- coeftest(
    model,
    vcov = vcovCL(model, cluster = model.frame(model)$country)
  )
  
  if (!(term_name %in% rownames(clustered))) {
    return(tibble(
      model = model_name,
      term = term_name,
      coef_display = "",
      se_display = ""
    ))
  }
  
  coef <- clustered[term_name, "Estimate"]
  se   <- clustered[term_name, "Std. Error"]
  p    <- clustered[term_name, "Pr(>|t|)"]
  
  stars <- case_when(
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
  
  tibble(
    model = model_name,
    term = term_name,
    coef_display = paste0(sprintf("%.3f", coef), stars),
    se_display = paste0("(", sprintf("%.3f", se), ")")
  )
}

# ---- Extract coefficient rows ----
rf_results <- bind_rows(
  get_clustered_coef(rf_1, "oil_growth",     "(1)"),
  get_clustered_coef(rf_1, "oil_growth_lag", "(1)"),
  
  get_clustered_coef(rf_2, "oil_growth",     "(2)"),
  get_clustered_coef(rf_2, "oil_growth_lag", "(2)"),
  
  get_clustered_coef(rf_3, "oil_growth",     "(3)"),
  get_clustered_coef(rf_3, "oil_growth_lag", "(3)"),
  
  get_clustered_coef(rf_4, "oil_growth",     "(4)"),
  get_clustered_coef(rf_4, "oil_growth_lag", "(4)")
)

coef_rows <- rf_results %>%
  mutate(
    variable = case_when(
      term == "oil_growth" ~ "Oil Growth",
      term == "oil_growth_lag" ~ "Oil Growth (t-1)"
    )
  ) %>%
  select(variable, model, coef_display) %>%
  pivot_wider(names_from = model, values_from = coef_display)

se_rows <- rf_results %>%
  mutate(
    variable = case_when(
      term == "oil_growth" ~ "",
      term == "oil_growth_lag" ~ " "
    )
  ) %>%
  select(variable, model, se_display) %>%
  pivot_wider(names_from = model, values_from = se_display)

reduced_form_coef_table <- bind_rows(
  coef_rows[1, ],
  se_rows[1, ],
  coef_rows[2, ],
  se_rows[2, ]
)

# ---- Model info rows ----
info_rows <- tibble(
  variable = c(
    "Country Fixed Effects",
    "Common Time Trend",
    "Observations",
    "R-squared"
  ),
  `(1)` = c(
    "No",
    "No",
    nobs(rf_1),
    sprintf("%.3f", summary(rf_1)$r.squared)
  ),
  `(2)` = c(
    "No",
    "No",
    nobs(rf_2),
    sprintf("%.3f", summary(rf_2)$r.squared)
  ),
  `(3)` = c(
    "Yes",
    "No",
    nobs(rf_3),
    sprintf("%.3f", summary(rf_3)$r.squared)
  ),
  `(4)` = c(
    "Yes",
    "Yes",
    nobs(rf_4),
    sprintf("%.3f", summary(rf_4)$r.squared)
  )
)

# ---- Final table ----
reduced_form_table <- bind_rows(
  reduced_form_coef_table,
  info_rows
)

# ---- Style and save ----
reduced_form_table_plot <- reduced_form_table %>%
  gt() %>%
  tab_header(
    title = "Oil Shocks and Civil Conflict",
    subtitle = "Reduced-form estimates, civil conflict with 25+ deaths, 1981–2008"
  ) %>%
  cols_label(variable = "") %>%
  cols_align(align = "center", columns = -variable) %>%
  cols_align(align = "left", columns = variable) %>%
  cols_width(
    variable ~ px(145),
    everything() ~ px(80)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = variable %in% c("Oil Growth", "Oil Growth (t-1)"))
  ) %>%
  tab_style(
    style = cell_text(color = "gray40"),
    locations = cells_body(rows = variable %in% c("", " "))
  ) %>%
  tab_source_note(
    source_note = "Notes: Dependent variable equals 1 for civil conflict with at least 25 battle deaths. Standard errors clustered at the country level are reported in parentheses to account for heteroskedasticity and within-country serial correlation. *** p < 0.01, ** p < 0.05, * p < 0.10."
  ) %>%
  opt_row_striping() %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = "medium",
    data_row.padding = px(3)
  )

gtsave(
  reduced_form_table_plot,
  filename = "reduced_form_oil_conflict_table.png",
  path = "Tables"
)

# ============================================================
# ---- TABLE 5: ECONOMIC GROWTH AND CIVIL CONFLICT ----
# Miguel-style OLS and IV table
# ============================================================

library(dplyr)
library(tidyr)
library(gt)
library(lmtest)
library(sandwich)
library(AER)

dir.create("Tables", showWarnings = FALSE)

# ---- Models ----

ols_growth_lag <- lm(
  conflict_25 ~ gdp_growth + gdp_growth_lag,
  data = panel_iv
)

ols_growth_lag_trend <- lm(
  conflict_25 ~ gdp_growth + gdp_growth_lag + year,
  data = panel_iv
)

ols_growth_lag_fe_trend <- lm(
  conflict_25 ~ gdp_growth + gdp_growth_lag + country + year,
  data = panel_iv
)

iv_growth_lag_fe_trend <- iv_panel_current_lag

# ---- Function: clustered coefficient + SE ----
get_clustered_coef <- function(model, term_name, model_name) {
  
  clustered <- coeftest(
    model,
    vcov = vcovCL(model, cluster = model.frame(model)$country)
  )
  
  coef <- clustered[term_name, "Estimate"]
  se   <- clustered[term_name, "Std. Error"]
  p    <- clustered[term_name, "Pr(>|t|)"]
  
  stars <- case_when(
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
  
  tibble(
    model = model_name,
    term = term_name,
    coef_display = paste0(sprintf("%.3f", coef), stars),
    se_display = paste0("(", sprintf("%.3f", se), ")")
  )
}

# ---- Extract coefficient rows ----
results <- bind_rows(
  get_clustered_coef(ols_growth_lag,          "gdp_growth",     "(1) OLS"),
  get_clustered_coef(ols_growth_lag,          "gdp_growth_lag", "(1) OLS"),
  
  get_clustered_coef(ols_growth_lag_trend,    "gdp_growth",     "(2) OLS"),
  get_clustered_coef(ols_growth_lag_trend,    "gdp_growth_lag", "(2) OLS"),
  
  get_clustered_coef(ols_growth_lag_fe_trend, "gdp_growth",     "(3) OLS"),
  get_clustered_coef(ols_growth_lag_fe_trend, "gdp_growth_lag", "(3) OLS"),
  
  get_clustered_coef(iv_growth_lag,           "gdp_growth",     "(4) IV-2SLS"),
  get_clustered_coef(iv_growth_lag,           "gdp_growth_lag", "(4) IV-2SLS"),
  
  get_clustered_coef(iv_growth_lag_fe_trend,  "gdp_growth",     "(5) IV-2SLS"),
  get_clustered_coef(iv_growth_lag_fe_trend,  "gdp_growth_lag", "(5) IV-2SLS")
)

coef_rows <- results %>%
  mutate(
    variable = case_when(
      term == "gdp_growth" ~ "GDP Growth",
      term == "gdp_growth_lag" ~ "GDP Growth (t-1)"
    )
  ) %>%
  select(variable, model, coef_display) %>%
  pivot_wider(names_from = model, values_from = coef_display)

se_rows <- results %>%
  mutate(
    variable = case_when(
      term == "gdp_growth" ~ "",
      term == "gdp_growth_lag" ~ " "
    )
  ) %>%
  select(variable, model, se_display) %>%
  pivot_wider(names_from = model, values_from = se_display)

coef_table <- bind_rows(
  coef_rows[1, ],
  se_rows[1, ],
  coef_rows[2, ],
  se_rows[2, ]
)

# ---- Model info rows ----
info_rows <- tibble(
  variable = c(
    "Country Fixed Effects",
    "Linear Time Trend",
    "Observations",
    "R-squared"
  ),
  `(1) OLS` = c(
    "No",
    "No",
    nobs(ols_growth_lag),
    sprintf("%.3f", summary(ols_growth_lag)$r.squared)
    ),
  `(2) OLS` = c(
    "No",
    "Yes",
    nobs(ols_growth_lag_trend),
    sprintf("%.3f", summary(ols_growth_lag_trend)$r.squared)
      ),
  `(3) OLS` = c(
    "Yes",
    "Yes",
    nobs(ols_growth_lag_fe_trend),
    sprintf("%.3f", summary(ols_growth_lag_fe_trend)$r.squared)
      ),
  `(4) IV-2SLS` = c(
    "No",
    "No",
    nobs(iv_growth_lag),
    "---"
   ),
  `(5) IV-2SLS` = c(
    "Yes",
    "Yes",
    nobs(iv_growth_lag_fe_trend),
    "---"
  )
)

final_growth_conflict_table <- bind_rows(
  coef_table,
  info_rows
)

# ---- Style and save ----
growth_conflict_table_plot <- final_growth_conflict_table %>%
  gt() %>%
  tab_header(
    title = "Economic Growth and Civil Conflict",
    subtitle = "OLS and IV-2SLS estimates, civil conflict with 25+ deaths, 1981–2008"
  ) %>%
  cols_label(variable = "") %>%
  cols_align(align = "center", columns = -variable) %>%
  cols_align(align = "left", columns = variable) %>%
  cols_width(
    variable ~ px(120),
    everything() ~ px(75)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = variable %in% c("GDP Growth", "GDP Growth (t-1)"))
  ) %>%
  tab_style(
    style = cell_text(color = "gray40"),
    locations = cells_body(rows = variable %in% c("", " "))
  ) %>%
  tab_source_note(
    source_note = "Notes: Standard errors clustered at the country level are reported in parentheses to account for heteroskedasticity and within-country serial correlation. *** p < 0.01, ** p < 0.05, * p < 0.10. Dependent variable equals 1 for civil conflict with at least 25 battle deaths. IV-2SLS models instrument current and lagged GDP growth with current and lagged oil-price growth. Country-specific time trends are included in the appendix."
  ) %>%
  opt_row_striping() %>%
  tab_options(
    table.font.size = px(10),
    heading.title.font.size = "medium",
    data_row.padding = px(3)
  )

gtsave(
  growth_conflict_table_plot,
  filename = "economic_growth_conflict_table.png",
  path = "Tables"
)

# ============================================================
# ---- Create Figures for Paper ----
# Images Will Save in Tables Folder
# ============================================================

# ----Figure 1: Conflict Over Time Plot----
#Plot Conflict Over Time
conflict_plot <- panel_iv %>%
  filter(!is.na(conflict_25)) %>%
  ggplot(aes(x = year, y = country_label, fill = factor(conflict_25))) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c("0" = "gray90", "1" = "firebrick"),
    labels = c("No conflict", "Conflict"),
    name = ""
  ) +
  labs(
    title = "Civil Conflict Incidence by Country-Year",
    subtitle = "Binary indicator: 1 = civil conflict with at least 25 battle deaths",
    x = "Year",
    y = NULL
  ) +
  theme_minimal()

conflict_plot

#Save to Tables folder
ggsave(
  "conflict_incidence_plot.png",
  plot = conflict_plot,
  path = "Tables",
  width = 1920,
  height = 1080,
  units = "px"
)

# ----Figure 2: Oil Rents as % of GDP Over Time ----
#Plot Oil Rents by Country Over Time
#Plot
gdp_per %>%
  filter(!is.na(oil_rents)) %>%
  ggplot(aes(x = year, y = oil_rents, color = country)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Oil Rents as % of GDP Over Time",
    x = "Year",
    y = "Percent of GDP"
  ) +
  theme_minimal()
ggsave("oil_plot.png", path = "Tables", width = 10, height = 5.625, dpi = 300)
ggsave("oil_plot_paper.png", path = "Tables", width = 8, height = 6, dpi = 300)

# ============================================================
# ---- Appendix ----
# Other Tests Run for Robustness
# ============================================================

# ---- A1A. Country-specific Time Trends (OLS) ---- 
# Create centered year trend so coefficients are easier to read
panel_iv <- panel_iv %>%
  mutate(year_c = year - min(year, na.rm = TRUE))

# First stage with country-specific time trends
first_stage_cstrend <- lm(
  gdp_growth ~ oil_growth + oil_growth_lag + country + country:year_c,
  data = panel_iv
)

summary(first_stage_cstrend)

# Reduced form with country-specific time trends
reduced_form_cstrend <- lm(
  conflict_25 ~ oil_growth + oil_growth_lag + country + country:year_c,
  data = panel_iv
)

summary(reduced_form_cstrend)

# ---- A1B. Country-specific Time Trends (IV) ----
library(AER)

iv_cstrend <- ivreg(
  conflict_25 ~ gdp_growth + gdp_growth_lag + country + country:year_c |
    oil_growth + oil_growth_lag + country + country:year_c,
  data = panel_iv
)

summary(iv_cstrend, diagnostics = TRUE)

# ---- 13. Probit baseline ----
probit_model <- glm(
  conflict_25 ~ gdp_growth + country + year,
  data = panel_iv,
  family = binomial(link = "probit")
)

summary(probit_model)

# ============================================================
# ---- Appendix Tables ----
# Country-specific trends and serial correlation diagnostics
# ============================================================

library(dplyr)
library(tidyr)
library(gt)
library(lmtest)
library(sandwich)
library(AER)
library(car)
library(plm)

dir.create("Tables", showWarnings = FALSE)

# ------------------------------------------------------------
# Create appendix-only dataset so main paper models are untouched
# ------------------------------------------------------------

panel_iv_appendix <- panel_iv %>%
  mutate(year_c = year - min(year, na.rm = TRUE))

# ------------------------------------------------------------
# Helper function: clustered coefficient + SE
# ------------------------------------------------------------

get_app_clustered_coef <- function(model, term_name, model_name) {
  
  clustered <- coeftest(
    model,
    vcov = vcovCL(model, cluster = model.frame(model)$country)
  )
  
  if (!(term_name %in% rownames(clustered))) {
    return(tibble(
      model = model_name,
      term = term_name,
      coef_display = "",
      se_display = ""
    ))
  }
  
  coef <- clustered[term_name, "Estimate"]
  se   <- clustered[term_name, "Std. Error"]
  p    <- clustered[term_name, "Pr(>|t|)"]
  
  stars <- case_when(
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
  
  tibble(
    model = model_name,
    term = term_name,
    coef_display = paste0(sprintf("%.3f", coef), stars),
    se_display = paste0("(", sprintf("%.3f", se), ")")
  )
}

# ------------------------------------------------------------
# Helper function: clustered F-test
# ------------------------------------------------------------

get_app_clustered_f <- function(model, terms_to_test) {
  
  available_terms <- terms_to_test[terms_to_test %in% names(coef(model))]
  
  if (length(available_terms) == 0) {
    return("---")
  }
  
  restrictions <- paste0(available_terms, " = 0")
  
  test <- car::linearHypothesis(
    model,
    restrictions,
    vcov. = vcovCL(model, cluster = model.frame(model)$country),
    test = "F"
  )
  
  sprintf("%.2f", test$F[2])
}

# ============================================================
# APPENDIX TABLE A1: FIRST STAGE WITH COUNTRY-SPECIFIC TRENDS
# ============================================================

app_a1_first_stage <- lm(
  gdp_growth ~ oil_growth + oil_growth_lag + country + country:year_c,
  data = panel_iv_appendix
)

app_a1_results <- bind_rows(
  get_app_clustered_coef(app_a1_first_stage, "oil_growth",     "(1)"),
  get_app_clustered_coef(app_a1_first_stage, "oil_growth_lag", "(1)")
)

app_a1_coef_rows <- app_a1_results %>%
  mutate(
    variable = case_when(
      term == "oil_growth" ~ "Oil Growth",
      term == "oil_growth_lag" ~ "Oil Growth (t-1)"
    )
  ) %>%
  select(variable, model, coef_display) %>%
  pivot_wider(names_from = model, values_from = coef_display)

app_a1_se_rows <- app_a1_results %>%
  mutate(
    variable = case_when(
      term == "oil_growth" ~ "",
      term == "oil_growth_lag" ~ " "
    )
  ) %>%
  select(variable, model, se_display) %>%
  pivot_wider(names_from = model, values_from = se_display)

app_a1_table <- bind_rows(
  app_a1_coef_rows[1, ],
  app_a1_se_rows[1, ],
  app_a1_coef_rows[2, ],
  app_a1_se_rows[2, ],
  tibble(
    variable = c(
      "Country Fixed Effects",
      "Country-Specific Linear Trends",
      "Observations",
      "R-squared",
      "Clustered F-stat"
    ),
    `(1)` = c(
      "Yes",
      "Yes",
      nobs(app_a1_first_stage),
      sprintf("%.3f", summary(app_a1_first_stage)$r.squared),
      get_app_clustered_f(app_a1_first_stage, c("oil_growth", "oil_growth_lag"))
    )
  )
)

app_a1_plot <- app_a1_table %>%
  gt() %>%
  tab_header(
    title = "Appendix Table A1",
    subtitle = "First Stage with Country-Specific Linear Time Trends"
  ) %>%
  cols_label(variable = "") %>%
  cols_align(align = "center", columns = -variable) %>%
  cols_align(align = "left", columns = variable) %>%
  cols_width(
    variable ~ px(220),
    everything() ~ px(100)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = variable %in% c("Oil Growth", "Oil Growth (t-1)"))
  ) %>%
  tab_style(
    style = cell_text(color = "gray40"),
    locations = cells_body(rows = variable %in% c("", " "))
  ) %>%
  tab_source_note(
    source_note = "Notes: Dependent variable is annual GDP growth. Standard errors clustered at the country level are reported in parentheses. Country-specific linear trends are country-by-year trend interactions. *** p < 0.01, ** p < 0.05, * p < 0.10."
  ) %>%
  opt_row_striping() %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = "medium",
    data_row.padding = px(3)
  )

gtsave(
  app_a1_plot,
  filename = "appendix_A1_first_stage_cstrend.png",
  path = "Tables"
)

# ============================================================
# APPENDIX TABLE A2: REDUCED FORM WITH COUNTRY-SPECIFIC TRENDS
# ============================================================

app_a2_reduced_form <- lm(
  conflict_25 ~ oil_growth + oil_growth_lag + country + country:year_c,
  data = panel_iv_appendix
)

app_a2_results <- bind_rows(
  get_app_clustered_coef(app_a2_reduced_form, "oil_growth",     "(1)"),
  get_app_clustered_coef(app_a2_reduced_form, "oil_growth_lag", "(1)")
)

app_a2_coef_rows <- app_a2_results %>%
  mutate(
    variable = case_when(
      term == "oil_growth" ~ "Oil Growth",
      term == "oil_growth_lag" ~ "Oil Growth (t-1)"
    )
  ) %>%
  select(variable, model, coef_display) %>%
  pivot_wider(names_from = model, values_from = coef_display)

app_a2_se_rows <- app_a2_results %>%
  mutate(
    variable = case_when(
      term == "oil_growth" ~ "",
      term == "oil_growth_lag" ~ " "
    )
  ) %>%
  select(variable, model, se_display) %>%
  pivot_wider(names_from = model, values_from = se_display)

app_a2_table <- bind_rows(
  app_a2_coef_rows[1, ],
  app_a2_se_rows[1, ],
  app_a2_coef_rows[2, ],
  app_a2_se_rows[2, ],
  tibble(
    variable = c(
      "Country Fixed Effects",
      "Country-Specific Linear Trends",
      "Observations",
      "R-squared"
    ),
    `(1)` = c(
      "Yes",
      "Yes",
      nobs(app_a2_reduced_form),
      sprintf("%.3f", summary(app_a2_reduced_form)$r.squared)
    )
  )
)

app_a2_plot <- app_a2_table %>%
  gt() %>%
  tab_header(
    title = "Appendix Table A2",
    subtitle = "Reduced Form with Country-Specific Linear Time Trends"
  ) %>%
  cols_label(variable = "") %>%
  cols_align(align = "center", columns = -variable) %>%
  cols_align(align = "left", columns = variable) %>%
  cols_width(
    variable ~ px(220),
    everything() ~ px(100)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = variable %in% c("Oil Growth", "Oil Growth (t-1)"))
  ) %>%
  tab_style(
    style = cell_text(color = "gray40"),
    locations = cells_body(rows = variable %in% c("", " "))
  ) %>%
  tab_source_note(
    source_note = "Notes: Dependent variable equals 1 for civil conflict with at least 25 battle deaths. Standard errors clustered at the country level are reported in parentheses. Country-specific linear trends are country-by-year trend interactions. *** p < 0.01, ** p < 0.05, * p < 0.10."
  ) %>%
  opt_row_striping() %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = "medium",
    data_row.padding = px(3)
  )

gtsave(
  app_a2_plot,
  filename = "appendix_A2_reduced_form_cstrend.png",
  path = "Tables"
)

# ============================================================
# APPENDIX TABLE A3: OLS AND IV WITH COUNTRY-SPECIFIC TRENDS
# ============================================================

app_a3_ols_cstrend <- lm(
  conflict_25 ~ gdp_growth + gdp_growth_lag + country + country:year_c,
  data = panel_iv_appendix
)

app_a3_iv_cstrend <- ivreg(
  conflict_25 ~ gdp_growth + gdp_growth_lag + country + country:year_c |
    oil_growth + oil_growth_lag + country + country:year_c,
  data = panel_iv_appendix
)

app_a3_results <- bind_rows(
  get_app_clustered_coef(app_a3_ols_cstrend, "gdp_growth",     "(1) OLS"),
  get_app_clustered_coef(app_a3_ols_cstrend, "gdp_growth_lag", "(1) OLS"),
  
  get_app_clustered_coef(app_a3_iv_cstrend,  "gdp_growth",     "(2) IV-2SLS"),
  get_app_clustered_coef(app_a3_iv_cstrend,  "gdp_growth_lag", "(2) IV-2SLS")
)

app_a3_coef_rows <- app_a3_results %>%
  mutate(
    variable = case_when(
      term == "gdp_growth" ~ "GDP Growth",
      term == "gdp_growth_lag" ~ "GDP Growth (t-1)"
    )
  ) %>%
  select(variable, model, coef_display) %>%
  pivot_wider(names_from = model, values_from = coef_display)

app_a3_se_rows <- app_a3_results %>%
  mutate(
    variable = case_when(
      term == "gdp_growth" ~ "",
      term == "gdp_growth_lag" ~ " "
    )
  ) %>%
  select(variable, model, se_display) %>%
  pivot_wider(names_from = model, values_from = se_display)

app_a3_table <- bind_rows(
  app_a3_coef_rows[1, ],
  app_a3_se_rows[1, ],
  app_a3_coef_rows[2, ],
  app_a3_se_rows[2, ],
  tibble(
    variable = c(
      "Country Fixed Effects",
      "Country-Specific Linear Trends",
      "Observations",
      "R-squared"
    ),
    `(1) OLS` = c(
      "Yes",
      "Yes",
      nobs(app_a3_ols_cstrend),
      sprintf("%.3f", summary(app_a3_ols_cstrend)$r.squared)
    ),
    `(2) IV-2SLS` = c(
      "Yes",
      "Yes",
      nobs(app_a3_iv_cstrend),
      "---"
    )
  )
)

app_a3_plot <- app_a3_table %>%
  gt() %>%
  tab_header(
    title = "Appendix Table A3",
    subtitle = "Economic Growth and Civil Conflict with Country-Specific Linear Time Trends"
  ) %>%
  cols_label(variable = "") %>%
  cols_align(align = "center", columns = -variable) %>%
  cols_align(align = "left", columns = variable) %>%
  cols_width(
    variable ~ px(220),
    everything() ~ px(110)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = variable %in% c("GDP Growth", "GDP Growth (t-1)"))
  ) %>%
  tab_style(
    style = cell_text(color = "gray40"),
    locations = cells_body(rows = variable %in% c("", " "))
  ) %>%
  tab_source_note(
    source_note = "Notes: Dependent variable equals 1 for civil conflict with at least 25 battle deaths. Standard errors clustered at the country level are reported in parentheses. IV-2SLS instruments current and lagged GDP growth with current and lagged oil-price growth. Country-specific linear trends are country-by-year trend interactions. *** p < 0.01, ** p < 0.05, * p < 0.10."
  ) %>%
  opt_row_striping() %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = "medium",
    data_row.padding = px(3)
  )

gtsave(
  app_a3_plot,
  filename = "appendix_A3_ols_iv_cstrend.png",
  path = "Tables"
)

# ============================================================
# APPENDIX TABLE B1: WOOLDRIDGE TEST FOR SERIAL CORRELATION
# ============================================================

pdata_appendix <- pdata.frame(
  panel_iv_appendix,
  index = c("country", "year")
)

wooldridge_test <- pwartest(
  conflict_25 ~ gdp_growth,
  data = pdata_appendix
)

format_p_value <- function(p) {
  ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
}

app_b1_table <- tibble(
  statistic = c(
    "Test",
    "Null hypothesis",
    "Alternative hypothesis",
    "F-statistic",
    "p-value",
    "Conclusion"
  ),
  value = c(
    "Wooldridge test for serial correlation in fixed-effects panels",
    "No first-order serial correlation",
    "Serial correlation present",
    sprintf("%.2f", unname(wooldridge_test$statistic)),
    format_p_value(wooldridge_test$p.value),
    "Reject null; cluster standard errors at the country level"
  )
)

app_b1_plot <- app_b1_table %>%
  gt() %>%
  tab_header(
    title = "Appendix Table B1",
    subtitle = "Serial Correlation Diagnostic"
  ) %>%
  cols_label(
    statistic = "",
    value = ""
  ) %>%
  cols_align(align = "left", columns = everything()) %>%
  cols_width(
    statistic ~ px(180),
    value ~ px(360)
  ) %>%
  tab_source_note(
    source_note = "Notes: The Wooldridge test is applied to the country-year panel. Rejection of the null indicates within-country serial correlation, motivating country-clustered standard errors in the main tables."
  ) %>%
  opt_row_striping() %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = "medium",
    data_row.padding = px(4)
  )

gtsave(
  app_b1_plot,
  filename = "appendix_B_wooldridge_test.png",
  path = "Tables"
)
