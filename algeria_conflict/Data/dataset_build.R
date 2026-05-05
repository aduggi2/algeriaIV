# =========================
# REPLICATION DATASET BUILDER
# Algeria, 1980-2024
# =========================

# ---- 1. Install packages (run once) ----
install.packages(c("WDI", "fredr", "dplyr", "tidyr", "lubridate", "readr", "imfweo", "AER"))

# ---- 2. Load libraries ----
library(WDI)
library(AER)
library(fredr)
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(imfweo)
library(readxl)

# ---- 3. Set your FRED API key ----
# Get one free at: https://fred.stlouisfed.org/
# Then either paste it here or store it in your .Renviron
fredr_set_key("16e9a33f4ff54715fc9eef51e702d492")

# ---- 4. Pull World Bank data ----
# Algeria = DZA
# GDP growth = NY.GDP.MKTP.KD.ZG
# Population growth = SP.POP.GROW

wb <- WDI(
  country = "DZA",
  indicator = c(
    gdp_growth = "NY.GDP.MKTP.KD.ZG",
    pop_growth = "SP.POP.GROW"
  ),
  start = 1980,
  end = 2024,
  extra = FALSE
)

wb <- wb %>%
  transmute(
    year = as.integer(year),
    gdp_growth = gdp_growth,
    pop_growth = pop_growth
  ) %>%
  arrange(year)

print("World Bank data:")
print(head(wb, 10))

# ---- 5. Pull Annual APSP crude oil prices from IMF manually ----
# API from IMF appears to be broken at this time#
#Source:https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.RES:WEO(9.0.0)&INDICATOR=POILAPSP

oil_imf <- read_csv("Data/oil_imf.csv")

# ---- 6. Clean IMF Data ----
oil_imf_clean <- oil_imf %>%
  transmute(
    year = as.integer(TIME_PERIOD),
    oil_price = OBS_VALUE
  ) %>%
  arrange(year)

# ---- 7. Create Oil Shock Variable  ----
oil_imf_clean <- oil_imf_clean %>%
  mutate(
    oil_growth = c(NA, diff(log(oil_price)))
  )

# ---- 8. Merge Oil & World Bank Datasets  ----
data <- wb %>%
  inner_join(oil_imf_clean, by = "year") %>%
  arrange(year)

data_clean <- data %>%
  filter(!is.na(oil_growth))

# ----First-Stage Regression to Test Relevance----
first_stage <- lm(gdp_growth ~ oil_growth, data = data_clean)

summary(first_stage)

#----Run a second first-stage to test adding pop growth
first_stage2 <- lm(gdp_growth ~ oil_growth + pop_growth, data = data_clean)

summary(first_stage2)

#----Lag the oil shock----
data_clean <- data_clean %>%
  mutate(oil_growth_lag = lag(oil_growth))

summary(lm(gdp_growth ~ oil_growth_lag, data = data_clean))

# ---- 7. Create conflict variable manually ----
#Download the UCDP/RIO Armed Conflict Dataset from:
#https://cdn.cloud.prio.org/files/0cecffa6-1603-4418-a361-e24ee0aeed40/UCDP%20PRIO%20Armed%20Conflict%20Dataset%20v4-2009.zip?inline=true

conflict_data <- read_excel("Data/Main Conflict Table.xls")

#Filter Algeria
conflict_dza <- conflict_data %>%
  filter(
    grepl("Algeria", Location),   # country
    Type %in% c(3, 4)             # intrastate conflicts only
  )

# Keep only year
conflict_dza_year <- conflict_dza %>%
  select(YEAR) %>%
  distinct()

#Create Dummy = 1
conflict_dza_year <- conflict_dza_year %>%
  mutate(conflict_25 = 1)

#Create full year panel
conflict_full <- data.frame(year = 1980:2024) %>%
  left_join(conflict_dza_year, by = c("year" = "YEAR")) %>%
  mutate(
    conflict_25 = ifelse(is.na(conflict_25), 0, 1)
  )

# ----8. Merge Datasets ----
data_iv <- data_clean %>%
  left_join(conflict_full, by = "year")

# ----9. Baseline OLS ----
ols <- lm(conflict_25 ~ gdp_growth, data = data_iv)
summary(ols)

# ----10. Run IV/2SLS ----
library(AER)

iv_model <- ivreg(conflict_25 ~ gdp_growth | oil_growth, data = data_iv)
summary(iv_model)

# --- 1. Countries to include ---
countries <- c("DZA", "LBY", "EGY", "TUN", "MAR")

# --- 2. Pull World Bank data for all 5 countries at once ---
wb_panel <- WDI(
  country = countries,
  indicator = c(
    gdp_growth = "NY.GDP.MKTP.KD.ZG",
    pop_growth = "SP.POP.GROW"
  ),
  start = 1980,
  end = 2024,
  extra = FALSE
) %>%
  transmute(
    country = iso3c,
    year = as.integer(year),
    gdp_growth = gdp_growth,
    pop_growth = pop_growth
  ) %>%
  arrange(country, year)

# --- 3. Expand your IMF oil data across all countries ---
oil_panel <- expand.grid(
  country = countries,
  year = oil_imf_clean$year
) %>%
  left_join(oil_imf_clean, by = "year") %>%
  arrange(country, year)

# --- 4. Build PRIO conflict dummy for all 5 countries ---
conflict_panel <- conflict_data %>%
  filter(
    Location %in% c("Algeria", "Libya", "Egypt", "Tunisia", "Morocco"),
    Type %in% c(3, 4)
  ) %>%
  transmute(
    country = case_when(
      Location == "Algeria" ~ "DZA",
      Location == "Libya"   ~ "LBY",
      Location == "Egypt"   ~ "EGY",
      Location == "Tunisia" ~ "TUN",
      Location == "Morocco" ~ "MAR"
    ),
    year = as.integer(Year),
    conflict_25 = 1
  ) %>%
  distinct(country, year, .keep_all = TRUE)

# --- 5. Create full country-year panel for conflict zeros ---
conflict_full_panel <- expand.grid(
  country = countries,
  year = 1980:2024
) %>%
  left_join(conflict_panel, by = c("country", "year")) %>%
  mutate(conflict_25 = ifelse(is.na(conflict_25), 0, conflict_25)) %>%
  arrange(country, year)

# --- 6. Merge everything ---
panel_iv <- wb_panel %>%
  left_join(oil_panel, by = c("country", "year")) %>%
  left_join(conflict_full_panel, by = c("country", "year")) %>%
  filter(!is.na(oil_growth)) %>%
  arrange(country, year)

# --- 7. Check it ---
head(panel_iv, 20)
table(panel_iv$country)
table(panel_iv$conflict_25)
summary(panel_iv)