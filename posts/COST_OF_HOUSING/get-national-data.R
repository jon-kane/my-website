source(here::here("_common.R"))
library(tidyquant)

# data_series = c(
#   "CPIAUCSL",     # Consumer Price Index for All Urban Consumers: All Items in U.S. City Average
#   "MSPUS",        # Median Sales Price of Houses Sold for the United States
#   "MORTGAGE30US", # 30-Year Fixed Rate Mortgage Average in the United States
#   "MEHOINUSA646N" # Median Household Income in the United States
# )

# ------------------------------------------------------------------------------
# Consumer Price Index for All Urban Consumers: All Items in U.S. City Average
# This data from FRED goes back to 1947
tb_n1 = tq_get(
  "CPIAUCSL", 
  get = "economic.data", 
  from = as.Date("1900-01-01"), 
  to = Sys.Date()
  )

# ------------------------------------------------------------------------------
# Median Household Income in the United States
# This data goes back to 1984
# I'll augment this data from some census sources
tb_n2 = tq_get(
  "MEHOINUSA646N", 
  get = "economic.data", 
  from = as.Date("1900-01-01"), 
  to = Sys.Date()
  )

# 1947:1965: https://www2.census.gov/library/publications/1999/compendia/statab/119ed/tables/sec31.pdf p.877
# 1967-1983: https://www2.census.gov/programs-surveys/demo/tables/p60/276/tableD1.xlsx
# data pulled on 2026/01/07

# data in nominal amounts
tb_n2.1 = tibble(
  symbol = rep("MEHOINUSA646N", times = 5),
  date = c(
    as_date("1947-01-01"), 
    as_date("1950-01-01"), 
    as_date("1955-01-01"), 
    as_date("1960-01-01"), 
    as_date("1965-01-01")
  ),
  price = c(3031, 3319, 4418, 5620, 6957)
)

tb_n2.2 = tibble(
  symbol = rep("MEHOINUSA646N", times = 17),
  date = as_date("1967-01-01") + years(0:16),
  price = c(
    7143, 7743, 8389, 8734, 9028, 
    9697, 10512, 11197, 11800, 12686, 
    13572, 15064, 16461, 17710, 19074, 
    20171, 20885
  )
)

# ------------------------------------------------------------------------------
# 30-Year Fixed Rate Mortgage Average in the United States
tb_n3 = tq_get(
  "MORTGAGE30US", 
  get = "economic.data", 
  from = as.Date("1900-01-01"), 
  to = Sys.Date()
)

# ------------------------------------------------------------------------------
# Median Sales Price of Houses Sold for the United States
# This data from FRED goes back to 1963
# I'll augment this data using data accessed via this website: https://dqydj.com/historical-home-prices/
tb_n4 = tq_get(
  "MSPUS", 
  get = "economic.data", 
  from = as.Date("1900-01-01"), 
  to = Sys.Date()
)

tb_n4.1 = read_csv(
  here::here("posts/COST_OF_HOUSING/mspus-dqydj.csv"),
  show_col_types = F
  ) %>% 
  select(1:2) %>% 
  rename(date = category) %>% 
  rename(price = `Median Home Price (NSA)`) %>% 
  mutate(date = as_date(date, format = "%a %b %d %Y")) %>% 
  mutate(symbol = "MSPUS") %>% 
  relocate(symbol) %>% 
  filter(year(date) < 1963)

# ------------------------------------------------------------------------------
tb_n = bind_rows(
  tb_n1,
  tb_n2,
  tb_n2.1,
  tb_n2.2,
  tb_n3,
  tb_n4,
  tb_n4.1
) %>% 
  arrange(symbol, date) 

# get the max-min date to consider
min_date = tb_n %>% 
  filter(symbol != "MORTGAGE30US") %>% 
  group_by(symbol) %>% 
  summarize(min_date = min(date)) %>% 
  pull(min_date) %>% 
  max()

max_date = tb_n %>% 
  pull(date) %>% 
  max()
# ------------------------------------------------------------------------------
# Model the data

library(timetk)     # for easy time series padding/imputation
library(fable)      # for tidy forecasting
library(tsibble)    # fhe time-series tibble structure
library(distributional) # to extract SD from forecast distributions

# 2. Define a function to harmonize frequencies
# This function detects the data type and standardizes it to Monthly
harmonizeFrequency = function(data) {
  
  # Group by symbol to handle each series individually
  data_grouped = data %>%
    group_by(symbol) 
  
  # Step A: Summarize to Monthly (Handles Weekly -> Monthly)
  # This downsamples high freq and sets the grid for low freq
  data_monthly = data_grouped %>%
    summarise_by_time(
      .date_var = date,
      .by       = "month",
      price     = mean(price, na.rm = TRUE)
    )
  
  # Step B: Pad and Impute (Handles Quarterly/Annual -> Monthly)
  # This fills in the missing months for Q and A series
  data_harmonized = data_monthly %>%
    pad_by_time(date, .by = "month") %>%
    mutate(price = zoo::na.spline(price)) %>%
    ungroup()
  
  return(data_harmonized)
}

# 3. Process the data
# Convert to tsibble (time-series tibble) required for fable
harmonized_data = harmonizeFrequency(tb_n) %>%
  mutate(date = yearmonth(date)) %>%
  as_tsibble(index = date, key = symbol) %>% 
  group_by(symbol) %>%
  mutate(
    price_smooth = stats::loess(price ~ as.numeric(date), span = 0.05)$fitted
  ) %>%
  ungroup()

# 4. Build Models
# We use an automatic ARIMA model for each series
models = harmonized_data %>%
  model(
    arima = ARIMA(price_smooth)
  )

# 5. Forecast
# Forecast 2 years (24 months) ahead
forecasts = models %>%
  forecast(h = "2 years")

# 6. Extract Mean and Prediction SD
# fable returns a distribution; we extract the SD from it
final_predictions = forecasts %>%
  hilo(level = 95) %>% # Optional: Calculate intervals if needed
  mutate(
    pred_mean = .mean,
    pred_sd   = sqrt(variance(price_smooth)) # Extract sigma from the distribution
  ) %>%
  select(symbol, date, pred_mean, pred_sd)

# join the final predictions with the harmonized data

modeled_data = harmonized_data %>% 
  as_tibble() %>% 
  mutate(pred_sd = 0) %>% 
  bind_rows(
    final_predictions %>% 
      rename(price_smooth = pred_mean)
  ) %>% 
  filter(between(date, min_date, max_date)) %>% 
  mutate(date = as_date(date)) %>% 
  rename(price_raw = price) %>% 
  rename(price = price_smooth)

data_list = list(
  tb_n  = tb_n, # national data, raw
  tb_nm = modeled_data # national data, modeled
)

saveRDS(data_list, here::here("posts/COST_OF_HOUSING/national-data.rds"))
