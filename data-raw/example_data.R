## code to prepare `DATASET` dataset goes here
library(tidyverse)

dat <- readRDS(here::here("data-raw", "combined_model_inputs.RDS"))

# use PST troll data (complete seasonal cycle)
comp_full <- dat$comp_long[[1]]
catch_full <- dat$catch_data[[1]]

# drop some months and years to decrease size
month_seq <- seq(from = 1, to = 12, by = 2)
year_seq <- seq(from = 2007, to = 2011, by = 1)

comp_ex <- comp_full %>% 
  filter(
    month %in% month_seq,
    year %in% year_seq
    ) %>% 
  droplevels() %>% 
  select(-gear, -region_c, -month)

# subset catch data
catch_ex <- catch_full %>% 
  filter(
    month %in% month_seq,
    year %in% year_seq
  ) %>% 
  droplevels() %>% 
  select(-area_n, -region_c, -month, -eff, -eff_z)

# save
usethis::use_data(comp_ex, overwrite = TRUE)
usethis::use_data(catch_ex, overwrite = TRUE)
