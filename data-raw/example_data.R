## code to prepare `DATASET` dataset goes here
library(tidyverse)

dat <- readRDS(here::here("data-raw", "combined_model_inputs.RDS"))

# use PST troll data (complete seasonal cycle)
comp_full <- dat$comp_long[[1]]
catch_full <- dat$catch_data[[1]]

month_seq <- seq(from = 1, to = 12, by = 2)
year_seq <- seq(from = 2007, to = 2012, by = 1)

#subset by genetics sampling to increase speed
# samp_events <- comp_full %>% 
#   pull(sample_id) %>% 
#   unique() %>% 
#   sample(., size = 150, replace = FALSE)

comp_ex <- comp_full %>% 
  filter(
    #sample_id %in% samp_events
    month %in% month_seq,
    year %in% year_seq
    ) %>% 
  droplevels()

# strata to retain from catch data 
# strata_key <- comp_ex %>% 
#   mutate(key = paste(region_c, month_n, as.character(year), sep = "_")) %>% 
#   pull(key)

# subset catch data
catch_ex <- catch_full %>% 
  # mutate(key = paste(region_c, month_n, as.character(year), sep = "_")) %>% 
  filter(
    #key %in% strata_key
    month %in% month_seq,
    year %in% year_seq
  ) %>% 
  droplevels() %>% 
  select(-eff_z#, -key
         )

# save
usethis::use_data(comp_ex, overwrite = TRUE)
usethis::use_data(catch_ex, overwrite = TRUE)
