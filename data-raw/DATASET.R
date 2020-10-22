## code to prepare `DATASET` dataset goes here
library(tidyverse)

dat <- readRDS(here::here("data-raw", "combined_model_inputs.RDS"))

# use PST troll data (complete seasonal cycle)
comp_dat <- dat$comp_long[[1]]
catch_dat <- dat$catch_data[[1]]

#subset by genetics sampling to increase speed
samp_events <- comp_dat %>% 
  pull(sample_id) %>% 
  unique() %>% 
  sample(., size = 150, replace = FALSE)

comp_ex <- comp_dat %>% 
  filter(sample_id %in% samp_events) 

# strata to retain from catch data 
strata_key <- comp_ex %>% 
  mutate(key = paste(region_c, month_n, as.character(year), sep = "_")) %>% 
  pull(key)

# subset catch data
catch_ex <- catch_dat %>% 
  mutate(key = paste(region_c, month_n, as.character(year), sep = "_")) %>% 
  filter(key %in% strata_key) %>% 
  select(-eff_z, -key)

# save
usethis::use_data(comp_ex, overwrite = TRUE)
usethis::use_data(catch_ex, overwrite = TRUE)
