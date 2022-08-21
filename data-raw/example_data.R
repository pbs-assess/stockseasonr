## code to prepare `DATASET` dataset goes here
library(dplyr)

dat <- readRDS(here::here("data-raw", "combined_model_inputs.RDS"))

# use PST troll data (complete seasonal cycle)
comp_full <- dat$comp_long[[1]]
catch_full <- dat$catch_data[[1]]

# drop some months and years to decrease size
month_seq <- seq(from = 1, to = 12, by = 2)
year_seq <- seq(from = 2007, to = 2014, by = 1)

comp_ex <- comp_full %>% 
  filter(
    month %in% month_seq,
    year %in% year_seq
    ) %>% 
  #aggregate some stocks for plotting purposes
  mutate(
    agg = case_when(
      grepl("_sp", agg) ~ "CR_spring",
      grepl("CR", agg) ~ "CR_su/fa",
      agg %in% c("CA_ORCST", "WACST", "NBC_SEAK", "WCVI") ~ "other",
      TRUE ~ agg
    )
  ) %>% 
  group_by(sample_id, region, year, month_n, agg, nn) %>% 
  summarize(prob = sum(agg_prob),
            .groups = "drop") %>% 
  ungroup() %>% 
  droplevels() 

# subset catch data
catch_ex <- catch_full %>% 
  # filter(
  #   month %in% month_seq,
  #   year %in% year_seq
  # ) %>% 
  # consolidate by region
  group_by(region, month_n, year) %>% 
  summarize(eff = sum(eff),
            catch = sum(catch),
            .groups = "drop") %>% 
  ungroup() %>% 
  droplevels() 

# save
usethis::use_data(comp_ex, overwrite = TRUE)
usethis::use_data(catch_ex, overwrite = TRUE)
