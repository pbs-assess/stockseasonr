### Sandbox for increasing function flexibility
## Example with dirichlet and random intercepts
## March 28, 2022
## NOTE: change working directories as needed

library(tidyverse)
library(TMB)
library(stockseasonr)

# tmb models (functionally this is equivalent to PeerJ pub but does look 
# different inside)
compile(here::here("src", "dev_tmb", "dirichlet_ri.cpp"))
dyn.load(dynlib(here::here("src", "dev_tmb", "dirichlet_ri")))


# utility functions for prepping smooths 
source(here::here("R", "smooth_utils.R"))
# data prep and model fitting functions
source(here::here("R", "versatile_fit.R"))


comp1 <- readRDS(here::here("data", "rec_ex_data.rds"))

stock_comp <- comp1 %>%
  group_by(sample_id, reg, reg_c = cap_region, month_n, year, sample_size,
           pst_agg) %>%
  summarize(prob = sum(prob), .groups = "drop") %>%
  ungroup() %>%
  droplevels() %>%
  filter(reg %in% c("JdFS", "SSoG")) %>%
  mutate(year = as.factor(year),
         reg = as.factor(reg)
  ) %>%
  rename(agg = pst_agg)

# generate a predictive dataframe to estimate changes in stock composition over 
# the seasonal cycle; months may vary by region so split then recombine to make
# sure they are representative
pred_dat_stock_comp <- group_split(stock_comp, reg) %>%
  map_dfr(., function (x) {
    expand.grid(
      reg = unique(x$reg),
      month_n = seq(min(x$month_n),
                    max(x$month_n),
                    by = 0.1
      )
    )
  }) %>%
  mutate(
    reg_month = paste(reg, month_n, sep = "_"),
    key_var = as.factor(reg_month)
  ) %>%
  arrange(reg, month_n) %>%
  droplevels()


## FIT MODEL -------------------------------------------------------------------

# generate model inputs for TMB
model_inputs <- make_inputs(
  # fixed effects
  comp_formula = agg ~ reg + s(month_n, bs = "cc", k = 4, by = reg, m = 2),
  # input data
  comp_dat = stock_comp,
  # random intercepts
  comp_rint = "year",
  # predictive data
  pred_comp = pred_dat_stock_comp,
  # model type (eventually will includedirichlet, negbin, integrated)
  model = "dirichlet",
  # exclude random intercepts from predictions (see note for details)
  include_re_preds = FALSE
)

# uses fit_random = FALSE to ignore random intercepts; fits quickly
mod_fe <- fit_model(
  tmb_data = model_inputs$tmb_data, 
  tmb_pars = model_inputs$tmb_pars, 
  tmb_map = model_inputs$tmb_map, 
  tmb_random  = model_inputs$tmb_random,
  model = "dirichlet",
  fit_random = FALSE
)

# random effects model
mod_re <- fit_model(
  tmb_data = model_inputs$tmb_data, 
  tmb_pars = model_inputs$tmb_pars, 
  tmb_map = model_inputs$tmb_map, 
  tmb_random  = model_inputs$tmb_random,
  model = "dirichlet",
  fit_random = TRUE,
  # sometimes using the fixed effect estimates as initial values can make 
  # model fitting more difficult. If that happens try using ignore_fix = TRUE
  ignore_fix = FALSE,
  include_re_preds = FALSE
)


## EVALUATE MODEL PREDS --------------------------------------------------------

# summarize predictions and parameter estimates
ssdr <- mod_re$ssdr

# look at what's available
unique(rownames(ssdr))


# pull out predicted stock composition (in link space)
logit_pred_ppn <- ssdr[rownames(ssdr) == "logit_pred_Pi_prop", ]

# convert to real space
link_preds <- data.frame(
  link_prob_est = logit_pred_ppn[ , "Estimate"],
  link_prob_se =  logit_pred_ppn[ , "Std. Error"]
) %>% 
  mutate(
    pred_prob_est = plogis(link_prob_est),
    pred_prob_low = plogis(link_prob_est + (qnorm(0.025) * link_prob_se)),
    pred_prob_up = plogis(link_prob_est + (qnorm(0.975) * link_prob_se))
  ) 

# bind output matrices with predictive datasets (iterate over to represent 
# different stocks)
stock_seq <- colnames(model_inputs$tmb_data$Y2_ik)
pred_comp <- purrr::map(stock_seq, function (x) {
  dum <- pred_dat_stock_comp
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) 

# plot predictions
p <- ggplot(data = pred_comp, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_wrap(~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est, colour = reg))
  
p_ribbon <- p +
  geom_ribbon(data = pred_comp,
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = reg),
              alpha = 0.2)



