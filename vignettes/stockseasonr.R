## ----setup, include = FALSE, cache = FALSE------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.asp = 0.618
)

## ----packages, echo = FALSE, message = FALSE----------------------------------
library(stockseasonr)

## ----echo=TRUE----------------------------------------------------------------
dat_in <- fit_stockseasonr(abund_formula = catch ~ 1 +
                             s(month_n, bs = "tp", k = 3, m = 2) +
                             region +
                             (1 | year),
                           abund_dat = catch_ex,
                           comp_formula = agg ~ 1 + region +
                             s(month_n, bs = "tp", k = 4, m = 2) +
                             (1 | year),
                           comp_dat = comp_ex,
                           model = "integrated",
                           fit = FALSE)

## ----echo = TRUE--------------------------------------------------------------
new_dat <- expand.grid(
  month_n = seq(1, 12, by = 0.1),
  region = unique(comp_ex$region)
)

m1 <- fit_stockseasonr(abund_formula = catch ~ 1 +
                             s(month_n, bs = "tp", k = 3, m = 2) +
                             region +
                             (1 | year),
                           abund_dat = catch_ex,
                           comp_formula = agg ~ 1 + region +
                             s(month_n, bs = "tp", k = 4, m = 2) +
                             (1 | year),
                           comp_dat = comp_ex,
                           model = "integrated",
                       pred_dat = new_dat,
                           random_walk = TRUE,
                           fit = TRUE)

## ----echo = TRUE--------------------------------------------------------------
theta <- as.list(m1$sdr, "Estimate")

# abundance estimates for regional effects
theta$b1_j
#composition fixed effect estimates
theta$B2_jk

## ----echo = TRUE--------------------------------------------------------------
ssdr <- m1$ssdr
unique(rownames(ssdr))

# predicted total abundance
head(ssdr[rownames(ssdr) %in% "log_pred_mu1", ])
# predicted stock composition estimates
head(ssdr[rownames(ssdr) %in% "logit_pred_Pi_prop", ])
# predicted stock specific abundance
head(ssdr[rownames(ssdr) %in% "log_pred_mu1_Pi", ])

## ----echo=TRUE----------------------------------------------------------------
plot_stacked_comp(x_var = "month_n", grouping_var = "agg", 
                  facet_var = "region", ssdr = m1$ssdr, comp_dat = comp_ex,
                  pred_dat = new_dat)

## ----echo=TRUE----------------------------------------------------------------
plot_ribbon_abund(x_var = "month_n", grouping_var = "agg", 
                 colour_var = "region", ssdr = m1$ssdr, comp_dat = comp_ex,
                 pred_dat = new_dat)

