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

