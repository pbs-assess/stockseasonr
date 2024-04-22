
<!-- README.md is generated from README.Rmd. Please edit that file -->
  
# stockseasonr <a href='https://github.com/pbs-assess/stockseasonr'><img src='man/figures/stockseasonr-logo.png' align="right" style="height:139px;"/></a>
  
> Seasonal predictions of stock composition and abundance

<!-- badges: start -->
[![R-CMD-check](https://github.com/pbs-assess/stockseasonr/workflows/R-CMD-check/badge.svg)](https://github.com/pbs-assess/stockseasonr/actions)
<!-- badges: end -->

stockseasonr is an R package that fits models integrating composition and abundance data using Template Model Builder ([TMB](https://github.com/kaskr/adcomp)). 

## Installation

You can install the development version of stockseasonr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pbs-assess/stockseasonr")
```

## Overview

Composition data are common in ecological research. For example,
mixed-stock fisheries require biologists to account for stock
composition when describing distributions or abundance. stockseasonr is
intended to generate predictions of composition and group-specific
abundance. It borrows significant functionality from glmmTMB and sdmTMB.
Smooth terms can be incorporated using `s()` notation from mgcv and
random intercepts following the familiar lme4 syntax.

## Disclaimer

stockseasonr development is ongoing and sporadic! Please let us know if
you encounter bugs or if you have feature requests, please post in the
[issue tracker](https://github.com/pbs-assess/stockseasonr/issues).

## Basic Use

Spatially explicit, seasonal changes in stock composition and abundance
can be estimated using `fit_stockseasonr()`. `fit_stockseasonr()`
automatically generates the necessary data and parameter structures to
pass to TMB, fits, the model, and returns parameter estimates and
predictions.

We’ll generate inputs for an integrated version of the model using two
built-in datasets. The first is an example of composition data, which is
a subsample of available genetic stock identification data from the west
coast Vancouver Island commercial Chinook salmon fishery. `comp_ex` is a
long dataframe where each row includes the number of samples within a
sampling event (here a given day and sampling location) that correspond
to each stock aggregate. These data are non-integer because they
represent the sum, among individuals, of genetic stock assignment
probabilities. Each sampling event is spatially and temporally explicit
because it is assigned to a specific region, month, and year.

The second dataset is example catch and effort data from the same
fishery. `catch_ex` is a long dataframe where each row represents the
observed catch (individual pieces) for a given region, month, and year.
Note that both datasets have the same spatio-temporal strata, which is
necessary to fit integrated models.

We’ll assume that both stock composition and abundance are influenced by
season (month) and location (region). We’ll also account for interannual
variation with IID random intercepts. We can examine the model inputs
that are generated by `fit_stockseasonr()` without fitting the actual
model by setting `fit = FALSE`

``` r
dat_in <- fit_stockseasonr(
  abund_formula = catch ~ 1 + s(month_n, bs = "tp", k = 3, m = 2) + region + 
    (1 | year),
  abund_dat = catch_ex,
  comp_formula = agg ~ 1 + region + s(month_n, bs = "tp", k = 4, m = 2) +
    (1 | year),
  comp_dat = comp_ex,
  model = "integrated",
  fit = FALSE
)
```

`fit_stockseasonr()` reshapes or ‘widens’ the composition data, then
in-fills zero observations for a given group with small values (0.00001)
to ensure model convergence. It next uses the `mgcv` and `sdmTMB`
packages to create model matrices including seasonal smooths and random
intercepts for years. The result is a list of five elements:

- `tmb_data` list of data inputs (observations, model matrices, new
  dataset for prediction, etc.) passed to the TMB model
- `tmb_pars` list of initial parameter values passed to the TMB model
- `tmb_map` ‘mapped’ parameters that are not estimated
- `tmb_map` parameters that are estimated as random effects within TMB
- `wide_comp_dat` a wide version of the input composition data that
  includes infilled values

Next we’ll fit the model using `fit_stockseasonr()`, specifying the
input composition data, catch data, and model type (`negbin`,
`composition` or `integrated`) as above. This time we’ll also include a
new dataframe for which we would like to generate predictions. Note that
the model may take several seconds to converge.

``` r
new_dat <- expand.grid(
  month_n = seq(1, 12, by = 0.1),
  region = unique(comp_ex$region)
)

m1 <- fit_stockseasonr(
  abund_formula = catch ~ 1 + s(month_n, bs = "tp", k = 3, m = 2) + region +
    (1 | year),
  abund_dat = catch_ex,
  comp_formula = agg ~ 1 + region + s(month_n, bs = "tp", k = 4, m = 2) + 
    (1 | year),
  comp_dat = comp_ex,
  model = "integrated",
  pred_dat = new_dat,
  random_walk = TRUE,
  fit = TRUE
)
```

The object `m1` contains the model inputs described above, as well as a
`TMB::sdreport()` object (named `sdr`), which includes parameter
estimates for fixed effects.

``` r
theta <- as.list(m1$sdr, "Estimate")

# abundance estimates for regional effects
theta$b1_j
#> [1]  8.3078968 -0.2538266
#composition fixed effect estimates
theta$B2_jk
#>            [,1]       [,2]       [,3]       [,4]        [,5]       [,6]
#> [1,]  2.4356208  0.8879778  0.7411771  1.3658388  1.33768534  0.3426446
#> [2,] -0.5198176 -0.3716649  0.1682076 -1.0971537  0.79555707 -0.7513592
#> [3,] -0.6756046 -1.4294062 -0.4154278 -1.1222102 -0.05763734 -0.3148308
#> [4,]  0.7934574  0.2896438  0.4433785  0.2192940 -0.29938866 -0.2501597
#> [5,] -0.6468891 -1.4776849 -0.1728395 -0.7231197 -0.58343915 -0.5169122
#>             [,7]
#> [1,] -0.08099128
#> [2,] -0.20559936
#> [3,]  1.09572340
#> [4,]  1.41793680
#> [5,]  0.88044950
```

The `summary(TMB::sdreport))` (named `ssdr`) includes model derived
quantities, such as random effect estimates and integrated predictions
(i.e. estimates of stock-specific abundance). Note that all predictions
are provided in link space.

``` r
ssdr <- m1$ssdr
unique(rownames(ssdr))
#>  [1] "b1_j"               "ln_phi"             "bs"                
#>  [4] "ln_smooth_sigma"    "ln_sigma_re1"       "B2_jk"             
#>  [7] "ln_sigma_re2"       "b_smooth"           "re1"               
#> [10] "re2"                "log_pred_mu1"       "log_pred_Mu2"      
#> [13] "logit_pred_Pi_prop" "log_pred_mu1_Pi"

# predicted total abundance
head(ssdr[rownames(ssdr) %in% "log_pred_mu1", ])
#>              Estimate Std. Error
#> log_pred_mu1 6.041446  0.3293760
#> log_pred_mu1 6.154826  0.3258372
#> log_pred_mu1 6.268191  0.3225586
#> log_pred_mu1 6.381516  0.3195485
#> log_pred_mu1 6.494765  0.3168147
#> log_pred_mu1 6.607893  0.3143643
# predicted stock composition estimates
head(ssdr[rownames(ssdr) %in% "logit_pred_Pi_prop", ])
#>                      Estimate Std. Error
#> logit_pred_Pi_prop -0.9688726  0.1623525
#> logit_pred_Pi_prop -0.9425013  0.1546157
#> logit_pred_Pi_prop -0.9165025  0.1474827
#> logit_pred_Pi_prop -0.8908707  0.1410387
#> logit_pred_Pi_prop -0.8655975  0.1353599
#> logit_pred_Pi_prop -0.8406715  0.1305071
# predicted stock specific abundance
head(ssdr[rownames(ssdr) %in% "log_pred_mu1_Pi", ])
#>                 Estimate Std. Error
#> log_pred_mu1_Pi 4.750844  0.3497701
#> log_pred_mu1_Pi 4.883271  0.3443097
#> log_pred_mu1_Pi 5.015276  0.3393272
#> log_pred_mu1_Pi 5.146844  0.3348321
#> log_pred_mu1_Pi 5.277947  0.3308292
#> log_pred_mu1_Pi 5.408554  0.3273186
```

stockseasonr includes built-in plotting functions as example
visualizations for predictions. For example, predictions of seasonal
trends in mean stock composition can be visualized with
`plot_stacked_comp()`.

``` r
plot_stacked_comp(
  x_var = "month_n", grouping_var = "agg", facet_var = "region", 
  ssdr = m1$ssdr, comp_dat = comp_ex, pred_dat = new_dat
)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="60%" />

Estimates of stock-specific abundance can be used with the same inputs
and `plot_ribbon_abund()`.

``` r
plot_ribbon_abund(
  x_var = "month_n", grouping_var = "agg", colour_var = "region", 
  ssdr = m1$ssdr, comp_dat = comp_ex, pred_dat = new_dat
)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="60%" />
