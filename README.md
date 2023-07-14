
# stockseasonr <a href='https://github.com/pbs-assess/stockseasonr'><img src='man/figures/stockseasonr-logo.png' align="right" height="139" /></a>

> Seasonal predictions of stock composition and abundance

<!-- badges: start -->
[![R-CMD-check](https://github.com/pbs-assess/stockseasonr/workflows/R-CMD-check/badge.svg)](https://github.com/pbs-assess/stockseasonr/actions)
<!-- badges: end -->

## Installation

Assuming you have a [C++
compiler](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites)
installed, you can install *stockseasonr*:

``` r
# install.packages("remotes")
remotes::install_github("pbs-assess/stockseasonr")
```

## Summary 

Composition data are common in ecological research. For example, mixed-stock fisheries require biologists to account for stock composition when describing distributions or abundance. stockseasonr is intended to generate predictions of composition and group-specific abundance. It borrows significant functionality from glmmTMB and sdmTMB. Smooth terms can be incorporated using `s()` notation from mgcv and random intercepts following the familiar lme4 syntax.
