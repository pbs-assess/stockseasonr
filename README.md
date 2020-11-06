
# stockseasonr <a href='https://github.com/pbs-assess/stockseasonr'><img src='man/figures/stockseasonr-logo.png' align="right" height="139" /></a>

> Seasonal predictions of stock composition and abundance

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/pbs-assess/stockseasonr.svg?branch=master)](https://travis-ci.org/pbs-assess/stockseasonr)
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

*stockseasonr* is intended to generate seasonal predictions of stock composition and stock-specific abundance using spatially-stratified data. Currently interannual variation is accounted for using random intercepts. Model matrices are generated using *mgcv* and models are fit using *TMB*.