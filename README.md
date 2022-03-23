
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

Mixed-stock fisheries are common, requiring biologists to account for stock composition when describing distributions or abundance. In the case of migratory populations (e.g., salmon) the contribution of different stocks can change seasonally. stockseasonr is intended to generate seasonal predictions of stock composition and stock-specific abundance using spatially-stratified data. Currently interannual variation is accounted for using random intercepts. Briefly, it uses splines implemented in mgcv to account for seasonal effects, then fits a model in TMB assuming composition data follow a Dirichlet-multinomial distribution and catch data follow a negative binomial distribution. 

**NOTE**
*stockseasonr is under active development. As currently written it is well suited to analyses similar to [Freshwater et al. 2021](https://peerj.com/articles/11163/), but has not been extensively tested with other datasets. We are currently working to make its functions more flexible.*