---
title: "Introduction to stockseasonr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to stockseasonr}
  %\VignetteEngine{knitr::knitr}
  \usepackage[utf8]{inputenc}
---





Spatially explicit, seasonal changes in stock composition and abundance can be estimated using `fit_stockseason()`. Although `fit_stockseason()` automatically generates the necessary data and parameters to pass to TMB, the `gen_tmb()` function can be used to generate these and other objects externally, which are useful for data visualization and to understand the data pre-processing.

We'll generate inputs for an integrated version of the model using two built-in datasets. The first is example composition data, which is a subsample of available genetic stock identification data from the west coast Vancouver Island commercial Chinook salmon fishery. `comp_ex` is a long dataframe where each row includes the number of samples within an event that correspond to a given stock aggregate (i.e., summed individual assignment probabilities). Each sampling event is spatially and temporally explicit because it is assigned to a given region, month, and year.

The second dataset is example catch and effort data from the same fishery, which are available at a finer scale (statistical areas) than the composition data. `catch_ex` is a long dataframe where each row represents the observed catch (individual pieces) and effort (boat days) for a given statistical area, month, and year. The offset variable is log(effort).













