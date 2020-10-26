% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make-model-inputs.R
\name{gen_tmb}
\alias{gen_tmb}
\title{Prep TMB inputs}
\usage{
gen_tmb(
  comp_dat,
  catch_dat = NULL,
  model_type = c("composition", "integrated"),
  random_walk = TRUE
)
}
\arguments{
\item{comp_dat}{Composition data}

\item{catch_dat}{Catch data}

\item{model_type}{Model type}

\item{random_walk}{Random walk or IID intercepts by year?}
}
\value{
List of inputs (data and parameter lists) to pass to TMB and dataframes to
use in subsequent plotting
}
\description{
Prep TMB inputs
}
\examples{
dat_list <- gen_tmb(catch_ex, comp_ex, model_type = "integrated")
dat_list$data #data inputs for TMB
dat_list$pars #parameter inputs for TMB
}