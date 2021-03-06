% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_stacked_comp.R
\name{plot_stacked_comp}
\alias{plot_stacked_comp}
\title{Stacked composition plot}
\usage{
plot_stacked_comp(comp_dat, pred_dat_comp, mod_fit)
}
\arguments{
\item{comp_dat}{Composition dataframe; each row represents a sampling event
(e.g. genetic samples collected from a given day and region)}

\item{pred_dat_comp}{Dataframe of fixed effects composition predictions (i.e.
months and regions); generated using stockseasonr::gen_tmb()}

\item{mod_fit}{Parameter estimates generated by TMB::sdreport() and returned
by stockseasonr::fit_stockseason()}
}
\value{
Stacked ribbon plot representing seasonal trends in mean stock composition
and faceted by catch region.
}
\description{
Stacked composition plot
}
\examples{
in1 <- gen_tmb(comp_dat = comp_ex, model_type = "composition")
# model fitting will take several seconds
m1 <- fit_stockseason(comp_dat = comp_ex, model_type = "composition")
plot_stacked_comp(comp_dat = comp_ex, pred_dat_comp = in1$pred_dat_comp,
                  mod_fit = m1)

}
