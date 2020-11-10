#' @rawNamespace useDynLib(stockseasonr, .registration=TRUE); useDynLib(stockseasonr_TMBExports)
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".",
    "agg",
    "agg_prob",
    "area",
    "month",
    "month_n",
    "nn",
    "reg_month",
    "region",
    "sample_id",
    "strata",
    "total_obs",
    "link_abund_est",
    "link_abund_se",
    "link_abund_low",
    "link_abund_up",
    "pred_prob_est",
    "link_prob_est",
    "link_prob_se",
    "link_abund_up",
    "stock"
  ))
}