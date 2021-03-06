#' Example catch data
#'
#' Subset of catch and effort data from west-coast Vancouver Island commercial 
#' troll fishery for Chinook salmon.
#'
#' @format A data frame with 77 rows and six variables.
#' \describe{
#'   \item{region}{abbreviated catch region}
#'   \item{area}{Pacific Fishery Management Area (factor)}
#'   \item{month_n}{month (numeric)}
#'   \item{year}{year}
#'   \item{catch}{catch in pieces}
#'   \item{offset}{log(effort) used as an offset when fitting model (i.e. 
#'     coefficient fixed at one); effort units are vessel days}
#' }
#' @source Fisheries and Oceans Canada South Coast Area Office
"catch_ex"

#' Example composition data
#'
#' Subset of genetic stock composition data from west-coast Vancouver Island 
#' commercial troll fishery for Chinook salmon.
#'
#' @format A data frame with 1456 rows and seven variables.
#' \describe{
#'   \item{sample_id}{sampling event id; here unique to a year day-catch region-
#'     year}
#'   \item{region}{abbreviated catch region}
#'   \item{year}{year}
#'   \item{month_n}{month (numeric)}
#'   \item{agg}{stock aggregate; here corresponds to regional aggregates 
#'     identified by Pacific Salmon Treaty (e.g., PSD = Puget Sound)}
#'   \item{agg_prob}{summed individual stock assignment probabilities associated
#'     with a sampling event}
#'   \item{nn}{number of samples (individual fish) collected during a sampling 
#'     event}
#' }
#' @source Fisheries and Oceans Canada South Coast Area Office
"comp_ex"