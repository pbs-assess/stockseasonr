#' Example catch data
#'
#' Subset of catch and effort data from west-coast Vancouver Island commercial 
#' troll fishery for Chinook salmon.
#'
#' @format A data frame with 77 rows and ten variables.
#' \describe{
#'   \item{region}{abbreviated catch region}
#'   \item{region_c}{unabbreviated catch region}
#'   \item{area}{Pacific Fishery Management Area (factor)}
#'   \item{area_n}{Pacific Fishery Management Area (numeric)}
#'   \item{month}{month (factor)}
#'   \item{month}{month (numeric)}
#'   \item{year}{year}
#'   \item{catch}{catch in pieces}
#'   \item{eff}{effort in vessel days}
#'   \item{offset}{log(effort) used as an offset when fitting model (i.e. 
#'     coefficient fixed at one)}
#' }
#' @source Fisheries and Oceans Canada South Coast Area Office
"catch_ex"

#' Example composition data
#'
#' Subset of genetic stock composition data from west-coast Vancouver Island 
#' commercial troll fishery for Chinook salmon.
#'
#' @format A data frame with 1456 rows and ten variables.
#' \describe{
#'   \item{sample_id}{sampling event id; here unique to a year day-catch region-
#'     year}
#'   \item{gear}{gear type (troll or recreational)}
#'   \item{region}{abbreviated catch region}
#'   \item{region_c}{unabbreviated catch region}
#'   \item{year}{year}
#'   \item{month}{month (factor)}
#'   \item{month}{month (numeric)}
#'   \item{agg}{stock aggregate; here corresponds to regional aggregates 
#'     identified by Pacific Salmon Treaty (e.g., PSD = Puget Sound)}
#'   \item{agg_prob}{summed individual stock assignment probabilities associated
#'     with a sampling event}
#'   \item{nn}{number of samples (individual fish) collected during a sampling 
#'     event}
#' }
#' @source Fisheries and Oceans Canada South Coast Area Office
"comp_ex"