#' Fit model
#'
#' @param tmb_data List of tmb input data generated in 
#'   `stockseasonr::make_inputs()`.
#' @param tmb_pars List of tmb input parameters generated in 
#'   `stockseasonr::make_inputs()`.
#' @param tmb_map Optional list of tmb input parameters to be mapped (i.e. 
#'   ignored) generated in [stockseasonr::make_inputs()].
#' @param tmb_random Optional list of tmb input parameters treated as random
#'   effects (includes smooths and intercepts) generated in 
#'   `stockseasonr::make_inputs()`.
#' @param fit_random Logical (defaults to FALSE) specifying whether random 
#'   effects should be estimated.
#' @param ignore_fix Logical (defaults to FALSE) specifying whether a 
#'   fixed-effects only version of the model should be fit first to generate 
#'   initial parameters.
#' @param nlminb_loops How many times to run [stats::nlminb()] optimization.
#'   Sometimes restarting the optimizer at the previous best values aids
#'   convergence. If the maximum gradient is still too large,
#'   try increasing this to `2`.
#' @param newton_loops How many Newton optimization steps to try with 
#'   [stats::optimHess()] after running [stats::nlminb()]. Sometimes aids 
#'   convergence.
#'
#' @return
#' Output from [TMB::sdreport()] and [summary(TMB::sdreport())]
#' @export
#' 
#' @examples

fit_stockseasonr <- function(tmb_data, tmb_pars, tmb_map = NULL, 
                             tmb_random = NULL, fit_random = FALSE, 
                             ignore_fix = FALSE, model_specs,
                             nlminb_loops = 1L, newton_loops = 0L) {
  
  if (model_specs$model == "negbin") tmb_model <- "negbin_rsplines"
  # use MVN model if random effects predictions necessary
  if (model_specs$model == "dirichlet") {
    tmb_model <- ifelse(model_specs$include_re_preds == FALSE,
                        "dirichlet_ri",
                        "dirichlet_mvn")
  }
  if (model_specs$model == "integrated") {
    tmb_model <- ifelse(model_specs$include_re_preds == FALSE,
                        "negbin_rsplines_dirichlet_ri",
                        "negbin_rsplines_dirichlet_mvn")
  }
  
  ## fit fixed effects only 
  # map random effects
  if (fit_random == FALSE | ignore_fix == FALSE) {
    new_map_list <- tmb_pars[names(tmb_pars) %in% tmb_random]
    tmb_map_random <- c(
      tmb_map,
      map(new_map_list, function (x) factor(rep(NA, length(x))))
    )
    # fit
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_pars,
      map = tmb_map_random,
      DLL = tmb_model
    )
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                         control = list(eval.max = 1e4, iter.max = 1e4)
    )
  }
  
  ## fit with random effects 
  if (fit_random) {
    # pass parameter inits from above unless specified otherwise
    if (ignore_fix == TRUE) {
      pars_in <- tmb_pars   
    } else {
      pars_in <- obj$env$parList(opt$par)
    } 
    
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = pars_in,
      map = tmb_map,
      random = tmb_random,
      DLL = tmb_model
    )
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
  }
  
  if (nlminb_loops > 1) {
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
      temp <- opt[c("iterations", "evaluations")]
      opt <- stats::nlminb(
        start = opt$par, objective = obj$fn, gradient = obj$gr)
      opt[["iterations"]] <- opt[["iterations"]] + temp[["iterations"]]
      opt[["evaluations"]] <- opt[["evaluations"]] + temp[["evaluations"]]
    }
  }
  
  if (newton_loops > 0) {
    for (i in seq_len(newton_loops)) {
      g <- as.numeric(obj$gr(opt$par))
      h <- stats::optimHess(opt$par, fn = obj$fn, gr = obj$gr)
      opt$par <- opt$par - solve(h, g)
      opt$objective <- obj$fn(opt$par)
    }
  }
  
  sdr <- sdreport(obj)
  
  # derived quantities
  ssdr <- summary(sdr)
  
  return(list(sdr = sdr, ssdr = ssdr))
}
