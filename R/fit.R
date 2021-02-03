#' Fit model
#'
#' @param comp_dat Composition data
#' @param catch_dat Catch data
#' @param model_type Character specifying whether inputs for a composition-only
#'   or integrated model are being generated.
#' @param random_walk Logical (defaults to TRUE) specifying whether random 
#'   intercepts for year follow a random walk or are independent and identically 
#'   distributed. 
#' @param silent Logical
#' @param nlminb_loops How many times to run [stats::nlminb()] optimization.
#'   Sometimes restarting the optimizer at the previous best values aids
#'   convergence. If the maximum gradient is still too large,
#'   try increasing this to `2`.
#'
#' @return
#' Output from [TMB::sdreport()]
#' @export
#' 
#' @importFrom dplyr select distinct filter arrange mutate row_number mutate_if
#' @importFrom dplyr left_join group_by summarize group_split
#' @importFrom forcats fct_reorder
#' @importFrom tidyr complete nesting pivot_wider replace_na
#' @importFrom stats coef predict rnorm runif
#' @importFrom purrr map_dfr
#' @importFrom stats qnorm plogis
#'
#' @examples
#' m1 <- fit_stockseason(comp_dat = comp_ex, model_type = "composition")
#' theta <- as.list(m1, "Estimate")
#' theta_se <- as.list(m1, "Std. Error")
#' theta
#' 
#' m2 <- fit_stockseason(comp_dat = comp_ex, catch_dat = catch_ex, 
#'                       model_type = "integrated")

fit_stockseason <- function(comp_dat, catch_dat = NULL,
                            model_type = c("composition", "integrated"), 
                            random_walk = TRUE,
                            silent = TRUE, 
                            optim_fix_inits = TRUE,
                            nlminb_loops = 1) {
  
  model_type <- match.arg(model_type)
  #define model type
  if (model_type == "integrated" && is.null(catch_dat)) {
    stop("Cannot fit integrated model without catch data", call. = FALSE)
  }
  
  if (model_type == "integrated") {
    x <- gen_tmb(comp_dat = comp_dat, catch_dat = catch_dat, 
                 random_walk, model_type = model_type)
    random_ints <- c("z1_k", "z2_k")
    tmb_map1 <- c(
      x$tmb_map,
      list(
        z1_k = factor(rep(NA, length(x$pars$z1_k))),
        z2_k = factor(rep(NA, length(x$pars$z2_k))),
        log_sigma_zk1 = as.factor(NA),
        log_sigma_zk2 = as.factor(NA)
      )
    )
  } else if (model_type == "composition") {
    x <- gen_tmb(comp_dat = comp_dat, random_walk, model_type = model_type)
    random_ints <- "z2_k"
    tmb_map1 <- c(
      x$tmb_map,
      list(
        z2_k = factor(rep(NA, length(x$pars$z2_k))),
        log_sigma_zk2 = as.factor(NA)
      )
    )
  }
  
  #by default optimize for fixed effects first
  if (optim_fix_inits == TRUE) {
    if (!silent) cat("Optimizing for fixed effects\n")
    obj1 <- TMB::MakeADFun(
      data = c(list(model = "nb_dirchlet_1re"), x$data),
      parameters = x$pars,
      map = tmb_map1,
      DLL = "stockseasonr_TMBExports", silent = silent
    )
    opt1 <- stats::nlminb(obj1$par, obj1$fn, obj1$gr,
                          control = list(eval.max = 1e4, iter.max = 1e4)
    )
  }
  
  #use pars from optimizing for fixed effects as inits unless specified 
  #otherwise
  pars_in <- if (optim_fix_inits == TRUE) obj1$env$parList(opt1$par) else x$pars
  
  if (!silent) cat("Fitting fixed and random effects\n")
  obj <- TMB::MakeADFun(
    data = c(list(model = "nb_dirchlet_1re"), x$data),
    parameters = pars_in,#obj1$env$parList(opt1$par),
    map = x$tmb_map,
    random = random_ints,
    DLL = "stockseasonr_TMBExports", silent = silent
  )
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
  if (nlminb_loops > 1) {
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
      opt <- stats::nlminb(opt$par, obj$fn, obj$gr)
    }
  }
  TMB::sdreport(obj)
}
