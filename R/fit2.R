
## FIT MODELS  -----------------------------------------------------------------

# tmb_data = model_inputs_ri$tmb_data;
# tmb_pars = model_inputs_ri$tmb_pars;
# tmb_map = model_inputs_ri$tmb_map;
# tmb_random  = model_inputs_ri$tmb_random;
# fit_random = TRUE;
# ignore_fix = FALSE;
# model_specs = list(model = "dirichlet",
#                    include_re_preds = FALSE)
# nlminb_loops = 2
# newton_loops = 2

fit_model <- function(tmb_data, tmb_pars, tmb_map = NULL, tmb_random = NULL,
                      model = c("negbin", "dirichlet", "integrated"),
                      fit_random = TRUE, ignore_fix = FALSE,
                      include_re_preds = FALSE) {
  
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
