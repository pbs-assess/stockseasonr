## Updated version of make_inputs and fit_model intended to increase flexibility
#

# Note that two prediction datasets are necessary since abundance and composition
# data may be differently stratified

# Dev TO DO:
# 1) Make aggregate levels flexible (i.e. turn on/off based on input data)
# 2) Make RE intercepts flexible (i.e. turn on/off or remove MVN component)
# 3) Map 0 observations in composition component

library(tidyverse)
library(mgcv)
library(TMB)

# 
# abund_formula = catch ~ 0 + area +
#   s(month_n, bs = "tp", k = 3, by = area) +
#   s(month_n, by = year, bs = "tp", m = 1, k = 3) +
#   offset
# abund_dat = catch
# abund_rint = "year"
# pred_abund = pred_dat_catch
# comp_formula = coarse_agg ~  s(month_n, bs = "cc", k = 4#, by = reg
# )
# comp_dat = stock_comp
# comp_rint = "year"
# pred_comp = pred_dat_stock_comp
# model = "dirichlet"


## MAKE INPUTS  ----------------------------------------------------------------

make_inputs <- function(abund_formula = NULL, comp_formula = NULL, 
                        abund_dat = NULL, comp_dat = NULL,
                        abund_rint = NULL, comp_rint = NULL,
                        pred_abund = NULL, pred_comp = NULL,
                        model = c("negbin", "dirichlet", "integrated"),
                        include_re_preds = FALSE) {
  
  # make sure necessary components are present
  if (model %in% c("dirichlet", "integrated") & is.null(comp_dat) | 
      model %in% c("dirichlet", "integrated") & is.null(pred_comp)) {
    stop("Missing model inputs to fit integrated model")
  }
  
  if (model %in% c("dirichlet", "integrated") & is.null(comp_dat$prob)) {
    stop("Composition data not identified. Name vector in comp_dat 'prob' to 
         indicate proportions data.")
  }
  
  # initialize empty lists to fill with data and initial parameters
  tmb_data <- list()
  tmb_pars <- list()
  tmb_map <- list()
  tmb_random <- NULL
  
  if (model %in% c("negbin", "integrated")) {
    # generate inputs for negbin component
    formula_no_sm <- remove_s_and_t2(abund_formula)
    X_ij <- model.matrix(formula_no_sm, data = abund_dat)
    sm <- parse_smoothers(abund_formula, data = abund_dat)
    pred_X_ij <- predict(gam(formula_no_sm, data = abund_dat), 
                         pred_abund, type = "lpmatrix")
    pred_X_ij <- predict(gam(formula_no_sm, data = abund_dat), 
                         pred_abund, type = "lpmatrix")
    sm_pred <- parse_smoothers(abund_formula, data = abund_dat, 
                               newdata = pred_abund)
    
    # identify response
    mf <- model.frame(formula_no_sm, abund_dat)
    y_i <- model.response(mf, "numeric")
    
    # random intercepts
    rfac1 <-  as.numeric(as.factor(abund_dat[[abund_rint]])) - 1
    pred_rfac1 <- as.numeric(pred_abund[[abund_rint]]) - 1
    
    # aggregate key for pooling estimates
    # NOTE -- CHANGE TO MAKE FLEXIBLE (e.g conditional that changes key variable
    # to default stratification when missing)
    pred_rfac_agg <- as.numeric(pred_abund$key_var) - 1
    pred_rfac_agg_levels <- unique(pred_rfac_agg)
    
    # make abundance tmb inputs
    abund_tmb_data <- list(
      y1_i = y_i,
      X1_ij = X_ij,
      rfac1 = rfac1,
      n_rfac1 = length(unique(rfac1)),
      b_smooth_start = sm$b_smooth_start,
      Zs = sm$Zs, # optional smoother basis function matrices
      Xs = sm$Xs, # optional smoother linear effect matrix
      pred_X1_ij = pred_X_ij,
      pred_Zs = sm_pred$Zs,
      pred_Xs = sm_pred$Xs,
      pred_rfac1 = pred_rfac1,
      pred_rfac_agg = pred_rfac_agg,
      pred_rfac_agg_levels = pred_rfac_agg_levels
    )
    tmb_data <- c(tmb_data, abund_tmb_data)
    
    abund_tmb_pars <- list(
      b1_j = rep(0, ncol(X_ij)),
      ln_phi = log(1.5),
      bs = rep(0, ncol(sm$Xs)),
      ln_smooth_sigma = rnorm(length(sm$sm_dims), 0, 0.5), 
      b_smooth = rnorm(sum(sm$sm_dims), 0, 0.5),
      a1 = rnorm(length(unique(rfac1)), 0, 0.5), 
      ln_sigma_a1 = log(0.25)
    )
    tmb_pars <- c(tmb_pars, abund_tmb_pars)
    
    # mapped parameters
    offset_pos <- grep("^offset$", colnames(X_ij))
    tmb_pars$b1_j[offset_pos] <- 1
    b1_j_map <- seq_along(tmb_pars$b1_j)
    b1_j_map[offset_pos] <- NA
    tmb_map <- c(tmb_map, list(b1_j = as.factor(b1_j_map)))
    
    # random parameters
    tmb_random <- c(tmb_random, c("b_smooth", "a1"))
  }
  
  if (model %in% c("dirichlet", "integrated")) {
    # adjust composition model formula
    comp_formula_split <- str_split(comp_formula, "~")
    comp_formula_new <- formula(paste("dummy", 
                                      comp_formula_split[[3]], sep = "~"))
    group_var <- comp_formula_split[[2]]
    
    # adjust input data to wide and convert observations to matrix
    comp_wide <- comp_dat %>%
      pivot_wider(., names_from = as.name(group_var), values_from = prob) %>%
      mutate_if(is.numeric, ~ replace_na(., 0.00001)) %>% 
      #add dummy response variable
      mutate(dummy = 0)
    group_names <- unique(comp_dat[[group_var]])
    obs_comp <- comp_wide[ , group_names] %>%
      as.matrix()
    
    # dummy model
    dummy_comp <- mgcv::gam(comp_formula_new, data = comp_wide)
    
    
    X2_ij <- predict(dummy_comp, type = "lpmatrix")
    pred_X2_ij <- predict(dummy_comp, pred_comp, type = "lpmatrix")
    
    # NOTE: REPLACE WITH FLEXIBLE STRUCTURE WHEN NULL
    rfac2 <- as.numeric(as.factor(as.character(comp_wide[[comp_rint]]))) - 1
    n_rint2 <- length(unique(rfac2))
    
    #make composition tmb inputs 
    comp_tmb_data <- list(
      Y2_ik = obs_comp,
      X2_ij = X2_ij,
      rfac2 = rfac2,
      n_rfac2 = n_rint2,
      pred_X2_ij = pred_X2_ij#,
      # pred_rfac2 = pred_rfac2
    )
    tmb_data <- c(tmb_data, comp_tmb_data)
    
    comp_tmb_pars <- list(
      B2_jk = matrix(0,
                     nrow = ncol(X2_ij),
                     ncol = ncol(obs_comp)
      ),
      ln_sigma_A2 = log(0.25)
    )
    tmb_pars <- c(tmb_pars, comp_tmb_pars)
    
    # adjust data and paramets when RI predictions being made
    if (include_re_preds == TRUE) {
      #vector of predicted random intercepts
      pred_rand_ints <- list(
        pred_rfac2 = as.numeric(pred_comp[[comp_rint]]) - 1
      )
      # mvn matrix of REs
      mvn_rand_inits <- list(
        A2_hk = matrix(rnorm(n_rint2 * ncol(obs_comp), 0, 0.5), 
                       nrow = n_rint2,
                       ncol = ncol(obs_comp))
      )
      
      tmb_data <- c(tmb_data, pred_rand_ints)
      tmb_pars <- c(tmb_pars, mvn_rand_inits)
      tmb_random <- c(tmb_random, "A2_hk")
    } else if (include_re_preds == FALSE) {
      # vector of random intercepts
      rand_inits <- list(
        A2_h = rnorm(n_rint2, 0, 0.5)
        )
      
      tmb_pars <- c(tmb_pars, rand_inits)
      tmb_random <- c(tmb_random, "A2_h")
    }
  }
  
  out_list <- list(tmb_data = tmb_data, tmb_pars = tmb_pars, tmb_map = tmb_map, 
                   tmb_random = tmb_random)
  if (model %in% c("dirichlet", "integrated")) {
    out_list <- c(out_list, list(wide_comp_dat = comp_wide))
  }
  
  return(out_list)
}


## FIT MODELS  -----------------------------------------------------------------

# tmb_data = model_inputs$tmb_data; 
# tmb_pars = model_inputs$tmb_pars; 
# tmb_map = model_inputs$tmb_map; 
# tmb_random  = model_inputs$tmb_random;
# model = "dirichlet";
# fit_random = FALSE

fit_model <- function(tmb_data, tmb_pars, tmb_map = NULL, tmb_random = NULL,
                      model = c("negbin", "dirichlet", "integrated"),
                      fit_random = TRUE, ignore_fix = FALSE,
                      include_re_preds = FALSE) {
  
  if (model == "negbin") tmb_model <- "negbin_rsplines"
  # use MVN model if random effects predictions necessary
  if (model == "dirichlet") {
    tmb_model <- ifelse(include_re_preds == FALSE,
                        "dirichlet_ri",
                        "dirichlet_mvn")
  }
  if (model == "integrated") {
    tmb_model <- ifelse(include_re_preds == FALSE,
                        "negbin_rsplines_dirichlet_ri",
                        "negbin_rsplines_dirichlet_mvn")
  }
  ## fit fixed effects only 
  # map random effects
  if (ignore_fix == FALSE) {
    new_map_list <- tmb_pars[names(tmb_pars) %in% tmb_random]
    tmb_map_random <- c(
      tmb_map,
      map(new_map_list, function (x) factor(rep(NA, length(x))))
    )
    # fit
    obj1 <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_pars,
      map = tmb_map_random,
      DLL = tmb_model
    )
    opt1 <- stats::nlminb(obj1$par, obj1$fn, obj1$gr,
                          control = list(eval.max = 1e4, iter.max = 1e4)
    )
    sdr <- sdreport(obj1)
  }
  
  ## fit with random effects 
  if (fit_random) {
    
    # pass parameter inits from above unless specified otherwise
    if (ignore_fix == TRUE) {
      pars_in <- tmb_pars   
    } else {
      pars_in <- obj1$env$parList(opt1$par)
    } 
    
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = pars_in,#obj1$env$parList(opt1$par),
      map = tmb_map,
      random = tmb_random,
      DLL = tmb_model
    )
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
    nlminb_loops = 2
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
      opt <- stats::nlminb(opt$par, obj$fn, obj$gr)
    }
    
    sdr <- sdreport(obj)
  }
  
  # derived quantities
  ssdr <- summary(sdr)
  
  return(list(sdr = sdr, ssdr = ssdr))
}
