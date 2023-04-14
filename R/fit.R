#' Fit model
#'
#' @param abund_formula Model formula for abundance component of model (fit with
#'   negative binomial distribution (nbinom2)). Random intercepts (IID or RW) 
#'   are possible using lme4 syntax, e.g., + (1 | g) where g is a column of 
#'   class factor representing groups. Penalized splines are possible via mgcv 
#'   with s() (te not currently supported).
#' @param comp_formula Model formula for composition component of model (fit with
#'   dirichlet distribution). Random intercepts (IID or RW) are possible using 
#'   lme4 syntax, e.g., + (1 | g) where g is a column of class factor 
#'   representing groups. Splines are possible via mgcv with s(), but unlike 
#'   abundance component these are unpenalized so knots must be specified with 
#'   k=.
#' @param abund_dat Dataframe containing abundance response (e.g. catch) and 
#'   covariates.
#' @param comp_dat Dataframe containing composition response in long format (i.e
#'   column of class factor) and covariates.
#' @param abund_offset A numeric vector representing the model offset. Usually a
#'   log transformed variable. *Not included in any prediction.*
#' @param pred_dat Optional dataframe used to generate predictions. Note that 
#'   in integrated model this applies to both abundance and composition 
#'   components. 
#' @param model String specifying whether abundance (negbin), composition 
#'   (dirichlet), or integrated model is fit. 
#' @param random_walk Logical specifying whether random intercepts are IID or 
#'   random walk. 
#' @param fit Logical specifying whether model should be fit in TMB. If FALSE 
#'   list of model inputs is returned. 
#' @param silent Silent or include optimization details? Helpful to set to FALSE
#'   for models that take a while to fit.
#' @param nlminb_loops How many times to run [stats::nlminb()] optimization.
#'   Sometimes restarting the optimizer at the previous best values aids
#'   convergence. If the maximum gradient is still too large,
#'   try increasing this to `2`.
#' @param newton_loops Alternative Newton optimizer for [stats::optimHess()]. 
#'   Sometimes aids convergence.
#'
#' @return
#' List including model inputs as well as output from [TMB::sdreport()] 
#' (if fit = TRUE).
#' @export
#' 
#' @importFrom dplyr select mutate mutate_if
#' @importFrom tidyr pivot_wider replace_na
#' 
#' @examples 
#' dum_pred <- expand.grid(month_n = seq(1, 12, by = 0.1),
#'                         region = unique(comp_ex$region))
#'
#' # composition model with hierarchical random walk
#' m1 <- fit_stockseasonr(comp_formula = agg ~ 1 + region + 
#'                         s(month_n, bs = "tp", k = 4, m = 2) + 
#'                         (1 | year),
#'                        comp_dat = comp_ex,
#'                        pred_dat = dum_pred,
#'                        model = "dirichlet",
#'                        random_walk = TRUE,
#'                        fit = TRUE,
#'                        newton_loops = 1)
#'                                  
#' # group-specific integrated model with hierarchical random walk
#' m2 <- fit_stockseasonr(abund_formula = catch ~ 1 +
#'                        s(month_n, bs = "tp", k = 3, m = 2) +
#'                          region +
#'                          (1 | year),
#'                        abund_dat = catch_ex,
#'                        abund_offset = "offset",
#'                        comp_formula = agg ~ 1 + region + 
#'                          s(month_n, bs = "tp", k = 4, m = 2) + 
#'                          (1 | year),
#'                        comp_dat = comp_ex,
#'                        pred_dat = dum_pred,
#'                        model = "integrated",
#'                        random_walk = TRUE,
#'                        fit = TRUE,
#'                        silent = FALSE,
#'                        newton_loops = 1)
#'                                  

fit_stockseasonr <- function(abund_formula = NULL, comp_formula = NULL, 
                             abund_dat = NULL, comp_dat = NULL,
                             abund_offset = NULL,
                             pred_dat = NULL,
                             model = c("negbin", "dirichlet", "integrated"),
                             random_walk = FALSE,
                             fit = TRUE,
                             silent = TRUE,
                             nlminb_loops = 1L, newton_loops = 0L) {
  
  # make sure necessary components are present
  if (model != "negbin") {
    if (is.null(comp_dat)) {
      stop("Missing model inputs to fit integrated model")
    }
    if (is.null(comp_dat$prob)) {
      stop("Composition data not identified. Name vector in comp_dat 'prob' to 
       indicate proportions data.")
    }
  }
  
  # initialize empty lists to fill with data and initial parameters
  tmb_data <- list()
  tmb_pars <- list()
  tmb_map <- list()
  tmb_random <- NULL
  
  # predictions present
  has_preds <- ifelse(is.null(pred_dat), 0, 1)
  
  # generate inputs for negbin component
  if (model %in% c("integrated", "negbin")) {
    # random intercepts (use sdmTMB to generate proper structure)
    sdmTMB_dummy <- sdmTMB::sdmTMB(
      abund_formula,
      data = abund_dat,
      spatial = "off",
      do_fit = FALSE
    )
    
    # generate offset separately
    if (is.null(abund_offset)) {
      offset <- rep(0, nrow(abund_dat))
    } else {
      offset <- abund_dat[[abund_offset]]
    }
    
    # smooths present (conditional for determining input structure)
    has_smooths <- as.integer(sdmTMB_dummy$tmb_data$has_smooths)
    
    if (has_preds == 1) {
      resp <- attr(stats::terms(abund_formula), which = "variables")[[2]] %>% 
        as.character()
      pred_dat[resp] <- 0
      
      sdmTMB_dummy_p <- sdmTMB::sdmTMB(
        # exclude REs from predictive dataset
        glmmTMB::splitForm(abund_formula)$fixedFormula,
        data = pred_dat,
        spatial = "off",
        do_fit = FALSE
      )
    } else {
      # copy original template inputs as placeholders (not used in cpp due to
      # conditional)
      sdmTMB_dummy_p <- sdmTMB_dummy
    }
    
    if (random_walk == TRUE) {
      rw_index <- rw_index_foo(sdmTMB_dummy$tmb_data$RE_indexes)
    } else {
      rw_index <- 0
    }
    
    # make abundance tmb inputs
    abund_tmb_data <- list(
      y1_i = sdmTMB_dummy$tmb_data$y_i %>% as.numeric(),
      X1_ij = sdmTMB_dummy$tmb_data$X_ij[[1]],
      re_index1 = sdmTMB_dummy$tmb_data$RE_indexes,
      ln_sigma_re_index1 = sdmTMB_dummy$tmb_data$ln_tau_G_index,
      nobs_re1 = sdmTMB_dummy$tmb_data$nobs_RE, # number of random intercepts
      Zs = sdmTMB_dummy$tmb_data$Zs, # optional smoother basis function matrices
      Xs = sdmTMB_dummy$tmb_data$Xs, # optional smoother linear effect matrix
      offset_i = offset,
      # random_walk = as.numeric(random_walk),
      rw_index1 = rw_index,
      has_smooths = has_smooths,
      b_smooth_start = sdmTMB_dummy$tmb_data$b_smooth_start,
      # has_preds = has_preds,
      pred_X1_ij = sdmTMB_dummy_p$tmb_data$X_ij[[1]],
      pred_Zs = sdmTMB_dummy_p$tmb_data$Zs,
      pred_Xs = sdmTMB_dummy_p$tmb_data$Xs
    )
    tmb_data <- c(tmb_data, abund_tmb_data)
    
    # TODO: adjust data when RI predictions being made
    # if (include_re_preds == TRUE) {
    #   #vector of predicted random intercepts
    #   pred_rand_ints <- list(
    #     pred_rfac1 = as.numeric(pred_dat[[abund_rint]]) - 1)
    # 
    #   tmb_data <- c(tmb_data, pred_rand_ints)
    # }
    
    abund_tmb_pars <- list(
      b1_j = rep(0, ncol(abund_tmb_data$X1_ij)),
      ln_phi = log(1.5),
      bs = if (has_smooths) {
        sdmTMB_dummy$tmb_params$bs %>% as.vector() 
      } else {
        0
      },
      ln_smooth_sigma = if (has_smooths) {
        sdmTMB_dummy$tmb_params$ln_smooth_sigma %>% as.vector()
      } else {
        0
      },
      b_smooth = if (has_smooths) {
        rep(0, nrow(sdmTMB_dummy$tmb_params$b_smooth))
      } else {
        0
      },
      #for some reason cannot fix as zeros when ignore_fix = FALSE; use sdmTMB
      #inputs for length but redefine with rnorm
      re1 = stats::rnorm(n = length(sdmTMB_dummy$tmb_params$RE %>% as.vector()),
                         0,
                         0.5),
      ln_sigma_re1 = sdmTMB_dummy$tmb_params$ln_tau_G %>% as.vector()
    )
    tmb_pars <- c(tmb_pars, abund_tmb_pars)
    
    # map parameters unless necessary
    if (abund_tmb_data$nobs_re1[[1]] > 0) {
      tmb_random <- c(tmb_random, "re1")
    } else {
      re_map <- map_foo(x = c("re1", "ln_sigma_re1"),
                        tmb_pars = tmb_pars)
      tmb_map <- c(tmb_map, 
                   re_map)
    }
    if (abund_tmb_data$has_smooths > 0) {
      tmb_random <- c(tmb_random, "b_smooth")
    } else {
      smooth_map <- map_foo(x = c("b_smooth", "ln_smooth_sigma", "bs"),
                            tmb_pars = tmb_pars)
      tmb_map <- c(tmb_map, 
                   smooth_map)
    }
  } 
  
  # as above but for composition component of model
  if (model %in% c("integrated", "dirichlet")) {
    # identify grouping variable for composition component (e.g. vector of 
    # stock names)
    comp_formula_new <- stats::formula(
      paste("dummy", deparse(comp_formula[[3]]), sep = "~")
    )
    group_var <- deparse(comp_formula[[2]])
    
    # adjust input data to wide and convert observations to matrix
    comp_wide <- comp_dat %>%
      pivot_wider(names_from = as.name(group_var), 
                  values_from = .data$prob) %>%
      mutate_if(is.numeric, ~ replace_na(., 0.00001)) %>% 
      #add dummy response variable
      mutate(dummy = 0)
    group_names <- unique(comp_dat[[group_var]])
    obs_comp <- comp_wide[ , group_names] %>%
      as.matrix()
    
    # generate main effects model matrix 
    # NOTE can't use sdmTMB because penalized smooths are not readily compatible
    # with multivariate response data 
    fixed_formula <- glmmTMB::splitForm(comp_formula_new)$fixedFormula
    dummy_comp <- mgcv::gam(fixed_formula, data = comp_wide)
    X2_ij <- stats::predict(dummy_comp, type = "lpmatrix")
    
    # generate RI inputs using sdmTMB 
    sdmTMB_dummy <- sdmTMB::sdmTMB(
      comp_formula_new,
      data = comp_wide,
      spatial = "off",
      do_fit = FALSE
    )
    
    # NOTE dummy sdmTMB not necessary because predictions based on main effects
    # only, but will be if random int predictions are incorporated (look at 
    # abundance version for reference)
    if (has_preds == 1) { 
      pred_X2_ij <- stats::predict(dummy_comp, pred_dat, type = "lpmatrix")
    } else {
      pred_X2_ij <- X2_ij
    } 
    
    # generate vector equal to length re2 indicating which random effects should
    # be estimated as random walk
    if (random_walk == TRUE) {
      rw_index <- rw_index_foo(sdmTMB_dummy$tmb_data$RE_indexes)
    } else {
      rw_index <- 0
    }
    
    # make composition tmb inputs
    comp_tmb_data <- list(
      Y2_ik = obs_comp,
      X2_ij = X2_ij,
      re_index2 = sdmTMB_dummy$tmb_data$RE_indexes,
      ln_sigma_re_index2 = sdmTMB_dummy$tmb_data$ln_tau_G_index,
      nobs_re2 = sdmTMB_dummy$tmb_data$nobs_RE, # number of random intercepts
      random_walk = as.numeric(random_walk),
      rw_index2 = rw_index,
      has_preds = has_preds,
      pred_X2_ij = pred_X2_ij
    )
    tmb_data <- c(tmb_data, comp_tmb_data)
    
    comp_tmb_pars <- list(
      B2_jk = matrix(0,
                     nrow = ncol(X2_ij),
                     ncol = ncol(obs_comp)
      ),
      #for some reason cannot fix as zeros when ignore_fix = FALSE; use sdmTMB
      #inputs for length but redefine with rnorm
      re2 = stats::rnorm(
        n = length(sdmTMB_dummy$tmb_params$RE %>% as.vector()),
        0,
        0.5
      ),
      ln_sigma_re2 = stats::rnorm(
        n = length(sdmTMB_dummy$tmb_params$ln_tau_G %>% 
                     as.vector()),
        0,
        0.5
      )
    )
    tmb_pars <- c(tmb_pars, comp_tmb_pars)
    
    # map random intercepts unless necessary
    if (comp_tmb_data$nobs_re2[[1]] > 0) {
      tmb_random <- c(tmb_random, "re2")
    } else {
      re_map <- map_foo(x = c("re2", "ln_sigma_re2"),
                        tmb_pars = tmb_pars)
      tmb_map <- c(tmb_map, 
                   re_map)
    }
    
    # TODO: adjust data and parameters when RI predictions being made
    # if (include_re_preds == TRUE) {
    #   #vector of predicted random intercepts
    #   #only added for dirichlet because generated in neg bin component for 
    #   #integrated model
    #   if (model == "dirichlet") {
    #     pred_rand_ints <- list(
    #       pred_rfac1 = as.numeric(pred_dat[[comp_rint]]) - 1
    #     )
    #     tmb_data <- c(tmb_data, pred_rand_ints)
    #   }
    #   # mvn matrix of REs
    #   mvn_rand_inits <- list(
    #     A2_hk = matrix(rnorm(n_rint2 * ncol(obs_comp), 0, 0.5), 
    #                    nrow = n_rint2,
    #                    ncol = ncol(obs_comp))
    #   )
    #   
    #   tmb_pars <- c(tmb_pars, mvn_rand_inits)
    #   tmb_random <- c(tmb_random, "A2_hk")
    # } else if (include_re_preds == FALSE) {
    #   # vector of random intercepts
    #   rand_inits <- list(
    #     A2_h = rep(0, n_rint2)
    #   )
    #   
    #   tmb_pars <- c(tmb_pars, rand_inits)
    #   tmb_random <- c(tmb_random, "A2_h")
    # }
  }
  
  # dimensions of predictions identical for both model components?
  if (model == "integrated" & has_preds == 1) {
    if(nrow(comp_tmb_data$pred_X2_ij) != nrow(abund_tmb_data$pred_X1_ij)) {
      stop("Composition and abundance predictions are not symmetrical.")
    }
  }
  
  
  # number of predictions 
  n_predX <- ifelse(model %in% c("dirichlet", "integrated"),
                    nrow(comp_tmb_data$pred_X2_ij),
                    nrow(abund_tmb_data$pred_X1_ij))
  
  # add data components that are present for both model types
  shared_tmb_data <- list(
    random_walk = as.numeric(random_walk),
    has_preds = has_preds,
    n_predX = n_predX
  )
  tmb_data <- c(tmb_data, shared_tmb_data)
  
  # combine model specs to pass as logicals to fitting function
  # 
  out_list <- list(tmb_data = tmb_data, tmb_pars = tmb_pars, tmb_map = tmb_map, 
                   tmb_random = tmb_random)
  
  if (model %in% c("dirichlet", "integrated")) {
    out_list <- c(out_list, list(wide_comp_dat = comp_wide))
  }
  
  # fit
  if (fit == TRUE) {
    if (model == "negbin") tmb_model <- "negbin"
    # use MVN model if random effects predictions necessary
    if (model == "dirichlet") {
      tmb_model <- #ifelse(include_re_preds == FALSE,
        "dirichlet"#,
      # "dirichlet_mvn")
    }
    if (model == "integrated") {
      tmb_model <- #ifelse(include_re_preds == FALSE,
        "integrated"#,
      # "negbin_rsplines_dirichlet_mvn")
    }
    
    obj <- TMB::MakeADFun(
      data = c(list(model = tmb_model), tmb_data), 
      parameters = tmb_pars, 
      map = tmb_map,
      random = tmb_random,
      DLL = "stockseasonr_TMBExports",
      silent = silent
    )
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
    
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
    
    sdr <- TMB::sdreport(obj)
    ssdr <- summary(sdr)
    
    out_list <- c(out_list, list(sdr = sdr, ssdr = ssdr))
  }
  
  return(out_list)
}
