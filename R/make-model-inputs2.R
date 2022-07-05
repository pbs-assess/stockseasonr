#' Prepare TMB inputs
#'
#' @param abund_formula Model formula (excluding RIs) for abundance component. 
#'   Penalized splines are possible via mgcv with s(). See examples and 
#'   details below.
#' @param comp_formula Model formula (excluding RIs) for composition component. 
#'   Splines are possible via mgcv, but are not penalized (i.e. k must be 
#'   specified). See examples and details below.
#' @param comp_knots Optional list specifying knots location for composition 
#'   formula.
#' @param abund_dat Dataframe used to fit abundance component.
#' @param comp_dat Dataframe used to fit composition component.
#' @param abund_rint Name of column in abund_dat containing random intercepts.
#' @param comp_rint Name of column in comp_dat containing random intercepts.
#' @param pred_dat Dataframe used to generate predictions. Applies to 
#'   both abundance and composition components in integrated model. If not 
#'   supplied then defaults to observations for "negbin" or "dirichlet", but 
#'   must be supplied for "integrated".
#' @param model Character specifying whether abundance ("negbin"), composition
#'   ("dirichlet"), or both ("integrated") model components should be estimated.
#' @param include_re_preds Logical specifying whether predictions should include
#'   random intercepts or not. In case of composition or integrated models this
#'   requires an alternative model with a multivariate normal distribution that 
#'   can be unstable with large numbers of parameters.
#'
#' @return
#' List of inputs (model_specs, tmb_data, tmb_pars, tmb_map and tmb_random) to 
#'   pass to TMB and dataframes to use in subsequent plotting.
#' @export
#' 
#' @importFrom dplyr mutate mutate_if
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom stats coef formula predict rnorm runif
#'
#' @examples

abund_formula = catch ~ s(month_n, k = 4) + region + offset
abund_dat = catch_ex 
abund_rint = "year"
model = "negbin"
include_re_preds = FALSE
abund_formula = catch ~ s(month_n, k = 4) + region + offset
abund_dat = catch_ex 
abund_rint = "year"
model = "negbin"

make_inputs <- function(abund_formula = NULL, comp_formula = NULL, 
                        comp_knots = NULL,
                        abund_dat = NULL, comp_dat = NULL,
                        abund_rint = NULL, comp_rint = NULL,
                        pred_dat = NULL,
                        model = c("negbin", "dirichlet", "integrated")) {
  
  # make sure necessary components are present
  if (model %in% c("dirichlet", "integrated")) {
    if (is.null(comp_dat)) {
      stop("Missing model inputs to fit integrated model")
    }
    if (is.null(comp_dat$prob)) {
      stop("Composition data not identified. Name vector in comp_dat 'prob' to 
         indicate proportions data.")
    }
  }
  
  # check predictions
  if (is.null(pred_dat)) {
    if (model == "integrated") stop("Must supply dataframe of predictions.")
    if (model == "negbin") pred_dat <- abund_dat
    if (model == "dirichlet") pred_dat <- comp_dat
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
                         pred_dat, type = "lpmatrix")
    sm_pred <- parse_smoothers(abund_formula, data = abund_dat, 
                               newdata = pred_dat)
    
    # identify response
    mf <- model.frame(formula_no_sm, abund_dat)
    y_i <- model.response(mf, "numeric")
    
    # random intercepts
    rfac1 <-  as.numeric(as.factor(abund_dat[[abund_rint]])) - 1
    
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
      pred_Xs = sm_pred$Xs
    )
    tmb_data <- c(tmb_data, abund_tmb_data)
    
    # adjust data when RI predictions being made
    if (include_re_preds == TRUE) {
      #vector of predicted random intercepts
      pred_rand_ints <- list(
        pred_rfac1 = as.numeric(pred_dat[[abund_rint]]) - 1)
      
      tmb_data <- c(tmb_data, pred_rand_ints)
    } 
    
    abund_tmb_pars <- list(
      b1_j = rep(0, ncol(X_ij)),
      ln_phi = log(1.5),
      bs = rep(0, ncol(sm$Xs)),
      ln_smooth_sigma = rnorm(length(sm$sm_dims), 0, 0.5), 
      b_smooth = rnorm(sum(sm$sm_dims), 0, 0.5),
      ln_sigma_a1 = log(0.25),
      a1 = rnorm(length(unique(rfac1)), 0, 0.5)
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
    comp_formula_split <- stringr::str_split(comp_formula, "~")
    comp_formula_new <- formula(paste("dummy", 
                                      comp_formula_split[[3]], sep = "~"))
    group_var <- comp_formula_split[[2]]
    
    # adjust input data to wide and convert observations to matrix
    comp_wide <- comp_dat %>%
      pivot_wider(names_from = as.name(group_var), 
                  values_from = prob) %>%
      mutate_if(is.numeric, ~ replace_na(., 0.00001)) %>% 
      #add dummy response variable
      mutate(dummy = 0)
    group_names <- unique(comp_dat[[group_var]])
    obs_comp <- comp_wide[ , group_names] %>%
      as.matrix()
    
    # dummy model
    dummy_comp <- mgcv::gam(comp_formula_new, data = comp_wide, 
                            knots = comp_knots)
    X2_ij <- predict(dummy_comp, type = "lpmatrix")
    pred_X2_ij <- predict(dummy_comp, pred_dat, type = "lpmatrix")
    
    # check to make sure predictive dataframes for composition and abundance
    # are same length
    if (model == "integrated"){
      if (nrow(pred_X2_ij) != nrow(pred_X_ij)) {
        stop("Dimensions of abundance and composition predictions are not
           compatible.")
      }}
    
    # NOTE: REPLACE WITH FLEXIBLE STRUCTURE WHEN NULL
    rfac2 <- as.numeric(as.factor(as.character(comp_wide[[comp_rint]]))) - 1
    n_rint2 <- length(unique(rfac2))
    
    #make composition tmb inputs 
    comp_tmb_data <- list(
      Y2_ik = obs_comp,
      X2_ij = X2_ij,
      rfac2 = rfac2,
      n_rfac2 = n_rint2,
      pred_X2_ij = pred_X2_ij
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
  }
  
  # combine model specs to pass as logicals to fitting function
  model_specs <- list(model = model, include_re_preds = include_re_preds)
  
  out_list <- list(model_specs = model_specs,
                   tmb_data = tmb_data, tmb_pars = tmb_pars, tmb_map = tmb_map, 
                   tmb_random = tmb_random)
  if (model %in% c("dirichlet", "integrated")) {
    out_list <- c(out_list, list(wide_comp_dat = comp_wide))
  }
  
  return(out_list)
}
