#' Prepare TMB inputs
#'
#' @param comp_dat Composition dataframe; each row represents a sampling event 
#'   (e.g. genetic samples collected from a given day and region).
#' @param catch_dat Catch data
#' @param model_type Character specifying whether inputs for a composition-only
#'   or integrated model are being generated.
#' @param random_walk Logical (defaults to TRUE) specifying whether random 
#'   intercepts for year follow a random walk or are independent and identically 
#'   distributed. 
#'
#' @return
#' List of inputs (data and parameter lists) to pass to TMB and dataframes to 
#'   use in subsequent plotting.
#' @export
#' 
#' @importFrom dplyr select distinct filter arrange mutate row_number mutate_if
#' @importFrom dplyr left_join group_by summarize group_split
#' @importFrom forcats fct_reorder
#' @importFrom tidyr complete nesting pivot_wider replace_na
#' @importFrom stats coef predict rnorm runif
#' @importFrom purrr map_dfr
#'
#' @examples
#' dat_list <- gen_tmb(comp_dat = comp_ex, catch_dat = catch_ex, 
#'                     model_type = "integrated")
#' dat_list$data #data inputs for TMB
#' dat_list$pars #parameter inputs for TMB

gen_tmb <- function(comp_dat, catch_dat = NULL,  
                    model_type = c("composition", "integrated"),
                    random_walk = TRUE) {
  
  model_type <- match.arg(model_type)
  
  if (model_type == "integrated" && is.null(catch_dat)) {
    stop("Cannot fit integrated model without catch data", call. = FALSE)
  }
  
  if (model_type == "integrated") {
    ## catch data
    yr_vec_catch <- as.numeric(as.factor(as.character(catch_dat$year))) - 1
    
    # generate model matrix based on GAM with spline type a function of months
    # retained
    months1 <- unique(catch_dat$month_n)
    spline_type <- if (max(months1) == 12) "cc" else "tp"
    n_knots <- if (max(months1) == 12) 4 else 3
    
    m1 <- mgcv::gam(catch ~ area + s(month_n, bs = spline_type, k = n_knots, 
                                     by = area) + offset,
                    data = catch_dat,
                    family = mgcv::nb()
    )
    fix_mm_catch <- predict(m1, type = "lpmatrix")
    
    # position of offset vector
    offset_pos <- grep("^offset$", colnames(fix_mm_catch))
    
    # generic data frame that is skeleton for both comp and catch predictions
    # conditional allows different regions to have different ranges of months
    pred_dat <- group_split(comp_dat, region) %>%
      map_dfr(., function(x) {
        expand.grid(
          month_n = seq(min(x$month_n),
                        max(x$month_n),
                        by = 0.1
          ),
          region = unique(x$region)
        )
      })
    
    
    # list of strata from composition dataset to retain in catch predictive 
    # model to ensure aggregate abundance can be calculated
    comp_strata <- unique(paste(pred_dat$month_n, pred_dat$region, sep = "_"))
    
    # make predictive model matrix including null values for effort
    pred_dat_catch <- expand.grid(
      month_n = seq(min(catch_dat$month_n),
                    max(catch_dat$month_n),
                    by = 0.1
      ),
      area = unique(catch_dat$area),
      offset = mean(catch_dat$offset)
    ) %>%
      left_join(.,
                catch_dat %>% select(region, area) %>% distinct(),
                by = "area"
      ) %>%
      mutate(strata = paste(month_n, region, sep = "_")) %>%
      filter(strata %in% comp_strata) %>%
      arrange(region) %>%
      # convoluted to ensure ordering is correct for key
      mutate(
        month = as.factor(month_n),
        order = row_number(),
        reg_month = paste(region, month_n, sep = "_"),
        reg_month_f = fct_reorder(factor(reg_month), order)
      ) %>%
      select(-order, -reg_month, -strata) %>%
      distinct()
    pred_mm_catch <- predict(m1, pred_dat_catch, type = "lpmatrix")
    
    # construct factor key for regional aggregates associated with areas
    grouping_vec <- as.numeric(pred_dat_catch$reg_month_f) - 1
    grouping_key <- unique(grouping_vec)
    ## composition data
    gsi_wide <- comp_dat %>%
      pivot_wider(., names_from = agg, values_from = agg_prob) %>%
      mutate_if(is.numeric, ~ replace_na(., 0.00001))
    
    obs_comp <- gsi_wide %>%
      select(-c(sample_id:nn)) %>%
      as.matrix()
    
    yr_vec_comp <- as.numeric(gsi_wide$year) - 1
    
    months2 <- unique(gsi_wide$month_n)
    spline_type <- if (max(months2) == 12) "cc" else "tp"
    n_knots <- if (max(months2) == 12) 4 else 3
    # response variable doesn't matter, since not fit
    m2 <- mgcv::gam(rep(0, length.out = nrow(gsi_wide)) ~
                      region + s(month_n, bs = spline_type, k = n_knots, 
                                 by = region),
                    data = gsi_wide
    )
    fix_mm_comp <- predict(m2, type = "lpmatrix")
    
    pred_mm_comp <- predict(m2, pred_dat, type = "lpmatrix")
    
    ## region-stock combinations with 0 observations to map
    comp_map <- comp_dat %>%
      group_by(region, agg) %>%
      complete(., region, nesting(agg)) %>%
      summarize(total_obs = sum(agg_prob), .groups = "drop") %>%
      filter(is.na(total_obs)) %>%
      select(agg, region)
    
    # input data
    data <- list(
      # abundance input data
      y1_i = catch_dat$catch,
      X1_ij = fix_mm_catch,
      factor1k_i = yr_vec_catch,
      nk1 = length(unique(yr_vec_catch)),
      X1_pred_ij = pred_mm_catch,
      pred_factor2k_h = grouping_vec,
      pred_factor2k_levels = grouping_key,
      # composition input data
      y2_ig = obs_comp,
      X2_ij = fix_mm_comp,
      factor2k_i = yr_vec_comp,
      nk2 = length(unique(yr_vec_comp)),
      X2_pred_ij = pred_mm_comp,
      # conditionals
      abundance_component = as.integer(1),
      random_walk = as.integer(random_walk)
    )
    
    # input parameter initial values
    pars <- list(
      # abundance parameters
      b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
      log_phi = log(1.5),
      z1_k = rep(0, length(unique(yr_vec_catch))),
      log_sigma_zk1 = log(0.25),
      # composition parameters
      b2_jg = matrix(0,#runif(ncol(fix_mm_comp) * ncol(obs_comp), 0.1, 0.9),
                     nrow = ncol(fix_mm_comp),
                     ncol = ncol(obs_comp)
      ),
      z2_k = rep(0, times = length(unique(yr_vec_comp))),
      log_sigma_zk2 = log(0.5)
    )
    
    # mapped values (not estimated)
    # offset for effort
    pars$b1_j[offset_pos] <- 1
    b1_j_map <- seq_along(pars$b1_j)
    b1_j_map[offset_pos] <- NA
    tmb_map <- list(b1_j = as.factor(b1_j_map))
    # offset for composition parameters
    if (!is.na(comp_map$agg[1])) {
      temp_betas <- pars$b2_jg
      for (i in seq_len(nrow(comp_map))) {
        offset_stock_pos <- grep(
          paste(comp_map$agg[i], collapse = "|"),
          colnames(obs_comp)
        )
        offset_region_pos <- grep(
          paste(comp_map$region[i], collapse = "|"),
          colnames(fix_mm_comp)
        )
        for (j in seq_len(ncol(fix_mm_comp))) {
          for (k in seq_len(ncol(obs_comp))) {
            if (j %in% offset_region_pos & k %in% offset_stock_pos) {
              # pars$b2_jg[j, k] <- 0.00001
              temp_betas[j, k] <- NA
            }
          }
        }
      }
      comp_map_pos <- which(is.na(as.vector(temp_betas)))
      b2_jg_map <- seq_along(pars$b2_jg)
      b2_jg_map[comp_map_pos] <- NA
      tmb_map <- c(tmb_map, list(b2_jg = as.factor(b2_jg_map)))
    }
    
    out_list <- list(
      "data" = data, "pars" = pars, "comp_wide" = gsi_wide,
      "pred_dat_catch" = pred_dat_catch, "pred_dat_comp" = pred_dat,
      "tmb_map" = tmb_map
    )
  }
  
  # stripped down version of above to only generate data for composition only 
  # version of mdoel
  if (model_type == "composition") {
    ## composition data
    gsi_wide <- comp_dat %>%
      pivot_wider(., names_from = agg, values_from = agg_prob) %>%
      mutate_if(is.numeric, ~ replace_na(., 0.00001))
    
    obs_comp <- gsi_wide %>%
      select(-c(sample_id:nn)) %>%
      as.matrix()
    
    yr_vec_comp <- as.numeric(gsi_wide$year) - 1
    
    # generic data frame that is skeleton for predictions
    # conditional allows different regions to have different ranges of months
    pred_dat <- group_split(comp_dat, region) %>%
      map_dfr(., function(x) {
        expand.grid(
          month_n = seq(min(x$month_n),
                        max(x$month_n),
                        by = 0.1
          ),
          region = unique(x$region)
        )
      })
    
    months2 <- unique(gsi_wide$month_n)
    spline_type <- if (max(months2) == 12) "cc" else "tp"
    n_knots <- if (max(months2) == 12) 4 else 3
    # response variable doesn't matter, since not fit
    m2 <- mgcv::gam(rep(0, length.out = nrow(gsi_wide)) ~
                      region + s(month_n, bs = spline_type, k = n_knots, 
                                 by = region),
                    data = gsi_wide
    )
    fix_mm_comp <- predict(m2, type = "lpmatrix")
    
    pred_mm_comp <- predict(m2, pred_dat, type = "lpmatrix")
    
    ## region-stock combinations with 0 observations to map
    comp_map <- comp_dat %>%
      group_by(region, agg) %>%
      complete(., region, nesting(agg)) %>%
      summarize(total_obs = sum(agg_prob), .groups = "drop") %>%
      filter(is.na(total_obs)) %>%
      select(agg, region)
    
    # input data
    data <- list(
      #fake catch data to feed to template
      y1_i = rep(1, length = 2),
      X1_ij = matrix(rep(1, length = 2), ncol = 1),
      factor1k_i = rep(as.integer(1), length = 2),
      nk1 = as.integer(1),
      X1_pred_ij = matrix(rep(1, length = 2), ncol = 1),
      pred_factor2k_h = rep(as.integer(1), length = 2),
      pred_factor2k_levels = rep(as.integer(1), length = 2),
      # comp data
      y2_ig = obs_comp,
      X2_ij = fix_mm_comp,
      factor2k_i = yr_vec_comp,
      nk2 = length(unique(yr_vec_comp)),
      X2_pred_ij = pred_mm_comp,
      abundance_component = as.integer(0),
      random_walk = as.integer(random_walk)
    )
    
    # input parameter initial values
    pars <- list(
      #fake abundance parameters to feed to template
      b1_j = rep(0, length = 2),
      log_phi = 0,
      z1_k = rep(0, length = 2),
      log_sigma_zk1 = 0,
      # composition parameters
      b2_jg = matrix(0,#runif(ncol(fix_mm_comp) * ncol(obs_comp), 0.1, 0.9),
                     nrow = ncol(fix_mm_comp),
                     ncol = ncol(obs_comp)
      ),
      z2_k = rep(0, times = length(unique(yr_vec_comp))),
      log_sigma_zk2 = log(0.5)
    )
    
    # mapped values (not estimated); includes all the fake abundance parameters
    # and levels with missing composition data
    tmb_map <- list(
      log_phi = factor(NA),
      log_sigma_zk1 = factor(NA),
      b1_j = factor(rep(NA, length(pars$b1_j))),
      z1_k = factor(rep(NA, length(pars$z1_k)))
    )
    # offset for composition parameters
    if (!is.na(comp_map$agg[1])) {
      temp_betas <- pars$b2_jg
      for (i in seq_len(nrow(comp_map))) {
        offset_stock_pos <- grep(
          paste(comp_map$agg[1], collapse = "|"),
          colnames(obs_comp)
        )
        offset_region_pos <- grep(
          paste(comp_map$region[1], collapse = "|"),
          colnames(fix_mm_comp)
        )
        for (j in seq_len(ncol(fix_mm_comp))) {
          for (k in seq_len(ncol(obs_comp))) {
            if (j %in% offset_region_pos & k %in% offset_stock_pos) {
              # pars$b2_jg[j, k] <- 0.00001
              temp_betas[j, k] <- NA
            }
          }
        }
      }
      comp_map_pos <- which(is.na(as.vector(temp_betas)))
      b2_jg_map <- seq_along(pars$b2_jg)
      b2_jg_map[comp_map_pos] <- NA
      tmb_map <- c(tmb_map, list(b2_jg = as.factor(b2_jg_map)))
    }
    
    out_list <- list(
      "data" = data, "pars" = pars, "comp_wide" = gsi_wide,
      "pred_dat_comp" = pred_dat, "tmb_map" = tmb_map
    )
  }
  return(out_list)
}
