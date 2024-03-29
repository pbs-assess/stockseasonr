#' Group-specific abundance plot
#'
#' @param x_var Explanatory variable within pred_dat to be plotted on x-axis 
#'   (must be continuous). 
#' @param grouping_var Character or factor variable within pred_dat that
#'   was used to generate composition estimates.
#' @param colour_var Optional character or factor variable within pred_dat 
#'   corresponding to categorical main effects in fitted model.
#' @param ssdr Parameter estimates generated by 
#'   stockseasonr::fit_stockseasonr()$ssdr.
#' @param comp_dat Composition dataframe; each row represents a sampling event 
#'   (e.g. genetic samples collected from a given day and region)
#' @param pred_dat Dataframe of fixed effects composition predictions
#'
#' @return
#' Ribbon plot representing fixed effects predictions of trends in group-
#'   specific abundance, faceted by grouping variable, with 95% confidence 
#'   intervals.
#' @export
#' 
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_ribbon  
#' @importFrom ggplot2 scale_colour_brewer labs theme scale_x_continuous 
#' @importFrom ggplot2 facet_wrap coord_cartesian element_text unit 
#' @importFrom ggplot2 theme_classic scale_fill_brewer
#' @importFrom rlang .data
#' 
#' @examples
#' dum_pred <- expand.grid(
#'   month_n = seq(1, 12, by = 0.1),
#'   region = unique(comp_ex$region)
#' )
#' # model fitting will take several seconds
#' m2 <- fit_stockseasonr(
#'   abund_formula = catch ~ 1 +
#'     s(month_n, bs = "tp", k = 3, m = 2) +
#'     region +
#'     (1 | year),
#'   abund_dat = catch_ex,
#'   abund_offset = "offset",
#'   comp_formula = agg ~ 1 + region + 
#'     s(month_n, bs = "tp", k = 4, m = 2) + 
#'     (1 | year),
#'   comp_dat = comp_ex,
#'   pred_dat = dum_pred,
#'   model = "integrated",
#'   random_walk = TRUE,
#'   fit = TRUE,
#'   silent = FALSE,
#'   nlminb_loops = 2, newton_loops = 1
#' )
#' plot_ribbon_abund(x_var = "month_n", grouping_var = "agg", 
#'                   colour_var = "region", ssdr = m2$ssdr, comp_dat = comp_ex,
#'                   pred_dat = dum_pred)
#' 

plot_ribbon_abund <- function(x_var, grouping_var, colour_var = NULL,
                              ssdr, comp_dat, pred_dat) {
  
  comp_abund_pred <- ssdr[rownames(ssdr) %in% "log_pred_mu1_Pi", ]
  group_names <- unique(comp_dat[[grouping_var]])
  dum <- purrr::map_dfr(seq_along(group_names), ~ pred_dat)
  
  dum2 <- data.frame(
    group = as.character(rep(group_names, each = nrow(pred_dat))),
    link_abund_est = comp_abund_pred[ , "Estimate"],
    link_abund_se =  comp_abund_pred[ , "Std. Error"]
  ) %>% 
    cbind(dum, .) %>%
    mutate(
      comp_abund_est = exp(.data$link_abund_est),
      comp_abund_low = exp(.data$link_abund_est + 
                             (stats::qnorm(0.025) * .data$link_abund_se)),
      comp_abund_up = exp(.data$link_abund_est + 
                            (stats::qnorm(0.975) * .data$link_abund_se))
    ) 
  
  ggplot(dum2, aes_string(x = x_var)) +
    geom_line(aes_string(y = dum2$comp_abund_est, colour = colour_var)) +
    geom_ribbon(aes_string(ymin = dum2$comp_abund_low, ymax = dum2$comp_abund_up, 
                           fill = colour_var), 
                alpha = 0.5) +
    scale_fill_brewer(name = colour_var, palette = "Dark2") +
    scale_colour_brewer(name = colour_var, palette = "Dark2") +
    labs(y = "Predicted Abundance") +
    facet_wrap(~group, scales = "free_y") +
    theme_classic() +
    theme(legend.position = "top",
          axis.text = element_text(size=9),
          plot.margin = unit(c(5.5, 10.5, 5.5, 5.5), "points")) +
    coord_cartesian(expand = FALSE, ylim = c(0, NA))
}
