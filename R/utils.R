## Utility functions for generating TMB inputs

map_foo <- function(x, tmb_pars) {
  out_list <- vector(mode = "list", length = length(x))
  names(out_list) <- x
  for (i in seq_along(x)) {
    out_list[[i]] <- as.factor(rep(NA, length(tmb_pars[[x[i]]])))
  }
  return(out_list)
}

# helper function to generate random walk index when multiple random intercept
# variables are present; re_indexes should be matrix 
# (i.e. sdmTMB_dummy$tmb_data$RE_indexes)
rw_index_foo <- function(re_indexes) {
  rw_index <- NULL
  for (g in 1:ncol(re_indexes)) {
    rw_index <- c(rw_index, sort(unique(re_indexes[ , g])))
  }
  return(rw_index)
}
