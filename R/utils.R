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


# functions from sdmTMB to replace old glmmTMB::splitForm 
safe_deparse <- function(x, collapse = " ") {
  paste(deparse(x, 500L), collapse = collapse)
}

barnames <- function (bars) {
  vapply(bars, function(x) safe_deparse(x[[3]]), "")
}

split_form <- function (f) {
  b <- lme4::findbars(f)
  bn <- barnames(b)
  fe_form <- lme4::nobars(f)
  list(bars = b, barnames = bn, form_no_bars = fe_form, n_bars = length(bn))
}


## helpful for parsing smoothers (from sdmTMB but currently unused)
# from brms:::rm_wsp()
# rm_wsp <- function (x) {
#   out <- gsub("[ \t\r\n]+", "", x, perl = TRUE)
#   dim(out) <- dim(x)
#   out
# }
# # from brms:::all_terms()
# all_terms <- function (x) {
#   if (!length(x)) {
#     return(character(0))
#   }
#   if (!inherits(x, "terms")) {
#     x <- terms(stats::as.formula(x))
#   }
#   rm_wsp(attr(x, "term.labels"))
# }
# 
# get_smooth_terms <- function(terms) {
#   x1 <- grep("s\\(", terms)
#   x2 <- grep("t2\\(", terms)
#   c(x1, x2)
# }
# 
# parse_smoothers <- function(formula, data, knots = NULL, newdata = NULL, basis_prev = NULL) {
#   terms <- all_terms(formula)
#   smooth_i <- get_smooth_terms(terms)
#   basis <- list()
#   basis_out <- list()
#   Zs <- list()
#   Xs <- list()
#   labels <- list()
#   classes <- list()
#   if (length(smooth_i) > 0) {
#     has_smooths <- TRUE
#     smterms <- terms[smooth_i]
#     ns <- 0
#     ns_Xf <- 0
#     for (i in seq_along(smterms)) {
#       if (grepl('bs\\=\\"re', smterms[i])) cli_abort("bs = 're' is not currently supported for smooths")
#       if (grepl('fx\\=T', smterms[i])) cli_abort("fx = TRUE is not currently supported for smooths")
#       # if (grepl('m\\=3', smterms[i])) cli_abort("m > 2 is not currently supported for smooths")
#       obj <- eval(str2expression(smterms[i]))
#       labels[[i]] <- obj$label
#       classes[[i]] <- attr(obj, "class")
#       if (is.null(newdata)) {
#         basis[[i]] <- mgcv::smoothCon(
#           object = obj, data = data,
#           knots = knots, absorb.cons = TRUE,
#           diagonal.penalty = FALSE
#         )
#         basis_out[[i]] <- mgcv::smoothCon( # to be used on prediction
#           object = obj, data = data,
#           knots = knots, absorb.cons = TRUE,
#           diagonal.penalty = FALSE#  modCon = 3 # modCon set differently as per brms
#         )
#       } else {
#         basis[[i]] <- basis_prev[[i]] # predicting on new data
#       }
#       for (j in seq_along(basis[[i]])) { # elements > 1 with `by` terms
#         ns_Xf <- ns_Xf + 1
#         rasm <- mgcv::smooth2random(basis[[i]][[j]], names(data), type = 2)
#         if (!is.null(newdata)) {
#           rasm <- s2rPred(basis[[i]][[j]], rasm, newdata)
#         }
#         for (k in seq_along(rasm$rand)) { # elements > 1 with if s(x, y) or t2()
#           ns <- ns + 1
#           Zs[[ns]] <- rasm$rand[[k]]
#         }
#         Xs[[ns_Xf]] <- rasm$Xf
#       }
#     }
#     sm_dims <- unlist(lapply(Zs, ncol))
#     Xs <- do.call(cbind, Xs) # combine 'em all into one design matrix
#     b_smooth_start <- c(0, cumsum(sm_dims)[-length(sm_dims)])
#   } else {
#     has_smooths <- FALSE
#     sm_dims <- 0L
#     b_smooth_start <- 0L
#     Xs <- matrix(nrow = 0L, ncol = 0L)
#   }
#   list(Xs = Xs, Zs = Zs, has_smooths = has_smooths, labels = labels,
#        classes = classes, basis_out = basis_out,
#        sm_dims = sm_dims, b_smooth_start = b_smooth_start)
# }
