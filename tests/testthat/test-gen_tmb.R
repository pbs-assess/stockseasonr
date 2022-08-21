context("gen_tmb inputs")

test_that("gen_tmb works", {
  x <- fit_stockseasonr(abund_formula = catch ~ 1 +
                          s(month_n, bs = "tp", k = 3, m = 2) +
                          region ,
                        abund_dat = catch_ex,
                        comp_formula = agg ~ 1 + region +
                          s(month_n, bs = "tp", k = 4, m = 2),
                        comp_dat = comp_ex,
                        model = "integrated",
                        fit = FALSE)
  expect_true("pars" %in% names(x))
  expect_equal(class(x), "list")
})
