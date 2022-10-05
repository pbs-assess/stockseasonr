context("fit_stockseasonr test")

test_that("fit_stockseasonr works", {
  x <- fit_stockseasonr(comp_formula = agg ~ 1 + region + 
                          s(month_n, bs = "tp", k = 4, m = 2) + 
                          (1 | year),
                        comp_dat = comp_ex,
                        model = "dirichlet",
                        random_walk = TRUE,
                        fit = TRUE,
                        nlminb_loops = 2, newton_loops = 1)
  expect_true("sdr" %in% names(x))
  expect_equal(class(x), "list")
})
