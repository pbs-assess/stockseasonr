context("gen_tmb inputs")

test_that("gen_tmb works", {
  x <- gen_tmb(catch_ex, comp_ex)
  expect_true("pars" %in% names(x))
  expect_equal(class(x), "list")
})
