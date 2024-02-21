testthat::test_that("Test Truth 1", {
  result <- calc_roc(c(1,1,1,0), c(1,1,1,0))
  expect_equal(result@predictions[[1]][1], 1)
})

testthat::test_that("Test Truth 2", {
  result <- calc_roc(c(1,1,1,0), c(1,1,1,0))
  expect_equal(result@predictions[[1]][4], 0)
})
