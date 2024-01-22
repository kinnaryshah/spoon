# run examples from function documentation
example(weighted_nnSVG, echo = FALSE)

test_that("example object has correct class", {
  expect_s4_class(spe_results, "SpatialExperiment")
})

test_that("example objects have correct dimensions", {
  expect_equal(dim(spe_results), c(6, 3620))
  expect_equal(dim(results), c(6, 6))
})

