test_that("Checking same output for seed 42", {
  set.seed(42)
  test <- testsim(10, 15, 5, 0.5, 0.05, TRUE)
  expect_equal(test, data_for_testsim)
})

test_that("Output is a list of 4", {
  set.seed(42)
  test <- testsim(10, 15, 5, 0.5, 0.05, TRUE)
  expect_equal(class(test), "list")
  expect_equal(length(test), 4)
})
