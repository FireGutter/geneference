test_that("Standard matrix", {
  expect_equal(covmatrix(sib = 0, hsq = 0.5), 
               matrix(c(0.5, 0.5, 0.25, 0.25, 
                        0.5, 1, 0.25, 0.25, 
                        0.25, 0.25, 1, 0, 
                        0.25, 0.25, 0, 1), nrow = 4))
})

test_that("Check input", {
  expect_error(covmatrix(sib = -1, hsq = 0.5), "sib needs to be a non-negative integer")
  expect_error(covmatrix(sib = 0, hsq = -1), "hsq needs to be a number between 0 and 1")
})

test_that("Output", {
  expect_type(covmatrix(sib = 0, hsq = 0.5), "double")
  expect_equal(dim(covmatrix(sib = 0, hsq = 0.5)), c(4, 4))
})
