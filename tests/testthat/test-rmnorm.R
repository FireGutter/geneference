test_that("Check input", {
  expect_error(rmnorm(n = -1, S = covmatrix(sib = 0, hsq = 0.5)), 
               "n needs to be an integer greater than 0")
  expect_error(rmnorm(n = 10, S = matrix(c(1, 1))), 
               "S needs to be a valid covariance matrix")
})

test_that("Correct output for seed = 42", {
  set.seed(42)
  expect_equal(rmnorm(n = 10, S = covmatrix(sib = 0, hsq = 0.5)), 
               matrix(c(0.96941401,  1.8920962,  0.19787287,  0.94734819,
                        -0.39930191,  1.2176006, -1.86591242,  0.69093897,
                        0.25677056, -0.7253023, -0.03242868,  1.10967835,
                        0.44750144,  0.2503680,  1.35997486, -0.50232329,
                        0.28586087,  0.1915885,  1.91572159,  0.35717216,
                        -0.07504136,  0.3746435, -0.44018768, -1.56963797,
                        1.06880745,  0.8678103,  0.29375026, -0.15748512,
                        -0.06693405, -1.9453317, -1.68275507, -0.58564180,
                        1.42724110, -0.2984296,  1.14400221, -1.58298451,
                        -0.04434556,  0.8891155, -0.62083317,  0.09679317), 
                        nrow = 10, byrow = T))
})


test_that("Output", {
  expect_type(rmnorm(n = 10, S = covmatrix(sib = 0, hsq = 0.5)), "double")
  expect_equal(dim(rmnorm(n = 10, S = covmatrix(sib = 0, hsq = 0.5))), c(10, 4))
})


