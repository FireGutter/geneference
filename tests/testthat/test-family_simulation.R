test_that("Check input", {
  expect_error(family_simulation(n = -1, m = 10000, q = 100, hsq = 0.5, k = 0.05, sib = 3, path = ""), "n needs to be an integer greater than 0")
  expect_error(family_simulation(n = 10000, m = -1, q = 100, hsq = 0.5, k = 0.05, sib = 3, path = ""), "m needs to be an integer greater than 0")
  expect_error(family_simulation(n = 10000, m = 10000, q = 100000, hsq = 0.5, k = 0.05, sib = 3, path = ""), "q needs to be an integer greater than 0 and smaller than m")
  expect_error(family_simulation(n = 10000, m = 10000, q = 100, hsq = 100, k = 0.05, sib = 3, path = ""), "hsq needs to be a number between 0 and 1")
  expect_error(family_simulation(n = 10000, m = 10000, q = 100, hsq = 0.5, k = -10, sib = 3, path = ""), "k needs to be a number between 0 and 1")
  expect_error(family_simulation(n = 10000, m = 10000, q = 100, hsq = 0.5, k = 0.05, sib = -3, path = ""), "sib needs to be a non-negative integer")
})

