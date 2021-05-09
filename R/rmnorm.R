#'
#' @title Multivariate normal distributions
#' @description Simulate random values from a multivariate normal distribution.
#'
#'
#' @param S a valid covariance matrix.
#' @param n number of random vectors that is to be simulated.
#'
#' @return Returns a matrix with n rows and i columns, where i is determined by
#' the size of the specified covariance matrix.
#'
#' @export

rmnorm <- function(n, S) {

  n_liab <- nrow(S)
  C <- chol(S)
  emptytest <- matrix(n, nrow = n_liab)
  Z <- sapply(emptytest, rnorm)
  X_t <- Z %*% C

  return(X_t)
}
