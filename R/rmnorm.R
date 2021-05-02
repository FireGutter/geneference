#'
#' @title rmnorm
#' @description Simulate random values form a multi-variate normal distribution.
#'
#'
#' @param S is a valid covariance matrix.
#' @param n is the amount of random vectors that is to be simulated.
#'
#' @return This function returns a matrix with n rows and i columns according to the amount of columns in S.
#'
#' @export

rmnorm <- function(n, S){
  #output is of shape  individuals x liabilities

  #assumes that the given S is a valid covariance matrix
  n_liab <- nrow(S)
  C <- chol(S)
  emptytest <- matrix(n, nrow = n_liab)
  Z <- sapply(emptytest, rnorm)
  X_t <- Z %*% C

  return(X_t)
}
