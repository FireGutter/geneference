#'
#' @title Multivariate normal distributions
#' @description Simulate random values from a multivariate normal distribution.
#'
#' @param S a valid covariance matrix.
#' @param n number of random vectors that is to be simulated.
#'
#' @import stats
#'
#' @return Returns a matrix with n rows and i columns, where i is determined by
#' the size of the specified covariance matrix.
#'
#' @export
#' @examples
#' rmnorm(n = 10, S = covmatrix(sib = 0, hsq = 0.5))

rmnorm <- function(n, S) {
  stopifnot("n needs to be an integer greater than 0" =
              (n > 0 && class(n) == "numeric" && n == round(n)),
            "S needs to be a valid covariance matrix" =
              (nrow(S) == ncol(S) && all(is.numeric(S))) &&
                  all(S == t(S)) && all(S >= 0))

  n_liab <- nrow(S)
  C <- chol(S)
  emptytest <- matrix(n, nrow = n_liab)
  Z <- sapply(emptytest, rnorm)
  X_t <- Z %*% C

  return(X_t)
}
