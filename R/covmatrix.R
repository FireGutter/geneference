#' @title Create covariance matrix for liabilities
#'
#' @description Design the theoretical covariance matrix from the multivariate
#' normal distribution used to model the liabilities.
#'
#' @details
#' The covariance matrix is explained in detail in
#' \code{vignette("liability-distribution")}.\cr
#'
#' @param hsq heritability parameter.
#' @param sib number of siblings.
#'
#' @return A covariance matrix for the liabilities of a family with \code{sib}
#' number of siblings. 
#'
#' @examples
#' covmatrix(0.5, 2)
#'
#' @export

covmatrix <- function(hsq, sib = 0) {
  stopifnot("sib needs to be a non-negative integer" =
              (sib >= 0 && is.numeric(sib) && round(sib) == sib &&
                 length(sib) == 1),
            "hsq needs to be a number between 0 and 1" =
              (hsq > 0 && hsq < 1 && is.numeric(hsq) && length(hsq) == 1))

  s <- matrix(c(hsq, hsq, rep(0.5 * hsq, sib + 2),
                hsq, 1, rep(0.5 * hsq, sib + 2),
                0.5 * hsq, 0.5 * hsq, 1, 0, rep(0.5 * hsq, sib),
                0.5 * hsq, 0.5 * hsq, 0, 1, rep(0.5 * hsq, sib)),
              nrow = 4, ncol = 4 + sib, byrow = T)

  if (sib != 0) {
    for (i in 1:sib) {
      s <- rbind(s, c(rep(0.5 * hsq, 3 + i), 1, rep(0.5 * hsq, sib - i)))
    }
  }

  return(s)
}
