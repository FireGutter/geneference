#'
#' @title covmatrix
#' @description Design the theoretical covariance matrix for the multivariate
#' normal used to model the distribution of liabilities.
#' The first entry in the vector of liabilities is the genetic liability, the
#' next is the full liability and the two following are the liabilities of the
#' parents. Liabilities of siblings constitute the remaining entries in the
#' vector.
#'
#'
#' @param hsq squared heritability parameter.
#' @param sib number of siblings.
#'
#' @return A covariance matrix for an individual with sib number of siblings,
#' such that the number of rows and columns are equal to 4 + number of
#' siblings.
#'
#'
#' @examples
#' covmatrix(2, 0.5)
#'
#' @export

covmatrix <- function(sib = 0, hsq) {
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
