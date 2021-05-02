#'
#' @title covmatrix
#' @description This function create a covariance matrix for a entety regarding the heritability.
#' The first row is full liability for a entety, next is genetic liability.
#' Row 3 and 4 is the entetys parents liability.
#'
#'
#' @param hsq heritability parameter that is squard
#' @param sib is the amount siblings to the entety
#'
#' @return This function returns the covariance matrix for an entety with sib - siblings.
#'
#' @examples
#' covmatrix(2, 0.5)
#'
#' @export

covmatrix <- function(sib = 0, hsq){
  s <- matrix(c(hsq, hsq, rep(0.5*hsq, sib + 2), hsq, 1, rep(0.5*hsq, sib + 2), 0.5*hsq, 0.5*hsq, 1, 0, rep(0.5*hsq, sib), 0.5*hsq, 0.5*hsq, 0, 1, rep(0.5*hsq, sib)), nrow = 4, ncol = 4 + sib, byrow = T)

  if(sib != 0){
    for(i in 1:sib){
      s <- rbind(s, c(rep(0.5*hsq, 3 + i), 1, rep(0.5*hsq, sib - i)))
    }
  }

  return(s)
}
