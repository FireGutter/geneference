#' Conditional means and covariances of MVN
#'
#' @param S covariance matrix of a multivariate normal distribution
#' @param a numeric vector of liabilities
#' @param person integer indicating which liability (index) we want to
#' calculate the distribution for, conditional on the other liabilities
#' 
#' @noRd
condition <- function(S, a, person){
  s11 <- S[person,person]
  s12 <- matrix(S[person, -c(person)], nrow = 1)
  s21 <- matrix(S[-c(person),person], ncol = 1)
  s22 <- S[-c(person),-c(person)]
  s22_inv <- solve(s22)
  
  cond_mu <- s12 %*% s22_inv %*% a[-person]
  cond_S <- s11 - s12 %*% s22_inv %*% s21
  
  return(c(cond_mu, cond_S))
}