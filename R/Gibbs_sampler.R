#' Conditional means and covariances of MVN
#'
#' Helper function for our Gibbs sampler.
#'
#' @param S covariance matrix of a multivariate normal distribution
#' @param a numeric vector of liabilities
#' @param person integer indicating which liability (index) we want to
#' calculate the distribution for, conditional on the other liabilities
#'
#' @noRd
condition <- function(S, a, person) {
  s11 <- S[person, person]
  s12 <- matrix(S[person, -c(person)], nrow = 1)
  s21 <- matrix(S[-c(person), person], ncol = 1)
  s22 <- S[-c(person), -c(person)]
  s22_inv <- solve(s22)
  
  cond_mu <- s12 %*% s22_inv %*% a[-person]
  cond_S <- s11 - s12 %*% s22_inv %*% s21
  
  return(c(cond_mu, cond_S))
}

#' Gibbs sampler for truncated normal distribution
#'
#' @param conf configuration to calculate posterior mean genetic liability for
#' @param burn_in number of iterations before burn in of the sampler is
#' considered to be completed
#' @param start_value initial value for the Gibbs sampler
#' @param alpha significance level
#'
#' @import stats
#' 
#' @noRd
generate_pmgl <- function(conf, burn_in, start_value, alpha) {
  conf_list <- strsplit(conf, "")
  n <- length(conf_list[[1]]) + 1
  
  liabs <-  matrix(rep(start_value, n), nrow = n, ncol = 1)
  genliabs <- matrix(NA, nrow = 10000, ncol = 1)
  
  S <- covmatrix(n - 4, 0.5)
  
  crit <- qnorm(1 - alpha)
  
  i <- 0
  while (TRUE) {
    # Batch sampling with sample size being burn_in size
    for (s in 1:burn_in) {
      i <- i + 1
      # Updating each of the liabilities
      for (j in 1:n) {
        # Finding conditional mean and variance
        param <- condition(S, matrix(liabs, nrow = n, ncol = 1), j)
        mu <- param[1]
        var <- param[2]
        # If genetic liability, then sample from normal distribution
        if (j == 1) {
          liabs[j] <- rnorm(1, mu, sqrt(var))
        } else{ # Else make transformation to sample from truncated distribution
          crit_U <- pnorm(crit, mu, sqrt(var))
          # If case, then sample from critical value and above
          if (conf_list[[1]][j - 1] == "2") {
            U <- runif(1, crit_U, 1)
          } else{ # Else, sample from critical value and below
            U <- runif(1, 0, crit_U)
          }
          liabs[j] <- qnorm(U, mu, sqrt(var))
        }
      }
      # If burn in is over, then save genetic liability and check SEM
      if (i > burn_in) {
        genliabs[i - burn_in] <- liabs[1]
        # If SEM is less than 0.01, then return mean of sampled genetic data
        if (sd(genliabs[1:(i - burn_in), ]) / sqrt(i - burn_in) < 0.01) {
          return(mean(genliabs[1:(i - burn_in)]))
        }
      }
    }
  }
}
