#'
#' @title testsim
#' @description A quick way of simulating a small amount of data to be stored
#' in the R-session. Note that this function only simulates data for
#' individuals without family history, i.e. they are not simulated as
#' offspring. Hence, the simulation cannot be used in a study of family
#' history.
#'
#' @param n number of genotypes (individuals).
#' @param m number of SNPS per genotype.
#' @param q number of causal SNPs, i.e. SNPs that effect chances of having
#' the phenotype.
#' @param hsq squared heritability parameter.
#' @param k prevalence of phenotype.
#' @param to_ped TRUE/FALSE indicating if the output should be in a
#' PLINK-friendly format or not.
#'
#' @return Returns a list of 4 entries, containing a matrix of genotypes, a
#' vector specifying indexes of causal SNPs, a vector the liability of
#' individual i, and a vector specifying whether individual i has the
#' phenotype.
#'
#' @importFrom stats qnorm rnorm runif rbinom
#'
#' @export
#' @examples
#' testsim(10, 15, 5, 0.5, 0.05, TRUE)
#'
#'

testsim <- function(n, m, q, hsq, k, to_ped = T) {

  n <- as.numeric(n)
  m <- as.numeric(m)
  q <- as.numeric(q)
  hsq <- as.numeric(hsq)
  k <- as.numeric(k)

  h <- sqrt(hsq)

  if (n != as.integer(n) | n <= 0) {
    stop("The input of n must be a positive integer!")
  }
  else if (m != as.integer(m) | m <= 0) {
    stop("The input of m must be a positive integer!")
  }
  else if (q != as.integer(q) | q < 1 | q > m) {
    stop("The input of q must be a positive integer between 1 and m!")
  }

  if (hsq < 0 | hsq > 1) {
    stop("The heritability parameter 'hsq' must be between 0 and 1")
  }
  else if (k < 0 | k > 1) {
    stop("the prevalence of trait 'k' must be between 0 and 1")
  }

  # Calculate the Minor Allele Frequency. All individuals have the same MAFs.
  MAFs <- runif(m, 0.01, 0.49)

  # Random selection of the causal SNPs on the genome.
  causual_SNP <- sample.int(m, size = q, replace = F)

  # Create the effect-sizes for the causal SNPs
  beta <- matrix(0, nrow = m, ncol = 1)
  beta[causual_SNP] <- rnorm(q, 0, h / sqrt(q))

  # Determine number of risk-allelles for each genotype
  persons <- t(sapply(1:n, function(y) rbinom(m, 2, MAFs)))

  mu <- 2 * MAFs # Find the mean value per SNP
  sigma <- sqrt(2 * MAFs * (1 - MAFs)) # Find the variance

  # Calculate the genetic liabilities
  lg <- sweep(sweep(persons, 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta

  liability <- lg + rnorm(n, 0, sqrt(1 - h^2)) # Find the full liabilities

  critical <- qnorm(1 - k) # Find the threshold of disease liability

  # Turn the matrix into binary on the format 1 or 2.
  y <- sapply(liability, function(x) ifelse(x > critical, 2, 1))

  if (to_ped) {
    return(list("genotypes" = to_ped(persons, 0),
                "SNP-Turn" = beta,
                "liability" = liability,
                "Phenotype" = y))
  }
  else{
    return(list("genotypes" = persons,
                "SNP-Turn" = beta,
                "liability" = liability,
                "Phenotype" = y))
  }
}
