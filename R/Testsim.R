#'
#' @title testsim
#' @description A quick way of simulating a small amount of data in the R-session. Note that this function
#' only simulate data for individuals. So the simulation should not be used in a study with siblings/parents.
#'
#' @param n is the amount of enteties (humans)
#' @param m is the amount of SNPS per entety
#' @param q is the amount of causal SNPs
#' @param hsq heritability parameter that is squard
#' @param k is the prevalence of trait parameter
#' @param to_ped is a T/F parameter indicatating if the output should be in a plink-friendly format or in binary
#'
#' @return This function returns a list containing a matrix of genotypes, causal SNPs, the enteties liabilities and the phenotypes
#'
#' @importFrom stats qnorm rnorm runif rbinom
#'
#' @export
#' @examples
#' testsim(10, 15, 5, 0.5, 0.05, TRUE)
#'
#'

testsim <- function(n, m, q, hsq, k, to_ped = T){

  n <- as.numeric(n)
  m <- as.numeric(m)
  q <- as.numeric(q)
  hsq <- as.numeric(hsq)
  k <- as.numeric(k)

  h <- sqrt(hsq)

  if(n != as.integer(n) | n <= 0){
    stop("The input of n must be an positive integer!")
  }
  else if(m != as.integer(m) | m <= 0){
    stop("The input of m must be an positive integer!")
  }
  else if(q != as.integer(q) | q < 1 | q > m){
    stop("The input of q must be an positive integer between 1 and m!")
  }

  if(h < 0 | h > 1){
    stop("The heritability parameter 'h' must be between 0 and 1")
  }
  else if(k < 0 | k > 1){
    stop("the prevalence of trait 'k' must be between 0 and 1 (it is a precentage)")
  }


  MAFs <- runif(m, 0.01, 0.49) # Calculate the maf (Minor Allele Frequency). Alle individuals have the same MAFs.

  causual_SNP <- sample.int(m, size = q, replace = F) # Random selection of the SNPs on the genome.

  beta <- matrix(0, nrow = m, ncol = 1)
  beta[causual_SNP] <- rnorm(q, 0, h/sqrt(q)) # Create the effect-size for the causal SNPs

  persons <- t(sapply(1:n, function(y) rbinom(m, 2, MAFs))) # Randomly sample the postions

  mu <- 2*MAFs # Find the mean value per SNP
  sigma <- sqrt(2*MAFs*(1-MAFs)) # Find the variance

  lg <- sweep(sweep(persons, 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta # Find the genetic liabilities

  liability <- lg + rnorm(n, 0, sqrt(1 - h^2)) # Find the full liabilities

  critical <- qnorm(1 - k) # Find the threshold of disease liability

  y <- sapply(liability, function(x) ifelse(x > critical, 2, 1)) # Turn the matrix into binary on the format 1 or 2.

  # id <- matrix(1:n)

  if(to_ped){
    return(list("genotypes" = to_ped(persons, 0), "SNP-Turn" = beta, "liability" = liability, "Phenotype" = y))
  }
  else{
    return(list("genotypes" = persons, "SNP-Turn" = beta, "liability" = liability, "Phenotype" = y))
  }
}
