#'
#' @title Simulation directly in the R session without family history
#'
#' @description Simulate genetic data in R, including genotypes,
#' phenotype status and liabilities, for individuals.
#' 
#' @details
#' A quick way of simulating a small amount of data to be stored
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
#' @param to_ped if `TRUE`, write the genotype output in PLINK's .ped format.
#' Else, write the number of minor alleles for each genotype at each SNP.
#'
#' @return Returns a list with 5 entries, containing a matrix of genotypes, a
#' vector specifying effect sizes of SNPs, a vector with liability of
#' individual i, a vector specifying phenotype  status of individual i, and the
#' transformed phenotype status used for linear regression in PLINK.
#'
#' @section Warning:
#' Since this function saves the generated data directly in R (that is, on RAM),
#' the function should not be used for large simulations, but rather to get an
#' idea of the general framework.\cr
#' The function does not run in parallel.
#'
#' @import stats
#'
#' @export
#' @examples
#' sim_test(10, 15, 5, 0.5, 0.05, TRUE)
sim_test <- function(n, m, q, hsq, k, to_ped = T) {
  stopifnot("n needs to be an integer greater than 0" =
              (n > 0 && is.numeric(n) && n == round(n) && length(n) == 1),
            "m needs to be an integer greater than 0" =
              (m > 0 && is.numeric(m) && m == round(m) && length(m) == 1),
            "q needs to be an integer greater than 0 and smaller than m" =
              (q > 0 && is.numeric(q) && q == round(q) && length(q) == 1 
               && q <= m),
            "hsq needs to be a number between 0 and 1" =
              (hsq > 0 && hsq < 1 && is.numeric(hsq) && length(hsq) == 1),
            "k needs to be a number between 0 and 1" =
              (k > 0 && k < 1 && is.numeric(k) && length(k) == 1),
            "to_ped needs to be logical" = class(to_ped) == "logical")

  # Calculate the Minor Allele Frequency. All individuals have the same MAFs.
  MAFs <- runif(m, 0.01, 0.49)

  # Random selection of the causal SNPs on the genome.
  causual_SNP <- sample.int(m, size = q, replace = F)

  # Create the effect-sizes for the causal SNPs
  beta <- matrix(0, nrow = m, ncol = 1)
  beta[causual_SNP] <- rnorm(q, 0, sqrt(hsq / q))

  # Determine number of risk-allelles for each genotype
  persons <- t(vapply(1:n, function(y) rbinom(m, 2, MAFs), 
                      FUN.VALUE = numeric(m)))

  mu <- 2 * MAFs # Find the mean value per SNP
  sigma <- sqrt(2 * MAFs * (1 - MAFs)) # Find the variance

  # Calculate the genetic liabilities
  lg <- sweep(sweep(persons, 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta

  liability <- lg + rnorm(n, 0, sqrt(1 - hsq)) # Find the full liabilities

  critical <- qnorm(1 - k) # Find the threshold of disease liability

  # Turn the matrix into binary on the format 1 or 2.
  pheno <- vapply(liability, function(x) ifelse(x > critical, 2, 1),
                  FUN.VALUE = matrix(n))
  line_pheno <- pheno + 1

  if (to_ped) {
    return(list("genotypes" = to_ped(persons, 0),
                "effect_sizes" = beta,
                "liability" = liability,
                "Phenotype" = pheno,
                "PLINK_linear_phenotype" = line_pheno))
  }
  else{
    return(list("genotypes" = persons,
                "effect_sizes" = beta,
                "liability" = liability,
                "Phenotype" = pheno,
                "PLINK_linear_phenotype" = line_pheno))
  }
}
