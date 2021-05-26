#'
#' @title Simulation with random family sizes
#' @md
#' @description Simulate genotypes for individuals with family history, where
#' each individual is given a random number of siblings. Parents' genotypes are
#' simulated and used for simulating the genotypes of individuals and their
#' siblings. This function call the family_dist_simulaiton with the randomized
#' input.
#' \code{family_dist_simulaiton} makes use of parallel computation in order to
#' decrease the running time.
#'
#' @section Warning:
#' Simulating large datasets lead to very large files (eg. in the order of 100k
#' genotypes with 100k SNPs takes up approximately 40GB of space).
#' Please ensure that you have sufficient disk space at the disk where
#' \code{path} resides prior to running the simulation.
#' See FUNCTION-NAME DOCUMENTATION REFERENCE(this is a placeholder until
#' function is implemented) in order to convert the generated .ped to the
#' smaller .bed format.
#'
#' @param n number of genotypes (individuals).
#' @param m number of SNPS per genotype.
#' @param q number of causal SNPs, i.e. SNPs that effect chances of having
#' the phenotype.
#' @param hsq squared heritability parameter.
#' @param k prevalence of phenotype.
#' @param sib_fert either the distribution-vector of siblings-fertility rate.
#' See meaning below.
#' @param dist if sib_fert is a distribution-vector then dist can be used to
#' specify how many siblings the different individuals get.
#' @param path directory where the files will be stored. If nothing is
#' specified, \code{family_simulation} writes its files in the current
#' working directory.
#' @details Note that there are two types of input to sib_fert:
#' distribution-vector or siblings-fertility. If a distribution-vector is used
#' a dist vector must also be specified. E.g. sib_fert = c(1/6, 1/6, 2/3)
#' and dist = c(0, 1, 2). In this case the function is asked to randomly choose
#' 1/6 individuals to have 0 siblings and 2/3 to have 2 siblings. The other case
#' a value is given to sib_fert. E.g. sib_fert = 2. Here the function use the
#' poisson distribution to randomly select how many siblings the different
#' individuals has. Note that in this case, no dist has to be specified, since
#' the function overwrites it with the range of the simulated siblings.
#'
#' @return a list where first entry is the distribution of the simulated individuals
#' and the second entry is the number of siblings (dist) for these individuals.
#' the function also call family_dist_simulaiton that prints the following five files to
#' the \code{path} parameter specified in the function call:
#' * Three text files:
#'     * beta.txt - a file of \code{m} rows with one column. The i'th row is
#'     the true effect of the i'th SNP.
#'     * MAFs.txt - a file of \code{m} rows with one column. The i'th row is
#'     the true Minor Allelle Frequence of the i'th SNP.
#'     * phenotypes.txt - a file of \code{n} rows, number of columns depend on
#'     number of siblings. Note that some individuals has fewer siblings as
#'     others. Thus we fill the empty value-slots with "-9". The file
#'     contains the phenotype and liability of
#'     each individual as well as information on the liabilities and phenotype
#'     status of their parents and siblings.
#' * genotypes.map - a file created such that PLINK will work with the genotype
#' data.
#' * genotypes.ped - the simulated genotypes in a PLINK-readable format.
#'
#' @import stats
#'
#' @export

sim_random_family <- function(n, m, q, hsq, k, sib_fert, dist = 0, path = "") {


  if (length(sib_fert) > 1) {
    individual_distribution <- c(rmultinom(1, n, sib_fert))
  }
  else {
    vals <- rpois(n, sib_fert)
    ran <- range(vals)
    individual_distribution <- as.vector(table(vals))
    dist <- ran[1]:ran[2]
  }
  family_dist_simulaiton(n = individual_distribution, m = m, q = q, hsq = hsq, k = k, dist = dist, path = path)
  return(list("simulated_individuals(n)" = individual_distribution, "number_of_siblings(dist)" = dist))
}

