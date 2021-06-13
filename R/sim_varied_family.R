#'
#' @title Simulation with varying family sizes
#' @md
#' @description Simulate genetic data, including genotypes,
#' phenotype status and liabilities, for individuals and their family,
#' where each individual has a specified number of siblings.
#'
#' @details
#' Parents' genotypes are simulated and used for creating the genotypes of
#' the individuals and their siblings. For the methodology behind the
#' simulation, see `vignette("liability-distribution")`.\cr
#' Note: Each entry in `dist` denotes a number of siblings. Each entry in
#' `n` then denotes how many individuals have the corresponding number of
#' siblings.\cr
#' E.g., `n = c(100, 200, 300, 400)` and `dist = c(0, 2, 3, 5)` would
#' give a total of 100 + 200 + 300 + 400 = 1000 individuals, where 100
#' individuals have 0 siblings, 200 have 2 siblings, and so on.\cr
#' Since individuals have a different number of siblings, some entries in
#' `phenotypes.txt` will be missing, denoted by -9.\cr
#' E.g., an individual with 1 sibling in a dataset where the maximum number of
#' siblings is 3, would have -9 in all columns relating to sibling 2 and
#' sibling 3.\cr
#' \code{sim_varied_family} makes use of parallel computation in order to
#' decrease the running time. As one CPU core is left unused, the user
#' should be able to do other work while the simulation is running.
#'
#' @param n number of genotypes (individuals), given as a vector of same length
#' as `dist`.
#' @param m number of SNPs per genotype.
#' @param q number of causal SNPs, i.e. SNPs that effect chances of having
#' the phenotype.
#' @param hsq squared heritability parameter.
#' @param k prevalence of phenotype.
#' @param dist the distribution of siblings. Given as a vector with the same
#' length as `n`.
#' @param path directory where the files will be stored. If nothing is
#' specified, \code{sim_varied_family} writes its files in the current
#' working directory.
#'
#' @return Does not return any value, but prints the following five files to
#' the \code{path} parameter specified in the function call:
#' * Three text files:
#'     * beta.txt - a file of \code{m} rows with one column. The i'th row is
#'     the true effect of the i'th SNP.
#'     * MAFs.txt - a file of \code{m} rows with one column. The i'th row is
#'     the true Minor Allelle Frequency of the i'th SNP.
#'     * phenotypes.txt - a file of \code{n} rows, number of columns depend on
#'     number of siblings. The file contains the phenotype status and liability
#'     of each individual as well as information on the liabilities and
#'     phenotype status of their parents and siblings.
#' * genotypes.map - a file created such that PLINK will work with the genotype
#' data.
#' * genotypes.ped - the simulated genotypes in a PLINK-readable format.
#' Note: The function only saves genotype data for the target individual.
#'
#' @section Warning:
#' Simulating large datasets takes time and generates large files. For details
#' on time complexity and required disk space, see
#' `vignette("sim-benchmarks")`.\cr
#' The largest file generated is `genotypes.ped`. See `convert_geno_file()` to convert it
#' to another file format, thereby reducing its size significantly.
#'
#' @importFrom data.table :=
#' @import stats
#'
#' @export

sim_varied_family <- function(n, m, q, hsq, k, dist, path = ""){
  stopifnot("n and dist needs to have same length" = length(n) == length(dist),
            "n needs to be a vector of positive integers" =
              (all(n > 0) && is.numeric(n) && all(n == round(n))),
            "m needs to be an integer greater than 0" =
              (m > 0 && is.numeric(m) && m == round(m)),
            "q needs to be an integer greater than 0 and smaller than m" =
              (q > 0 && is.numeric(q) && q == round(q) && length(q) == 1
               && q <= m),
            "hsq needs to be a number between 0 and 1" =
              (hsq > 0 && hsq < 1 && is.numeric(hsq) && length(hsq) == 1),
            "k needs to be a number between 0 and 1" =
              (k > 0 && k < 1 && is.numeric(k) && length(k) == 1),
            "dist needs to be a vector of non-negative integers" =
              (all(dist >= 0) && is.numeric(dist) && all(dist == round(dist))),
            "path needs to be default or a valid path ending with '/' or '\\\\'"
            = (path == "" || (dir.exists(path))
               && (substr(path, nchar(path), nchar(path)) == "/" ||
                     substr(path, nchar(path), nchar(path)) == "\\")))


  path <- path_validation(path)

  # Set worker nodes:
  future::plan(future::multiprocess, workers = max(future::availableCores(logical = F) - 1, 1))

  parent_maker <- function(m, number, MAFs) {
    vapply(1:number, function(y) {rbinom(m, 2, MAFs)}, FUN.VALUE = numeric(m))
  }

  largest_sib <- max(dist)
  values <- numeric(length(dist))

  for (i in seq_len(length(dist))) {
    values[i] <- ceiling((n[i]*m)/10000000)
  }


  parts <- sum(values)
  splits <- numeric(parts)
  amount_of_sibs <- splits


  lower_bound <- 1
  upper_bound <- 0

  for(i in seq_len(length(values))) {
    x <- values[i]
    pers <- n[i]
    upper_bound <- upper_bound + x
    val <- rep(ceiling(pers / x), x - 1)
    splits[lower_bound:upper_bound] <- c(val, pers - sum(val))
    amount_of_sibs[lower_bound:upper_bound] <- rep(dist[i], upper_bound + 1 - lower_bound)
    lower_bound <- lower_bound + x
  }

  cusplits <- c(0, cumsum(splits))

  create_map(m, path) #Generate MAP-file

  MAFs <- runif(m, 0.01, 0.49) #Find MAFs for children and parents

  data.table::fwrite(data.table::as.data.table(MAFs), #Save the MAFs
         paste0(path, "MAFs.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)

  mu <- 2 * MAFs
  sigma <- sqrt(2 * MAFs * (1 - MAFs))
  critical <- qnorm(1 - k)

  causual_SNP <- sample.int(m, size = q, replace = F) #Create causal SNPs

  beta <- matrix(0, nrow = m, ncol = 1) #Make the beta matrix

  beta[causual_SNP] <- rnorm(q, 0, sqrt(hsq/q))


  data.table::fwrite(data.table::as.data.table(beta), #Save the beta
         paste0(path, "beta.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)


  header <- c("FID", "IID", "pheno", "child_lg", "child_liab", "par1_pheno", "par1_lg", "par1_liab", "par2_pheno", "par2_lg", "par2_liab")

  if(largest_sib != 0) {
    sib_header <- numeric(3*largest_sib)
    for(i in 0:(largest_sib - 1)) {
      sib_header[3 * i + 1] <- paste0("sib", i + 1, "_pheno", sep = "")
      sib_header[3 * i + 2] <- paste0("sib", i + 1, "_lg", sep = "")
      sib_header[3 * i + 3] <- paste0("sib", i + 1, "_liab", sep = "")
    }
    header <- c(header, sib_header, "line_pheno")
  }
  else {
    header <- c(header, "line_pheno")
  }

  #We create the header for the phenofile:
  data.table::fwrite(data.table::as.data.table(rbind(header)),
         paste0(path, "phenotypes.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)


  lock = tempfile()

  future.apply::future_lapply(1:parts, function(i) {

    parentmatrix <- parent_maker(m = m, number = 2*splits[i], MAFs)

    child <- t(vapply(seq(1, ncol(parentmatrix), 2), function(j) {
      child_mean <- rowMeans(parentmatrix[,j:(j+1), drop = FALSE])
      round_vec <- rbinom(n = m, 1, 1/2)
      snps_for_child <- rbinom(n = m, 2, 1/2)
      vals <- ifelse(child_mean == 1, 1, 0)
      child_mean[vals] <- ifelse(parentmatrix[,j][vals] == 1, snps_for_child, 1)
      return(dplyr::if_else(round_vec == 1, ceiling(child_mean), floor(child_mean)))},
      FUN.VALUE = numeric(m)))



    if(largest_sib != 0) {

      sibtable <- data.table::data.table("start" = numeric(splits[i]))

      if (amount_of_sibs[i] != 0) {
        for(k in seq_len(amount_of_sibs[i])) {
          sibs <- t(vapply(seq(1, ncol(parentmatrix), 2), function(j){
            sibs <- rowMeans(parentmatrix[,j:(j+1), drop = FALSE])
            round_vec <- rbinom(n = m, 1, 1/2)
            snps_for_sibs <- rbinom(n = m, 2, 1/2)
            vals <- ifelse(sibs == 1, 1, 0)
            sibs[vals] <- ifelse(parentmatrix[,j][vals] == 1, snps_for_sibs, 1)
            return(dplyr::if_else(round_vec == 1, ceiling(sibs), floor(sibs)))},
            FUN.VALUE = numeric(m)))

          sib_lg <- sweep(sweep(sibs, 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta
          sib_liab <- sib_lg + rnorm(splits[i], 0, sqrt(1 - hsq))
          sib_pheno <- sapply(sib_liab, function(x) ifelse(x > critical, 2, 1))


          sibtable <- cbind(sibtable, sib_pheno, sib_lg, sib_liab)
        }}

      if (largest_sib != amount_of_sibs[i]) {
        Naval <- matrix(nrow = splits[i], ncol = 3 * (largest_sib - amount_of_sibs[i]),
                        data = -9)
        sibtable <- cbind(sibtable, Naval)
        #This should "fill" the empty spaces for if the amount of sibs < largest_sib
      }
      sibtable <- sibtable[,start := NULL]
    }



    parlg <- sweep(sweep(t(parentmatrix), 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta
    parliab <- parlg + rnorm(2*splits[i], 0, sqrt(1 - hsq))

    c_lg <- sweep(sweep(child, 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta


    c_liab <- c_lg + rnorm(splits[i], 0, sqrt(1 - hsq))



    #Make phenotypes:
    c_pheno <- vapply(c_liab, function(x) ifelse(x >= critical, 2, 1), 
                      FUN.VALUE = matrix(splits[i]))

    p_pheno <- vapply(parliab, function(x) ifelse(x >= critical, 2, 1), 
                      FUN.VALUE = matrix(splits[i]))

    c_line_pheno <- c_pheno + 1

    #FID for the children/parents:

    id <- matrix(c((cusplits[i]+1):cusplits[i+1]))


    locked = flock::lock(lock) # locks file

    data.table::fwrite(data.table::as.data.table(to_ped(child, part = cusplits[i])), ## writes to locked file
           paste0(path, "genotypes.ped", sep = ""),
           quote = F,
           sep = " ",
           col.names = F,
           append = T)

    if (largest_sib != 0) {
      data.table::fwrite(data.table::as.data.table(cbind(id, rep(1, splits[i]), c_pheno, c_lg, c_liab, p_pheno[seq(1, 2*splits[i], 2)], parlg[seq(1, 2*splits[i], 2)], parliab[seq(1, 2*splits[i], 2)], p_pheno[seq(2, 2*splits[i], 2)], parlg[seq(2, 2*splits[i], 2)], parliab[seq(2, 2*splits[i], 2)], sibtable, c_line_pheno)),
             paste0(path, "phenotypes.txt", sep = ""),
             quote = F,
             sep = " ",
             col.names = F,
             append = T)}
    else {
      data.table::fwrite(data.table::as.data.table(cbind(id, rep(1, splits[i]), c_pheno, c_lg, c_liab, p_pheno[seq(1, 2*splits[i], 2)], parlg[seq(1, 2*splits[i], 2)], parliab[seq(1, 2*splits[i], 2)], p_pheno[seq(2, 2*splits[i], 2)], parlg[seq(2, 2*splits[i], 2)], parliab[seq(2, 2*splits[i], 2)], c_line_pheno)),
             paste0(path, "phenotypes.txt", sep = ""),
             quote = F,
             sep = " ",
             col.names = F,
             append = T)
    }



    flock::unlock(locked) #unlocks file

  }, future.seed = T)
}
