#'
#' @title Simulation with fixed family size
#' @md
#' @description Simulate genetic data, including genotypes, phenotype status
#' and liabilities, for individuals and their family, where each individual
#' has a fixed number of siblings specified by the `sib` parameter.
#' 
#'
#' @details
#' Parents' genotypes are simulated and used for creating the genotypes of
#' the individuals and their siblings. For the methodology behind the
#' simulation, see `vignette("liability-distribution")`.\cr
#' \code{sim_fixed_family} makes use of parallel computation in order to
#' decrease the running time. As one CPU core is left unused, the user
#' should be able to do other work while the simulation is running.
#'
#' @param n number of genotypes (individuals).
#' @param m number of SNPS per genotype.
#' @param q number of causal SNPs, i.e. SNPs that effect chances of having
#' the phenotype.
#' @param hsq squared heritability parameter.
#' @param k prevalence of phenotype.
#' @param path directory where the files will be stored. If nothing is
#' specified, \code{sim_fixed_family} writes its files in the current
#' working directory.
#' @param sib number of siblings per individual.
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


sim_fixed_family <- function(n, m, q, hsq, k, sib = 0, path = "") {

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
            "sib needs to be a non-negative integer" =
              (sib >= 0 && is.numeric(sib) && round(sib) == sib 
               && length(sib) == 1),
            "path needs to be default or a valid path ending with '/' or '\\\\'"
            = (path == "" || (dir.exists(path))
               && (substr(path, nchar(path), nchar(path)) == "/" ||
                     substr(path, nchar(path), nchar(path)) == "\\")))

  path <- path_validation(path)

  # Set worker nodes:
  future::plan(future::multiprocess, workers = max(future::availableCores(logical = F) - 1, 1))

  # Defining a function that creates genotypes for parents
  parent_maker <- function(m, number, MAFs) {
    vapply(1:number, function(y) {rbinom(m, 2, MAFs)}, FUN.VALUE = numeric(m))
  }

  # Here we hard code the function to simulate 10 million SNPs per session
  parts <- ceiling((n * m) / 10000000)

  # Splitting the work on the number of workers
  splits <- c(rep(ceiling(n / parts), parts - 1),
              n - sum(rep(ceiling(n / parts), parts - 1)))

  # Making ids for the different cores.
  cusplits <- c(0, cumsum(splits))

  create_map(m, path) # Generate MAP-file

  MAFs <- runif(m, 0.01, 0.49) # Find MAFs for children and parents
  data.table::fwrite(data.table::as.data.table(MAFs),
         paste0(path, "MAFs.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)

  mu <- 2 * MAFs # Find mean per SNP
  sigma <- sqrt(2 * MAFs * (1 - MAFs)) # Find standard deviation per SNP
  critical <- qnorm(1 - k) # Calculate the critical value

  causual_SNP <- sample.int(m, size = q, replace = F) # Create causal SNPs

  beta <- matrix(0, nrow = m, ncol = 1) # Make the vector of true effects

  beta[causual_SNP] <- rnorm(q, 0, sqrt(hsq / q))

  data.table::fwrite(data.table::as.data.table(beta),
         paste0(path, "beta.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)


  header <- c("FID", "IID", "pheno", "child_lg", "child_liab", "par1_pheno",
              "par1_lg", "par1_liab", "par2_pheno", "par2_lg", "par2_liab")

  if (sib != 0) {
    sib_header <- numeric(3 * sib)
    for (i in 0:(sib - 1)) {
      sib_header[3 * i + 1] <- paste0("sib", i + 1, "_pheno", sep = "")
      sib_header[3 * i + 2] <- paste0("sib", i + 1, "_lg", sep = "")
      sib_header[3 * i + 3] <- paste0("sib", i + 1, "_liab", sep = "")
    }
    header <- c(header, sib_header, "line_pheno")
  }
  else {
    header <- c(header, "line_pheno")
  }

  # Create the header for the phenofile:
  data.table::fwrite(data.table::as.data.table(rbind(header)),
         paste0(path, "phenotypes.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)


  lock <- tempfile()

  future.apply::future_lapply(1:parts, function(i) {
    # Create the parents
    parentmatrix <- parent_maker(m = m, number = 2 * splits[i], MAFs)

    # Using the parents' genotypes the genotypes can be calculated for children
    child <- t(vapply(seq(1, ncol(parentmatrix), 2), function(j){
      child_mean <- rowMeans(parentmatrix[,j:(j+1), drop = FALSE])
      round_vec <- rbinom(n = m, 1, 1/2)
      snps_for_child <- rbinom(n = m, 2, 1/2)
      vals <- ifelse(child_mean == 1, 1, 0)
      child_mean[vals] <- ifelse(parentmatrix[,j][vals] == 1, snps_for_child, 1)
      return(dplyr::if_else(round_vec == 1, ceiling(child_mean), floor(child_mean)))},
      FUN.VALUE = numeric(m)))

    if (sib != 0) {

      sibtable <- data.table::data.table("begin" = numeric(splits[i]))

      for (k in seq_len(sib)) {
        # We create genotypes for siblings in same way as for the individuals
        sibs <- t(vapply(seq(1, ncol(parentmatrix), 2), function(j){
          sibs <- rowMeans(parentmatrix[,j:(j+1), drop = FALSE])
          round_vec <- rbinom(n = m, 1, 1/2)
          snps_for_sibs <- rbinom(n = m, 2, 1/2)
          vals <- ifelse(sibs == 1, 1, 0)
          sibs[vals] <- ifelse(parentmatrix[,j][vals] == 1, snps_for_sibs, 1)
          return(dplyr::if_else(round_vec == 1, ceiling(sibs), floor(sibs)))},
          FUN.VALUE = numeric(m)))

        # Find the genetic liability, liability and phenotype for the
        # siblings and add this to the "sibtable"
        sib_lg <- sweep(sweep(sibs, 2, mu, FUN = "-"),
                        2, sigma, FUN = "/") %*% beta
        sib_liab <- sib_lg + rnorm(splits[i], 0, sqrt(1 - hsq))
        sib_pheno <- sapply(sib_liab, function(x) ifelse(x > critical, 2, 1))


        sibtable <- cbind(sibtable, sib_pheno, sib_lg, sib_liab)


      }
      # Remove the begin column in the table.
      sibtable <- sibtable[, "begin" := NULL]

    }



    parlg <- sweep(sweep(t(parentmatrix), 2, mu, FUN = "-"),
                   2, sigma, FUN = "/") %*% beta
    parliab <- parlg + rnorm(2 * splits[i], 0, sqrt(1 - hsq))

    c_lg <- sweep(sweep(child, 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta

    c_liab <- c_lg + rnorm(splits[i], 0, sqrt(1 - hsq))



    # Make phenotypes:
    c_pheno <- vapply(c_liab, function(x) ifelse(x >= critical, 2, 1), 
                      FUN.VALUE = matrix(splits[i]))
    
    p_pheno <- vapply(parliab, function(x) ifelse(x >= critical, 2, 1), 
                      FUN.VALUE = matrix(splits[i]))

    c_line_pheno <- c_pheno + 1

    # Create the ID per individual
    id <- matrix(c((cusplits[i] + 1):cusplits[i + 1]))


    locked <- flock::lock(lock) # locks file

    # writes to locked file
    data.table::fwrite(data.table::as.data.table(to_ped(child, part = cusplits[i])),
           paste0(path, "genotypes.ped", sep = ""),
           quote = F,
           sep = " ",
           col.names = F,
           append = T)


    if(sib != 0) {
      data.table::fwrite(data.table::as.data.table(cbind(id, rep(1, splits[i]), c_pheno, c_lg, c_liab,
                                 p_pheno[seq(1, 2 * splits[i], 2)],
                                 parlg[seq(1, 2 * splits[i], 2)],
                                 parliab[seq(1, 2 * splits[i], 2)],
                                 p_pheno[seq(2, 2 * splits[i], 2)],
                                 parlg[seq(2, 2 * splits[i], 2)],
                                 parliab[seq(2, 2 * splits[i], 2)],
                                 sibtable, c_line_pheno)),
             paste0(path, "phenotypes.txt", sep = ""),
             quote = F,
             sep = " ",
             col.names = F,
             append = T)
    }
    else {
      data.table::fwrite(data.table::as.data.table(cbind(id, rep(1, splits[i]), c_pheno, c_lg, c_liab,
                                 p_pheno[seq(1, 2 * splits[i], 2)],
                                 parlg[seq(1, 2 * splits[i], 2)],
                                 parliab[seq(1, 2 * splits[i], 2)],
                                 p_pheno[seq(2, 2 * splits[i], 2)],
                                 parlg[seq(2, 2 * splits[i], 2)],
                                 parliab[seq(2, 2 * splits[i], 2)],
                                 c_line_pheno)),
             paste0(path, "phenotypes.txt", sep = ""),
             quote = F,
             sep = " ",
             col.names = F,
             append = T)
    }


    flock::unlock(locked) #unlocks file

  }, future.seed = T)
}
