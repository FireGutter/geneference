#'
#' @title Simulation with a fixed family size
#' @md
#' @description Simulate genotypes for individuals with family history, where
#' each individual has a fixed number of siblings. Parents' genotypes are
#' simulated and used for simulating the genotypes of individuals and their
#' siblings.
#' \code{family_simulation} makes use of parallel computation in order to
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
#' @param path directory where the files will be stored. If nothing is
#' specified, \code{family_simulation} writes its files in the current
#' working directory.
#' @param sib number of siblings per individual.
#'
#' @importFrom data.table :=
#' @import stats
#'
#' @return Does not return any value, but prints the following five files to
#' the \code{path} parameter specified in the function call:
#' * Three text files:
#'     * beta.txt - a file of \code{m} rows with one column. The i'th row is
#'     the true effect of the i'th SNP.
#'     * MAFs.txt - a file of \code{m} rows with one column. The i'th row is
#'     the true Minor Allelle Frequence of the i'th SNP.
#'     * phenotypes.txt - a file of \code{n} rows, number of columns depend on
#'     number of siblings. The file contains the phenotype and liability of
#'     each individual as well as information on the liabilities and phenotype
#'     status of their parents and siblings.
#' * genotypes.map - a file created such that PLINK will work with the genotype
#' data.
#' * genotypes.ped - the simulated genotypes in a PLINK-readable format.
#'
#' @export
#' @examples
#' family_simulation(n = 10000, m = 10000, q = 100, hsq = 0.5, k = 0.05, sib = 3, path = "")

family_simulation <- function(n, m, q, hsq, k, sib = 0, path = "") {
  stopifnot("n needs to be an integer greater than 0" =
              (n > 0 && class(n) == "numeric" && n == round(n)),
            "m needs to be an integer greater than 0" =
              (m > 0 && class(m) == "numeric" && m == round(m)),
            "q needs to be an integer greater than 0 and smaller than m" =
              (q > 0 && class(q) == "numeric" && q == round(q) && q <= m),
            "hsq needs to be a number between 0 and 1" =
              (hsq > 0 && hsq < 1 && class(hsq) == "numeric"),
            "k needs to be a number between 0 and 1" =
              (k > 0 && k < 1 && class(k) == "numeric"),
            "sib needs to be a non-negative integer" =
              (sib >= 0 && class(sib) == "numeric" && round(sib) == sib),
            "path needs to be default or a valid path ending with '/' or '\\\\'"
            = (path == "" || (dir.exists(path))
               && (substr(path, nchar(path), nchar(path)) == "/" ||
                     substr(path, nchar(path), nchar(path)) == "\\")))

  # Set worker nodes:
  future::plan(future::multiprocess, workers = max(future::availableCores(logical = F) - 1, 1))
  
  # Defining a function that creates genotypes for parents
  parent_maker <- function(m, number, MAFs) {
    sapply(1:number, function(i) {rbinom(m, 2, MAFs)})
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
    header <- c(header, sib_header)
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
    child <- t(sapply(seq(1, ncol(parentmatrix), 2), function(i){
      child_mean <- rowMeans(parentmatrix[,i:(i+1), drop = FALSE])
      round_vec <- rbinom(n = m, 1, 1/2)
      snps_for_child <- rbinom(n = m, 2, 1/2)
      vals <- ifelse(child_mean == 1, 1, 0)
      child_mean[vals] <- ifelse(parentmatrix[,i][vals] == 1, snps_for_child, 1)
      return(dplyr::if_else(round_vec == 1, ceiling(child_mean), floor(child_mean)))}))

    if (sib != 0) {

      sibtable <- data.table::data.table("begin" = numeric(splits[i]))

      for (j in 1:sib) {
        # We create genotypes for siblings in same way as for the individuals
        sibs <- t(sapply(seq(1, ncol(parentmatrix), 2), function(i){
          sibs <- rowMeans(parentmatrix[,i:(i+1), drop = FALSE])
          round_vec <- rbinom(n = m, 1, 1/2)
          snps_for_sibs <- rbinom(n = m, 2, 1/2)
          vals <- ifelse(sibs == 1, 1, 0)
          sibs[vals] <- ifelse(parentmatrix[,i][vals] == 1, snps_for_sibs, 1)
          return(dplyr::if_else(round_vec == 1, ceiling(sibs), floor(sibs)))}))

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
    c_pheno <- sapply(c_liab, function(x) ifelse(x > critical, 2, 1))
    p_pheno <- sapply(parliab, function(x) ifelse(x > critical, 2, 1))

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
                                 sibtable)),
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
                                 parliab[seq(2, 2 * splits[i], 2)])),
             paste0(path, "phenotypes.txt", sep = ""),
             quote = F,
             sep = " ",
             col.names = F,
             append = T)
    }
    
    
    flock::unlock(locked) #unlocks file

  }, future.seed = T)
}
