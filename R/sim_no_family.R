#'
#' @title Simulation without family history
#' @md
#' @description Simulate genetic data, including genotypes, phenotype status
#' and liabilities, for individuals.
#'
#' @details
#' As this function does not include family history, its resulting data cannot
#' be used by \code{assign_ltfh_phenotype()} or
#' \code{assign_GWAX_phenotype()}.\cr
#' For the methodology behind the simulation, see
#' `vignette("liability-distribution")`.\cr
#' \code{sim_no_family()} makes use of parallel computation in order to
#' decrease the running time. As at least one CPU core is left unused, the user
#' should be able to do other work while the simulation is running.
#'
#'
#' @param n number of genotypes (individuals).
#' @param m number of SNPS per genotype.
#' @param q number of causal SNPs, i.e. SNPs that effect chances of having
#' the phenotype.
#' @param hsq squared heritability parameter.
#' @param k prevalence of phenotype.
#' @param path directory where the files will be stored. If nothing is
#' specified, \code{sim_no_family} writes its files in the current
#' working directory.
#'
#' @return Does not return any value, but prints the following five files to
#' the \code{path} parameter specified in the function call:
#' * Three text files:
#'     * beta.txt - a file of \code{m} rows with one column. The i'th row is
#'     the true effect of the i'th SNP.
#'     * MAFs.txt - a file of \code{m} rows with one column. The i'th row is
#'     the true Minor Allelle Frequency of the i'th SNP.
#'     * phenotypes.txt - a file of \code{n} rows. The file contains the
#'     phenotype status and liability of each individual.
#' * genotypes.map - a file created such that PLINK will work with the genotype
#' data.
#' * genotypes.ped - the simulated genotypes in a PLINK-readable format.
#'
#' @section Warning:
#' Simulating large datasets takes time and generates large files. For details
#' on time complexity and required disk space, see
#' `vignette("sim-benchmarks")`.\cr
#' The largest file generated is `genotypes.ped`. See `convert_geno_file()` to convert it
#' to another file format, thereby reducing its size significantly.
#'
#' @import stats
#'
#' @export



sim_no_family <- function(n, m, q, hsq, k, path){
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
            "path needs to be default or a valid path ending with '/' or '\\\\'"
            = (path == "" || (dir.exists(path))
               && (substr(path, nchar(path), nchar(path)) == "/" ||
                     substr(path, nchar(path), nchar(path)) == "\\")))
  
  path <- path_validation(path)

  future::plan(future::multisession(workers = 2))

  SNP_simulation <- function(n, m, MAFs){
    #Function that draws from a binomial distribution and uses the MAF vector of probabilities
    #Returns a list of length n with entires m integers.
    return(lapply(1:n, function(y) rbinom(m, 2, MAFs)))
  }

  simulation <- function(n, m, MAFs){
    #Detecting the number of cores of ones computer and making one core available
    cores <- parallel::detectCores(logical = FALSE)-1

    #If the cores are creater than 1 then we check if n is a multiple of cores
    if (cores>1){
      while (n %% cores != 0){
        cores <- cores-1
      }
    } else {
      cores <- 1
    }

    #Making a cluster of the cores
    cl <- parallel::makeCluster(cores)
    #Exporting the cores to the function we want them to run on
    parallel::clusterExport(cl, c("SNP_simulation"))

    #Making the multiprocessing simulation (note: because of the way rbindlist works we have to switch rows and columns when inputting it to 'SNP_simulation'). Can maybe also use bind_cols() instead
    geno <- parallel::clusterCall(cl, function(persons, SNPs, MAF) SNP_simulation(persons, SNPs, MAFs), persons = (n/cores), SNPs = m, MAF = MAFs)

    parallel::stopCluster(cl)

    #binding it all together
    return(t(dplyr::bind_cols(geno)))
  }

  create_map(n, path)

  #Calculating the number MAF probabilities so that each core has the same probabilities
  MAFs <- runif(m, 0.01, 0.49)

  parts <- ceiling((n*m)/10000000) #varaiable to control how many parts we divide the problem into

  data.table::fwrite(data.table::as.data.table(MAFs),
         paste(path,"MAFs.txt",sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)

  causual_SNP <- sample.int(m, size = q, replace = F)

  beta <- matrix(0, nrow = m, ncol = 1)
  beta[causual_SNP] <- rnorm(q, 0, sqrt(hsq/q))

  data.table::fwrite(data.table::as.data.table(beta),
         paste(path,"beta.txt",sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)

  data.table::fwrite(data.table::as.data.table(rbind(c("FID", "IID", "pheno", "line_pheno"))),
         paste(path,"phenotypes.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)

  lock = tempfile()

  future.apply::future_lapply(1:parts, function(i) {
    #for (i in 1:parts) {
    if (i == parts & n%%parts != 0){
      batch_size = n%%parts
    } else {
      batch_size = ceiling(n/parts)
    }

    locked = flock::lock(lock) # locks file

    persons <- simulation(batch_size, m, MAFs)

    mu <- 2*MAFs
    sigma <- sqrt(2*MAFs*(1-MAFs))

    pliab <- sweep(sweep(persons, 2, mu, FUN = "-"), 2, sigma, FUN = "/")

    liability <- pliab %*% beta + rnorm(batch_size, 0, sqrt(1-hsq))

    critical <- qnorm(1-k)

    y <- sapply(liability, function(x) ifelse(x >= critical, 2, 1))
    line_pheno <- y + 1

    id <- matrix(1:(batch_size))+(batch_size)*(i-1)

    data.table::fwrite(data.table::as.data.table(to_ped(persons, (i-1) * (batch_size))), ## writes to locked file
           paste(path,"genotypes.ped", sep = ""),
           quote = F,
           sep = " ",
           col.names = F,
           append = T)

    data.table::fwrite(data.table::as.data.table(cbind(id,rep(1,batch_size), y, line_pheno)),
           paste(path,"phenotypes.txt",sep = ""),
           quote = F,
           sep = " ",
           col.names = F,
           append = T)

    flock::unlock(locked) #unlocks file

  }, future.seed = T)
}
