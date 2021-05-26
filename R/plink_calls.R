#' @title Convert .ped files to .bed
#'
#' @description to be written. Essentially:
#' calls plink to convert for storing with smaller file size.
#' NOTE: .ped file can be viewed, .bed cannot be viewed
#'
#' @param ped_file path of the .ped file to be converted, excluding file
#' extension.
#' @param bed_file path of the .bed file to be written, excluding file
#' extension.
#' @param del if \code{TRUE}, deletes the \code{ped_file} after conversion.
#' @param plink_path \code{TRUE} if user has added PLINK to the system
#' variable "Path". Otherwise a string specifying the path to PLINK.
#'
#' @return Does not return anything, but writes the converted \code{ped_file}
#' to disk at the path \code{bed_file}.
#'
#' @export
p_2_b <- function(ped_file, bed_file=ped_file, del=TRUE, plink_path=TRUE) {
  stopifnot("ped_file needs to be a valid file without extension" =
              (file.exists(paste0(ped_file, ".ped")) 
              && file.exists(paste0(ped_file, ".map"))),
            "bed_file needs to be a valid file without extension" =
              (tools::file_ext(output_file) == ""),
            "del needs to be either TRUE or FALSE" = is.logical(del),
            "plink_path needs to be a valid path to plink" =
              (plink_path == TRUE || file.exists(paste0(plink_path, "/plink.exe"))))

  if (plink_path != TRUE) {
    tmp_path <- paste0("SET PATH=", path_to_path(plink_path), ";")
  } else {
    tmp_path <- ""
  }
  shell(paste(tmp_path, "&&", "plink --file", path_to_path(ped_file),
               "--make-bed --out", path_to_path(bed_file)))

  if (del == TRUE) {
    unlink(paste0(ped_file, ".ped"))
    unlink(paste0(ped_file, ".map"))
  }
}

#' @title Run association analysis
#'
#' @description to be written. calls PLINK.
#' (include specifying pheno_names from simulation:
#' GWAS use pheno, for LTFH use LTFH_pheno, GWAX use GWAX_pheno)
#' NOTE: creates a temporary .bed file if \code{bed} is \code{FALSE}. This
#' takes some additional time.
#' Also: File extension will be .assoc or .qassoc depending on the used model.
#'
#' @param geno_file string specifying path to genotypes file, including file
#' extension.
#' @param pheno_file string specifying path to phenotypes file, including file
#' extension.
#' @param pheno_name column name of phenotype to be used in analysis.
#' @param out_file string specifying path and name of the resulting output
#' file, excluding file extension.
#' @param bed logical indicating whether or not the genotypes file is a .bed
#' file. \code{FALSE} if the file is .ped.
#' @param plink_path \code{TRUE} if user has added PLINK to the system
#' variable "Path". Otherwise a string specifying the path to PLINK.
#'
#' @return Does not return anything, but PLINK writes \code{out_file} to disk,
#' containing results of the performed analysis. The columns in the printed
#' dataset depends on which phenotype is used for the analysis.
#'
#' @import tools
#'
#' @export
analysis_association <- function(geno_file, pheno_file, pheno_name, out_file,
                                 bed=TRUE, plink_path=TRUE) {
  stopifnot("geno_file needs to be a valid file" =
              (file.exists(geno_file) &&
                 ((tools::file_ext(geno_file) == "ped"
                   && file.exists(paste0(ped_file, ".map"))) ||
                    (tools::file_ext(geno_file) == "bed"
                    && file.exists(paste0(ped_file, ".bim"))
                    && file.exists(paste0(ped_file, ".fam"))))),
            "pheno_file needs to be a valid file" =
              (file.exists(pheno_file) &&
                 tools::file_ext(pheno_file) == "txt"),
            "pheno_name does not exist in 'pheno_file'" =
              (pheno_name %in% colnames(data.table::fread(pheno_file))),
            "out_file needs to be a valid file path without file extension" =
              tools::file_ext(output_file) == "",
            "bed needs to be either TRUE or FALSE" = is.logical(bed),
            "plink_path needs to be a valid path to plink" =
              (plink_path == TRUE ||
                 file.exists(paste0(plink_path, "/plink.exe"))))


  if (plink_path != TRUE) {
    tmp_path <- paste0("SET PATH=", path_to_path(plink_path), ";")
  } else{
    tmp_path <- ""
  }
  if (bed) {
    file_type <- "--bfile"
  } else{
    file_type <- "--file"
  }
  geno_file <- file_path_sans_ext(geno_file) #NOTE IF NOT NEW CONVENTION: this is from tools
  plink_command <- paste(tmp_path, "&&",
                         "plink", file_type, path_to_path(geno_file),
                         "--pheno", path_to_path(pheno_file),
                         "--pheno-name", pheno_name,
                         "--out", path_to_path(out_file),
                         "--assoc")

  shell(cmd = plink_command)
}

#' @title Run Lasso analysis
#'
#' @description to be written.
#' NOTE: File extension will be .lasso.
#'
#' @param geno_file string specifying path to genotypes file, including file
#' extension.
#' @param pheno_file string specifying path to phenotypes file, including file
#' extension.
#' @param pheno_name column name of phenotype to be used in analysis.
#' @param out_file string specifying path and name of the resulting output
#' file, excluding file extension.
#' @param bed logical indicating whether or not the genotypes file is a .bed
#' file. \code{FALSE} if the file is .ped.
#' @param plink_path \code{TRUE} if user has added PLINK to the system
#' variable "Path". Otherwise a string specifying the path to PLINK.
#' @param hsq heritability parameter (used for calibration of lasso).
#'
#' @return Does not return anything, but PLINK writes \code{out_file} to disk,
#' containing results of the performed lasso-regression. The columns in the
#' printed dataset depends on which phenotype is used for the analysis.
#'
#' @import tools
#'
#' @export
analysis_lasso <- function(geno_file, pheno_file, pheno_name,
                           out_file, bed=TRUE, plink_path=TRUE, hsq=0.5) {
  stopifnot("geno_file needs to be a valid file" =
              (file.exists(geno_file) &&
                 ((tools::file_ext(geno_file) == "ped"
                   && file.exists(paste0(ped_file, ".map"))) ||
                    (tools::file_ext(geno_file) == "bed"
                     && file.exists(paste0(ped_file, ".bim"))
                     && file.exists(paste0(ped_file, ".fam"))))),
            "pheno_file needs to be a valid file" =
              (file.exists(pheno_file) &&
                 tools::file_ext(pheno_file) == "txt"),
            "pheno_name does not exist in 'pheno_file'" =
              (pheno_name %in% colnames(data.table::fread(pheno_file))),
            "out_file needs to be a valid file path without file extension" =
              tools::file_ext(output_file) == "",
            "bed needs to be either TRUE or FALSE" = is.logical(bed),
            "plink_path needs to be a valid path to plink" =
              (plink_path == TRUE ||
                 file.exists(paste0(plink_path, "/plink.exe"))),
            "hsq needs to a number between 0 and 1" =
              (is.numeric(hsq) && length(hsq) == 1 && 0 < hsq && hsq < 1))

  if (plink_path != TRUE) {
    tmp_path <- paste0("set PATH=%PATH%;", plink_path, ";")
  } else{
    tmp_path <- ""
  }
  if (bed) {
    file_type <- "--bfile"
  } else{
    file_type <- "--file"
  }
  geno_file <- file_path_sans_ext(geno_file)
  plink_command <- paste(tmp_path, "plink", file_type, geno_file,
                         "--pheno", pheno_file,
                         "--pheno-name", pheno_name,
                         "--out", out_file,
                         "--lasso", hsq)

  system(command = plink_command)
}
