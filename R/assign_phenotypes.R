#'
#'



#'
#'





  }

  # Set GWAX phenotype to maximum of the phenotypes in the family
  pheno[, "GWAX_pheno"] <- apply(pheno[, pheno_cols], 1, max, na.rm = TRUE)

  # Write to file
  data.table::fwrite(x = pheno, file = output_file,
                     quote = F,
                     sep = " ",
                     col.names = T,
                     append = F)
}
