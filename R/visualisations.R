#' @title Plot a QQ-plot of p-values
#'
#' @description description to be written. Assumes p-values are uniformly
#' distributed.
#' \code{vignette("ggplot2-specs", package = "ggplot2")} details the usage of
#' aesthetic parameters in ggplot2.
#'
#' @param dataset data imported to R by \code{load_assoc_results()}.
#' @param line_color color of the identity line.
#' @param line_size size of the identity line.
#' @param qq_color color of the points showing the quantiles of the p-values.
#' @param qq_shape shape of the points showing the quantiles of the p-values.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#' ÆNDRES: Must be either pdf, png or jpeg
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @import stats
#'
#' @export
plot_pval_QQ <- function(dataset,
                         line_color="black",
                         line_size = 1,
                         qq_color = "cornflowerblue",
                         qq_shape = 19,
                         save_plot_path = FALSE,
                         plot_filename = "QQ-pvals.png") {
  
  stopifnot("dataset must have a column named 'P'" = "P" %in% colnames(dataset),
            "line_size needs to be a positive number" =
              (is.numeric(line_size) && line_size > 0 &&
                 length(line_size) == 1),
            "save_plot_path needs to be default or a valid path" =
              (save_plot_path == FALSE || dir.exists(save_plot_path)),
            "plot_filename must have either '.png', '.pdf' or '.jpeg' as extension" =
              (file_ext(plot_filename) == "png" ||
                 file_ext(plot_filename) == "pdf" ||
                 file_ext(plot_filename) == "jpeg"))
  
  plt <- ggplot2::ggplot(data = dataset) +
    ggplot2::geom_abline(size = line_size, color = line_color) +
    ggplot2::geom_qq(mapping = ggplot2::aes(sample = P),
                     distribution = qunif,
                     color = qq_color,
                     shape = qq_shape) +
    ggplot2::labs(x = "Uniform Quantiles",
                  y = "P-values",
                  title = "Uniform QQ-plot",
                  subtitle = deparse(substitute(dataset))) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(face = "bold",
                                                         hjust = 0.5))
  if (save_plot_path != FALSE) {
    ggplot2::ggsave(filename = plot_filename,
                    plot = plt,
                    path = save_plot_path)
  } else{
    return(plt)
  }
}

#' @title Plot a histogram of p-values
#'
#' @description write something
#' \code{vignette("ggplot2-specs", package = "ggplot2")} details the usage of
#' aesthetic parameters in ggplot2.
#'
#' @param dataset data imported to R by \code{load_assoc_results()}.
#' @param bins number of bins used.
#' @param line_color outline color of bins.
#' @param fill_color fill color of bins.
#' @param mean_color color of the line showing the mean value.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @export
plot_pval_hist <- function(dataset,
                           bins = 20,
                           line_color = "black",
                           fill_color = "cornflowerblue",
                           mean_color = "red",
                           save_plot_path = FALSE,
                           plot_filename = "pvalue_histogram.png") {
  
  stopifnot("dataset must have a column named 'P'" = "P" %in% colnames(dataset),
            "bins needs to be a positive integer" =
              (is.numeric(bins) && bins > 0 && bins == round(bins) &&
                 length(bins) == 1),
            "save_plot_path needs to be default or a valid path" =
              (save_plot_path == FALSE || dir.exists(save_plot_path)),
            "plot_filename must have either '.png', '.pdf' or '.jpeg' as extension" =
              (file_ext(plot_filename) == "png" ||
                 file_ext(plot_filename) == "pdf" ||
                 file_ext(plot_filename) == "jpeg"))
  
  plt <- ggplot2::ggplot(data = dataset) +
    ggplot2::geom_histogram(mapping = ggplot2::aes(P),
                            bins = bins + 1,
                            color = line_color,
                            fill = fill_color,
                            boundary = 0) +
    ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = length(P) / bins),
                        color = mean_color) +
    ggplot2::labs(x = "P-values",
                  y = "Frequency",
                  title = "P-value Histogram",
                  subtitle = deparse(substitute(dataset))) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(face = "bold",
                                                         hjust = 0.5))
  if (save_plot_path != FALSE) {
    ggplot2::ggsave(filename = plot_filename,
                    plot = plt,
                    path = save_plot_path)
  } else{
    return(plt)
  }
}

#' @title Manhattan plot of p-values
#'
#' @description description to be written. Plots a Manhattan plot with only the
#' SNPs that are significant at a 5%-significance level.
#' Write a few sentences about why the plot is different from real world plots
#' (linkage disequilibrium)
#'
#' @param dataset data imported to R by \code{load_assoc_results()}.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @export
plot_manhattan <- function(dataset,
                           save_plot_path = FALSE,
                           plot_filename = "manhattan_plot.png") {
  
  stopifnot("dataset must have a column named 'P', 'SNP' and 'causal'" =
              all(c("SNP", "P", "causal") %in% colnames(dataset)),
            "save_plot_path needs to be default or a valid path" =
              (save_plot_path == FALSE || dir.exists(save_plot_path)),
            "plot_filename must have either '.png', '.pdf' or '.jpeg' as extension" =
              (file_ext(plot_filename) == "png" ||
                 file_ext(plot_filename) == "pdf" ||
                 file_ext(plot_filename) == "jpeg"))
  
  data_plt <- dataset %>%
    dplyr::filter(P < 0.05)
  m <- nrow(dataset)
  plt <- ggplot2::ggplot(data_plt) +
    ggplot2::geom_point(mapping = ggplot2::aes(SNP,
                                               -log10(P),
                                               colour = causal)) +
    ggplot2::theme_light() +
    ggplot2::geom_segment(ggplot2::aes(x = 1, xend = m, y = -log10(0.05),
                                       yend = -log10(0.05)),
                          colour = "cornflowerblue",
                          linetype = "dashed",
                          size = 1.25) +
    ggplot2::geom_segment(ggplot2::aes(x = 1,
                                       xend = m,
                                       y = -log10(0.05 / m),
                                       yend = -log10(0.05 / m)),
                          colour = "cornflowerblue",
                          linetype = "dashed", size = 1) +
    ggplot2::labs(title = "Manhattan plot",
                  subtitle = deparse(substitute(dataset)),
                  colour = "SNP is causal")

  if (save_plot_path != FALSE) {
    ggplot2::ggsave(filename = plot_filename,
                    plot = plt,
                    path = save_plot_path)
  } else{
    return(plt)
  }
}

#' @title Plot true effects against estimated effects of SNPs
#'
#' @description description to be written.
#'
#' @param dataset data imported to R by \code{load_assoc_results()}.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @export
plot_estimates_vs_true <- function(dataset,
                                   save_plot_path = FALSE,
                                   plot_filename = "beta_comparison.png") {
  
  stopifnot("dataset must have a column named 'BETA', 'true_effect' and 'bonferroni'" =
              all(c("BETA", "true_effect", "bonferroni") %in% colnames(dataset)),
            "save_plot_path needs to be default or a valid path" =
              (save_plot_path == FALSE || dir.exists(save_plot_path)),
            "plot_filename must have either '.png', '.pdf' or '.jpeg' as extension" =
              (file_ext(plot_filename) == "png" ||
               file_ext(plot_filename) == "pdf" ||
               file_ext(plot_filename) == "jpeg"))
  
  tmpdataacc <- dataset %>%
    dplyr::filter(significant)

  plt <- ggplot2::ggplot(tmpdataacc) +
    ggplot2::geom_point(mapping = ggplot2::aes(BETA,
                                               true_effect,
                                               colour = bonferroni)) +
    ggplot2::geom_abline(slope = 1) +
    ggplot2::labs(x = "Estimated slope coefficient",
                  y = "True effect",
                  title = "Scatterplot of true effects against estimates",
                  subtitle = deparse(substitute(dataset)),
                  colour = "Significant with\n Bonferroni-correction") +
    ggplot2::theme_light() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(face = "bold",
                                                         hjust = 0.5))

  if (save_plot_path != FALSE) {
    ggplot2::ggsave(filename = plot_filename,
                    plot = plt,
                    path = save_plot_path)
  } else{
    return(plt)
  }
}

#' @title Plot true genetic liabilities against posterior mean genetic
#' liabilities
#'
#' @description write description.
#' \code{vignette("ggplot2-specs", package = "ggplot2")} details the usage of
#' aesthetic parameters in ggplot2.
#'
#' @param dataset data imported to R by \code{load_assoc_results()}.
#' @param line_color color of identity line.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @export
plot_pmgl_vs_true <- function(dataset,
                              line_color = "black",
                              save_plot_path = FALSE,
                              plot_filename = "posterior_liabilities.png") {
  
  stopifnot("dataset must have a column named 'LTFH_pheno' and 'child_lg'" =
              all(c("LTFH_pheno", "child_lg") %in% colnames(dataset)),
            "save_plot_path needs to be default or a valid path" =
              (save_plot_path == FALSE || dir.exists(save_plot_path)),
            "plot_filename must have either '.png', '.pdf' or '.jpeg' as extension" =
              (file_ext(plot_filename) == "png" ||
                 file_ext(plot_filename) == "pdf" ||
                 file_ext(plot_filename) == "jpeg"))
  
  plt <- ggplot2::ggplot(dataset) +
    ggplot2::geom_point(ggplot2::aes(LTFH_pheno, child_lg),
                        color = "cornflowerblue",
                        size = 0.6) +
    ggplot2::geom_abline(color = line_color, linetype = "dashed") +
    ggplot2::theme_light() +
    ggplot2::labs(x = "Posterior mean genetic liability",
                  y = "True genetic liability",
                  title = "Scatterplot of liabilities and LT-FH phenotype",
                  subtitle = deparse(substitute(dataset))) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(face = "bold",
                                                         hjust = 0.5))

  if (save_plot_path != FALSE) {
    ggplot2::ggsave(filename = plot_filename,
                    plot = plt,
                    path = save_plot_path)
  } else{
    return(plt)
  }
}

#' @title Plot transformed linear GWAS beta values against the LTFH_transformed
#' values.
#'
#' @description description to be written.
#'
#' @param dataset data imported to R.
#' @param transformed_beta the transformed beta values from LTFH.
#' @param linear_beta the beta values from standard GWAS.
#' @param bonferroni corrected significant values from LTFH.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @export
GWAS_beta_VS_LTFH_transfored_beta <- function(dataset,
                                              transformed_beta,
                                              linear_beta,
                                              bonferroni,
                                              save_plot_path = FALSE,
                                              plot_filename = "LTFH.png") {
  
  stopifnot("save_plot_path needs to be default or a valid path" =
              (save_plot_path == FALSE || dir.exists(save_plot_path)),
            "plot_filename must have either '.png', '.pdf' or '.jpeg' as extension" =
              (tools::file_ext(plot_filename) == "png" ||
                 tools::file_ext(plot_filename) == "pdf" ||
                 tools::file_ext(plot_filename) == "jpeg"))
  
  plt <- ggplot2::ggplot(dataset, ggplot2::aes(x = transformed_beta,
                                               y = linear_beta)) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = bonferroni)) +
    ggplot2::geom_smooth(se = FALSE, method = "lm")+
    ggplot2::geom_abline(intercept = 0, slope = 1, color="green", 
                         linetype="dashed", size = 1) +
    ggplot2::labs(x = "LTFH beta",
                  y = "GWAS beta",
                  title = "Comparing LTFH beta to GWAS beta",
                  subtitle = deparse(substitute(dataset)),
                  colour = "Significant with\n Bonferroni-correction") +
    ggplot2::theme_light() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(face = "bold",
                                                         hjust = 0.5))
  
  if (save_plot_path != FALSE) {
    ggplot2::ggsave(filename = plot_filename,
                    plot = plt,
                    path = save_plot_path)
  } else{
    return(plt)
  }
}

