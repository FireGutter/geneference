#' @title Plot a QQ-plot of p-values
#'
#' @description Creates a QQ-plot of the p-values specified.
#' 
#' @details
#' Use \code{analysis_association()} to obtain p-values for an analysis and
#' import these with \code{load_results()}.\cr
#' Assumes p-values are uniformly distributed, and hence uses the uniform
#' distribution for the theoretical quantiles.\cr
#' \code{vignette("ggplot2-specs", package = "ggplot2")} details the usage of
#' aesthetic parameters in ggplot2.
#'
#' @param dataset data imported to R by \code{load_results()}.
#' @param P column containing p-values to be plotted. State in the format
#' dataset$P_value.
#' @param line_color color of the identity line.
#' @param line_size size of the identity line.
#' @param qq_color color of the points showing the quantiles of the p-values.
#' @param qq_shape shape of the points showing the quantiles of the p-values.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#' Must be either .pdf, .png or .jpeg
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @import stats
#'
#' @export
plot_pval_QQ <- function(dataset,
                         P,
                         line_color="black",
                         line_size = 1,
                         qq_color = "cornflowerblue",
                         qq_shape = 20,
                         save_plot_path = FALSE,
                         plot_filename = "QQ-pvals.png") {
  
  stopifnot("'P' must be a column in 'dataset'"  = sub(".*\\$", "", deparse(substitute(P))) %in% colnames(dataset),
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
#' @description
#' Creates a histogram of the specified p-values.
#'
#' @details
#' Use \code{analysis_association()} to obtain p-values for an analysis and
#' import these with \code{load_results()}.\cr
#' The function also draws a line on the y-axis showing the mean of the
#' p-values.\cr
#' \code{vignette("ggplot2-specs", package = "ggplot2")} details the usage of
#' aesthetic parameters in ggplot2.
#'
#' @param dataset data imported to R by \code{load_results()}.
#' @param P column containing p-values to be plotted. State in the format
#' dataset$P_value.
#' @param bins number of bins used.
#' @param line_color outline color of bins.
#' @param fill_color fill color of bins.
#' @param mean_color color of the line showing the mean value.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#' Must be either .pdf, .png or .jpeg
#' 
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @export
plot_pval_hist <- function(dataset,
                           P,
                           bins = 20,
                           line_color = "black",
                           fill_color = "cornflowerblue",
                           mean_color = "red",
                           save_plot_path = FALSE,
                           plot_filename = "pvalue_histogram.png") {
  
  stopifnot("'P' must be a column in 'dataset'" = sub(".*\\$", "", deparse(substitute(P))) %in% colnames(dataset),
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
#' @description Plots a Manhattan plot with only the
#' SNPs that have a p-value of at most 0.1.
#' 
#' @details
#' All points above the upper dashed line show SNPs that are significant at
#' 0.05 significance level with Bonferroni-correction. Points above the lower
#' dashed line show SNPs that are significant at a 0.05 significance level.\cr
#' The color of the point show whether or not they are truly causal for the
#' phenotype.\cr
#' The user must run some analysis with \code{analysis_association()} 
#' and load the data with \code{load_results()}. Furthermore, \code{augment_results()}
#' must be run where the user should create a 'causal' column
#' from the true effect sizes. 
#'
#' @param dataset data imported to R by \code{load_results()}.
#' @param P column containing p-values to be plotted. State in the format
#' dataset$P_value.
#' @param SNP column containing the SNP values. If loaded with \code{load_results()}
#' can be left default otherwise state in the format dataset$SNP.
#' @param causal column containing the causal values. If created created with
#' \code{augment_results()} can be left default, otherwise state in the format
#' dataset$causal.
#' @param line_color color of the lines.
#' @param line_type the type of the lines.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#' Must be either .pdf, .png or .jpeg
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @export
plot_manhattan <- function(dataset,
                           P,
                           SNP = SNP,
                           causal = causal,
                           line_color = "cornflowerblue",
                           line_type = "dashed",
                           save_plot_path = FALSE,
                           plot_filename = "manhattan_plot.png") {
  
  stopifnot("'P', 'SNP' and 'causal' must be columns in 'dataset'" =
              all(c(sub(".*\\$", "", deparse(substitute(SNP))), 
                    sub(".*\\$", "", deparse(substitute(P))), 
                    sub(".*\\$", "", deparse(substitute(causal)))) 
                  %in% colnames(dataset)),
            "save_plot_path needs to be default or a valid path" =
              (save_plot_path == FALSE || dir.exists(save_plot_path)),
            "plot_filename must have either '.png', '.pdf' or '.jpeg' as extension" =
              (file_ext(plot_filename) == "png" ||
                 file_ext(plot_filename) == "pdf" ||
                 file_ext(plot_filename) == "jpeg"))
  
  data_plt <- dataset %>%
    dplyr::filter(P < 0.1)
  
  if(deparse(substitute(SNP)) != "SNP"){
    SNP <- SNP[P < 0.1]
  }
  
  if(deparse(substitute(causal)) != "causal"){
    causal <- causal[P < 0.1]
  }
  
  if(deparse(substitute(P)) != "P"){
    P <- P[P < 0.1]
  }
  
  m <- nrow(dataset)
  plt <- ggplot2::ggplot(data_plt) +
    ggplot2::geom_point(mapping = ggplot2::aes(SNP,
                                               -log10(P),
                                               colour = causal)) +
    ggplot2::geom_segment(ggplot2::aes(x = 1, xend = m, y = -log10(0.05),
                                       yend = -log10(0.05)),
                          colour = line_color,
                          linetype = line_type,
                          size = 1.25) +
    ggplot2::geom_segment(ggplot2::aes(x = 1,
                                       xend = m,
                                       y = -log10(0.05 / m),
                                       yend = -log10(0.05 / m)),
                          colour = line_color,
                          linetype = line_type, size = 1) +
    ggplot2::theme_light() +
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
#' @description Creates a scatterplot of true effect sizes against estimated effect sizes.
#' 
#' @details
#' In the plot an identity line is included, which helps visualise how well
#' the estimations compare to the true effect sizes. \cr
#' Only SNPs with a p-value smaller than 0.05 are included in the plot. \cr
#' The plot requires the user to run an analysis with
#' \code{analysis_association()}, load the results and merge with
#' \code{beta.txt} using \code{load_results()} and augment results with
#' \code{augment_results()}.
#'
#' @param dataset data imported to R by \code{load_results()}.
#' @param BETA column containing the estimated effect sizes from the analysis.
#' State in the format dataset$BETA.
#' @param P column containing the p-values from the analysis. 
#' State in the format dataset$P.
#' @param bonferroni the column containing significance at bonferroni correction.
#' State in the format dataset$bonferroni.
#' @param true_effect the column containing the true effect sizes. Default can be used
#' if loaded by \code{load_results()} together with \code{beta.txt}, else
#' state in the format dataset$true_effect.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#' Must be either .pdf, .png or .jpeg
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @export
plot_estimates_vs_true <- function(dataset,
                                   BETA,
                                   P,
                                   bonferroni,
                                   true_effect = beta,
                                   save_plot_path = FALSE,
                                   plot_filename = "beta_comparison.png") {
  
  stopifnot("'BETA', 'P', 'true_effect' and 'bonferroni' must be columns in 'dataset'" =
          all(c(sub(".*\\$", "", deparse(substitute(BETA))), 
                    sub(".*\\$", "", deparse(substitute(true_effect))), 
                    sub(".*\\$", "", deparse(substitute(bonferroni))),
                    sub(".*\\$", "", deparse(substitute(P)))) 
                  %in% colnames(dataset)),
    "save_plot_path needs to be default or a valid path" =
              (save_plot_path == FALSE || dir.exists(save_plot_path)),
            "plot_filename must have either '.png', '.pdf' or '.jpeg' as extension" =
              (file_ext(plot_filename) == "png" ||
                 file_ext(plot_filename) == "pdf" ||
                 file_ext(plot_filename) == "jpeg"))
  
  tmpdataacc <- dataset %>%
    dplyr::filter(P < 0.05)
  
  if(deparse(substitute(BETA)) != "BETA"){
    BETA <- BETA[P < 0.05]
  }
  
  if(deparse(substitute(true_effect)) != "true_effect"){
    true_effect <- true_effect[P < 0.05]
  }
  
  if(deparse(substitute(bonferroni)) != "bonferroni"){
    bonferroni <- bonferroni[P < 0.05]
  }
  
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
#' @description Creates a scatterplot of genetic liability against posterior mean genetic liability.
#' 
#' @details
#' Visualises how well the posterior mean genetic liabilities fit the genetic liabilities.
#' An identity line is included to help see if the genetic liabilities are balanced above
#' and below for each posterior mean genetic liability. Furthermore, the color and labels
#' show each configuration class, which is detailed in
#' \code{vignette("liability-distribution")}. \cr
#' This function requires the user to run \code{assign_ltfh_phenotype()}, and then
#' \code{association_analysis(pheno_name = "LTFH_pheno")}. Data is loaded with
#' \code{load_phenotypes()}. \cr
#' \code{vignette("ggplot2-specs", package = "ggplot2")} details the usage of
#' aesthetic parameters in ggplot2.
#'
#' @param dataset data imported to R by \code{load_phenotypes()}.
#' @param line_color color of identity line.
#' @param line_type type of identity line.
#' @param label_size the size of the boxes.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#' Must be either .pdf, .png or .jpeg
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @export
plot_pmgl_vs_true <- function(dataset,
                              line_color = "black",
                              line_type = "dashed",
                              label_size = 3,
                              save_plot_path = FALSE,
                              plot_filename = "posterior_liabilities.png") {
  
  stopifnot("dataset must have a columns named 'LTFH_pheno', 'child_lg' and 'conf_class'" =
              all(c("LTFH_pheno", "child_lg", "conf_class") %in% colnames(dataset)),
            "save_plot_path needs to be default or a valid path" =
              (save_plot_path == FALSE || dir.exists(save_plot_path)),
            "plot_filename must have either '.png', '.pdf' or '.jpeg' as extension" =
              (file_ext(plot_filename) == "png" ||
                 file_ext(plot_filename) == "pdf" ||
                 file_ext(plot_filename) == "jpeg"))
  
  repel_data <- dataset[, .SD[which.min(child_lg)], by = conf_class]
  
  plt <- ggplot2::ggplot(dataset) +
    ggplot2::geom_point(ggplot2::aes(LTFH_pheno, child_lg, color = conf_class),
                        size = 0.6,
                        show.legend = FALSE) +
    ggplot2::geom_abline(color = line_color, linetype = line_type) +
    ggplot2::theme_light() +
    ggplot2::labs(x = "Posterior mean genetic liability",
                  y = "True genetic liability",
                  title = "Scatterplot of liabilities and LT-FH phenotype",
                  subtitle = deparse(substitute(dataset))) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(face = "bold",
                                                         hjust = 0.5)) +
    ggrepel::geom_label_repel(data = repel_data,
                              ggplot2::aes(LTFH_pheno, child_lg, 
                                           label = gsub("x", ".", conf_class),
                                           color = as.factor(conf_class)),
                              nudge_y = -3 - repel_data$child_lg,
                              box.padding = 0.4,
                              label.padding = 0.1,
                              show.legend = FALSE,
                              size = label_size)
  
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
#' @description Creates a scatterplot of transformed linear GWAS beta values against the LTFH_transformed values.
#' 
#' @details 
#' In the plot an identity line is included, which helps visualise how well
#' the estimations compare to each other. \cr
#' The plot requires the user to run an analysis with
#' \code{analysis_association()}, load the results and merge with
#' \code{beta.txt} using \code{load_results()} and augment results with
#' \code{augment_results()} choosing to create the LTFH transformed column. \cr
#' \code{vignette("ggplot2-specs", package = "ggplot2")} details the usage of
#' aesthetic parameters in ggplot2.
#' 
#'
#' @param dataset data imported to R by \code{load_results()}.
#' @param beta_x column containing beta-values to be displayed on the x-axis.
#' E.g. the transformed LTFH beta values. State in the format dataset$beta_x.
#' @param beta_y column containing beta-values to be displayed on the y-axis.
#' E.g. the linear GWAS beta values. State in the format dataset$beta_y.
#' @param bonferroni column containing the corrected significant values from LTFH.
#' State in the format dataset$bonferroni.
#' @param line_color the color of the identity line.
#' @param line_type the type of the identity line.
#' @param save_plot_path if \code{FALSE}, return the function returns a ggplot
#' object. Else a path of the directory to save the plot to.
#' @param plot_filename name of the file to be saved, including file extension.
#' Must be either .pdf, .png or .jpeg
#'
#' @return Either returns a \code{ggplot} object or saves the plot to
#' \code{save_plot_path} and returns NULL.
#'
#' @export
compare_beta <- function(dataset, beta_x, beta_y, bonferroni,
                         line_color = "green", line_type = "dashed",
                         save_plot_path = FALSE, 
                         plot_filename = "Compare_beta.png") {
  
  stopifnot("save_plot_path needs to be default or a valid path" =
              (save_plot_path == FALSE || dir.exists(save_plot_path)),
            "plot_filename must have either '.png', '.pdf' or '.jpeg' as extension" =
              (tools::file_ext(plot_filename) == "png" ||
                 tools::file_ext(plot_filename) == "pdf" ||
                 tools::file_ext(plot_filename) == "jpeg"))
  
  plt <- ggplot2::ggplot(dataset, ggplot2::aes(x = beta_x,
                                               y = beta_y)) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = bonferroni)) +
    ggplot2::geom_smooth(se = FALSE, method = "lm")+
    ggplot2::geom_abline(intercept = 0, slope = 1, color=line_color, 
                         linetype=line_type, size = 1) +
    ggplot2::labs(title = "Comparing beta values",
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

