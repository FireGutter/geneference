% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisations.R
\name{plot_estimates_vs_true}
\alias{plot_estimates_vs_true}
\title{Plot true effects against estimated effects of SNPs}
\usage{
plot_estimates_vs_true(
  dataset,
  BETA,
  P,
  bonferroni,
  true_effect,
  save_plot_path = FALSE,
  plot_filename = "beta_comparison.png"
)
}
\arguments{
\item{dataset}{data imported to R by \code{load_results()}.}

\item{BETA}{column containing the estimated effect sizes from the analysis.
State in the format dataset$BETA.}

\item{P}{column containing the p-values from the analysis. 
State in the format dataset$P.}

\item{bonferroni}{the column containing significance at bonferroni correction.
State in the format dataset$bonferroni.}

\item{true_effect}{the column containing the true effect sizes.
State in the format dataset$true_effect. Note: standard called 'beta' if
loaded by \code{load_results()}.}

\item{save_plot_path}{if \code{FALSE}, return the function returns a ggplot
object. Else a path of the directory to save the plot to.}

\item{plot_filename}{name of the file to be saved, including file extension.
Must be either .pdf, .png or .jpeg}
}
\value{
Either returns a \code{ggplot} object or saves the plot to
\code{save_plot_path} and returns NULL.
}
\description{
Creates a scatterplot of true effect sizes against estimated effect sizes.
}
\details{
In the plot an identity line is included, which helps visualise how well
the estimations compare to the true effect sizes. \cr
Only SNPs with a p-value smaller than 0.05 are included in the plot. \cr
The plot requires the user to run an analysis with
\code{analysis_association()}, load the results and merge with
\code{beta.txt} using \code{load_results()} and augment results with
\code{augment_results()}.
}
