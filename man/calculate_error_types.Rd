% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/error_types.R
\name{calculate_error_types}
\alias{calculate_error_types}
\title{Calculate error types}
\usage{
calculate_error_types(causal, significant, bonferroni, analysis_name)
}
\arguments{
\item{causal}{column in dataset denoting causal SNPs.}

\item{significant}{column in dataset denoting which SNPs are significant.}

\item{bonferroni}{column in dataset denoting which SNPs are significant with 
Bonferroni correction.}

\item{analysis_name}{name of the analysis.}
}
\value{
A tibble containing the number of true positives, false positives,
false negatives and true negatives with each of the significance levels
having their own row.
}
\description{
Calculate the number of type 1 (false positive) and type 2 (false negative)
errors in an analysis.
}
\details{
To use this function, the user must run any test with \code{analysis_association()}
and use \code{load_results()} to merge this with \code{beta.txt} from the
simulation. Then, \code{augment_results()} needs to be run in order to obtain
columns for information on which SNPs are causal, significant, and significant
with Bonferroni correction. \cr
The standard significance level refers to the one used in 
\code{augment_results()}.
}
