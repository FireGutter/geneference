% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_map.R
\name{create_map}
\alias{create_map}
\title{Create .map-file for PLINK}
\usage{
create_map(n, path)
}
\arguments{
\item{n}{number of sequenced genotypes.}

\item{path}{location where the file is going to be stored.}
}
\value{
Does not return anything, but creates a file with the name
"genotypes.map" in the folder specified with 'path'-parameter.
}
\description{
As PLINK requires a .map-file in order to run its analysis,
we create such a file and fill it with some arbitrary information, since
it is not actually used for the analysis.
}
