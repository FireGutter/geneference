
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geneference <img src="man/figures/logo.png" align="right" width="120"/>

#### Causal inference with simulated genetic data

<!-- badges: start -->

[![R-CMD-check](https://github.com/FireGutter/geneference/workflows/R-CMD-check/badge.svg)](https://github.com/FireGutter/geneference/actions)
[![CodeFactor](https://www.codefactor.io/repository/github/firegutter/geneference/badge?s=d63fc08844421c6003bc29b6e159ab5051f1f619)](https://www.codefactor.io/repository/github/firegutter/geneference)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT/)
<!-- badges: end -->

geneference is an R package for simulating and analysing genetic data
for genome wide association studies.  
The implementations are primarily based on the work by Hujoel et al
(2020), in which the authors propose a method of including family
history in a representation of a case-control phenotype. The authors
show that their method, named Liability Threshold conditional on Family
History, or LT-FH in short, leads to a significant increase in
association power in genome wide association studies.  
geneference aims to provide a basis for investigating the merits of
LT-FH compared to traditional approaches in genome wide association
studies. With geneference, the user can control several parameters in
the simulations, carry out statistical analysis with different
association methods and use plotting functions to visualise certain
aspects of the analyses.  
geneference can be used as a stand-alone package, but those interested
in diving deeper and/or using geneference as a starting point for their
own implementations can freely browse our source code at our [repository
on Github](https://github.com/FireGutter/geneference/).

geneference was developed as an exam project for the Data Project
course, which is part of the Data Science bachelor programme at Aarhus
University.

## Usage

You can install the development version from
[GitHub](https://github.com/) with:

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("FireGutter/geneference")
```

Some of the functions implemented in geneference rely on the software
PLINK, which the user will need to install. However, geneference handles
the calls to PLINK, meaning that besides installation the user will not
need to interact directly with PLINK.

If you are interested in seeing what geneference is capable of, you can
read through the introductory vignette, `vignette("geneference")`. For a
complete overview of the functions in geneference, see our [package
documentation](https://firegutter.github.io/geneference/reference/index.html),
and for further details on the theoretical background of our
implementations, see [our
vignettes](https://firegutter.github.io/geneference/articles/index.html).

## References

-   Hujoel, M.L.A., Gazal, S., Loh, PR. et al. Liability threshold
    modeling of case–control status and family history of disease
    increases association power. Nat Genet 52, 541–547 (2020).
    <https://doi.org/10.1038/s41588-020-0613-6>

Additionally, the open source software Inkscape was used together with
the hexSticker package to create the logo for geneference.  
Finally, we want to thank our supervisor Emil Michael Pedersen for his
patience in answering all our questions.
