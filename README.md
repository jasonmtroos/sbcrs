
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sbcrs

<!-- badges: start -->

<!-- badges: end -->

This package provides tools to simplify the implementation of simulation
based calibration using rank statistics (Talts, Betancourt, Simpson,
Vehtari, and Gelman,
[arXiv:1804.06788](https://arxiv.org/abs/1804.06788)). It implements a
very similar validation procedure to that in `rstan::sbc` but using a
different workflow. In `rstan::sbc`, the Stan model must be modified by
the user to generate rank statistics during sampling. In this package,
the Stan model is left unmodified, and the code needed to calculate the
rank statistics is written in R. This provides a potentially faster
development workflow (Stan recompiles are not needed). It also allows
SBC to be used in cases where generating parameters from their prior
distributions, or the modeled variable from the likelihood function,
would be too complicated if written in
Stan.

## Installation

<!-- You can install the released version of sbcrs from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("sbcrs") -->

<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jasonmtroos/sbcrs")
```

To build the package vignettes, install the package using:

``` r
devtools::install_github("jasonmtroos/sbcrs", build_vignettes = TRUE)
```

The package vignettes are a useful starting point.

  - `intro-to-sbc` provides an overview to simulation-based calibration,
    and the features of the SBC package

  - `funnel` shows how SBC identified sampling problems in Neal’s funnel

  - `comparison-to-rstan-sbc` implements the example from `rstan::sbc`
    using this package, shows the rank statistics are the same, and
    provides a basis for understanding the different design philosophies
    behind the two approaches.
