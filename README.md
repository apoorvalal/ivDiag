
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ivDiag

<!-- badges: start -->
<!--
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stablel)
-->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**`ivDiag`** is toolkit for estimation and diagnostics with instrumental
variable (IV) designs. It provides `R` implementations of the guidelines
proposed in Lal et al. (2023), enabling researchers to obtain reliable
and robust estimates of their IV models.

**Examples:** `R` code used in the
[tutorial](https://yiqingxu.org/packages/ivDiag/articles/iv_tutorial.html)
can be downloaded from
[here](https://raw.githubusercontent.com/xuyiqing/ivDiag/main/pkgdown/ivDiag_examples.R).

**Reference:** : Lal, Apoorva, Mackenzie William Lockhart, Yiqing Xu,
and Ziwen Zu (2023). [How Much Should We Trust Instrumental Variable
Estimates in Political Science? Practical Advice Based on 67 Replicated
Studies](https://yiqingxu.org/papers/english/2021_iv/LLXZ.pdf), Mimeo,
Stanford University.

------------------------------------------------------------------------

## Installation

You can install the **ivDiag** package from CRAN:

``` r
install.packages('ivDiag') # not on CRAN yet 
```

You can also install the up-to-date development version from Github:

``` r
library(remotes)
install_github("apoorvalal/ivDiag")
```

**ivDiag** depends on the following packages, which will be installed
automatically when **ivDiag** is being installed; you can also install
them manually:

``` r
require(foreach) 
require(future)
require(lfe)
require(fixest)  
require(ggplot2)
require(ggfortify)
require(wCorr)
require(haven)
require(ggdistribution)
require(glue)
require(patchwork)
```

------------------------------------------------------------------------

## Report bugs

This package is a work in progress.
[Github](https://github.com/apoorvalal/ivDiag) issues and/or pull
requests are welcome.
