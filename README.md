# scREhurdle

scREhurdle is an R package for detecting differntially expressed genes in discrete single-cell RNA sequencing data. This package is built on top of [`rstan`](https://mc-stan.org/users/interfaces/rstan) and fits a mixed effect hurdle model on zero-inflated count data. The `rstan` package should be properly
installed before installing scREhurdle. See [Rstan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for details.

## Installation
To install this package from GitHub use:
```{r, warning=FALSE}
# Install from GitHub using devtools
require(devtools)
devtools::install_github("mnsekula/scREhurdle", build_vignettes = TRUE)
```
Setting `build_vignettes = FALSE` will build the package without the vignette.

## Getting Started
The vignette demonstrates how to use the scREhurdle package to perform a differential expression analysis and can be viewed online [here](http://htmlpreview.github.io/?https://github.com/mnsekula/scREhurdle/blob/master/scREhurdle.html).

If the vignette was built when installing scREhurdle from GitHub, it can also be obtained locally in R with:

```{r, warning=FALSE}
library(scREhurdle)
vignette("scREhurdle")
```

