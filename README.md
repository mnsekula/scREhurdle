# scREhurdle

scREhurdle is an R package for detecting differntially expressed genes in discrete single-cell RNA sequencing data. This package is built on top of [`rstan`](https://mc-stan.org/users/interfaces/rstan) and fits a mixed effect hurdle model on zero-inflated count data.

## Installation
To install this package from GitHub use:
```{r, warning=FALSE}
# Install from GitHub using devtools
require(devtools)
devtools::install_github("mnsekula/scREhurdle")
```

## Getting Started
The vignette demonstrates how to use the scREhurdle package to perform a differential expression analysis. Once installed, the vignette can be obtained with:

```{r, warning=FALSE}
library(scREhurdle)
vignette("scREhurdle")
```
