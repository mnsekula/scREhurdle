# scREhurdle

**scREhurdle** is an R package for detecting differntially expressed genes in discrete single-cell RNA sequencing data. This package interfaces with [`rstan`](https://mc-stan.org/users/interfaces/rstan) and fits a mixed effect hurdle model on zero-inflated count data. The `rstan` package should be properly
installed before installing **scREhurdle**. See [Rstan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for details.

## Installation
To install this package from GitHub use:
```{r, warning=FALSE}
# Install from GitHub using devtools
require(devtools)
devtools::install_github("mnsekula/scREhurdle")
```

Note: There will likely be `rstan` compiler warnings related to certain dependent packages (e.g., RcppEigen, StanHeaders). These compiler warnings can safely be ignored as long as the **scREhurdle** package is installed successfully. The code below should output `TRUE`:

```{r, warning=FALSE}
# Check if scREhurdle was successfully installed
is.element("scREhurdle", installed.packages()[,1])
```

For more information on Stan's warnings, please visit [Brief Guide to Stan’s Warnings](https://mc-stan.org/misc/warnings.html).

## Getting Started
The package vignette demonstrates how to use the **scREhurdle** package to perform a differential expression analysis. This vignette can be viewed online [here](http://htmlpreview.github.io/?https://github.com/mnsekula/scREhurdle/blob/master/scREhurdle.html).


