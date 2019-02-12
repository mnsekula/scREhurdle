# scREhurdle

**scREhurdle** is an R package for detecting differentially expressed genes in discrete single-cell RNA sequencing data. This package interfaces with [`rstan`](https://mc-stan.org/users/interfaces/rstan) and fits a mixed effect hurdle model on zero-inflated count data. The `rstan` package should be
installed and working properly before installing **scREhurdle**. See [Rstan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for details.

## Installation
For most users, this package can be installed from GitHub with:
```{r, warning=FALSE}
# Install from GitHub using devtools
require(devtools)
devtools::install_github("mnsekula/scREhurdle")
```
Note: There may be some `rstan` compiler warnings related to certain dependent packages (e.g., RcppEigen, StanHeaders, etc.). According to the [Brief Guide to Stan’s Warnings](https://mc-stan.org/misc/warnings.html), these warnings can safely be ignored as long as the **scREhurdle** package is installed successfully. The code below should output `TRUE` for a successful installation:

```{r, warning=FALSE}
# Check if scREhurdle was successfully installed
is.element("scREhurdle", installed.packages()[,1])
```

Windows users with both 32 and 64 bit versions of R installed on their machine may need to install this package from GitHub with:
```{r, warning=FALSE}
# Install only on the main architecture from GitHub using devtools 
require(devtools)
devtools::install_github("mnsekula/scREhurdle", INSTALL_opts = "--no-multiarch")
```
This will build and install the **scREhurdle** package for the version of R (either 32 bit or 64 bit) the code is run on.

## Getting Started
The package vignette demonstrates how to use the **scREhurdle** package to perform a differential expression analysis. This vignette can be viewed online [here](http://htmlpreview.github.io/?https://github.com/mnsekula/scREhurdle/blob/master/scREhurdle.html).

