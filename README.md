# scREhurdle

**scREhurdle** is an R package for detecting differentially expressed genes in discrete single-cell RNA sequencing data. This package interfaces with [`rstan`](https://mc-stan.org/users/interfaces/rstan) and fits a mixed effect hurdle model on zero-inflated count data. The `rstan` package should be
installed and working properly before installing **scREhurdle**. Note: the optional step of configuring the C++ toolchain in the `rstan` installation instructions is required for compiling the code in **scREhurdle**. See [Rstan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for more details on the `rstan` installation process.

### Windows Users
Please verify that `CXX14 = g++ -m$(WIN) -std=c++1y` is included in the Makevars file. This file is created during the optional step of configuring the C++ toolchain in the `rstan` installation process. The following code will open the Makevars file from R.

```{r, warning=FALSE}
M <- file.path(Sys.getenv("HOME"), ".R", ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
file.edit(M)
```

If this line is not present, manually add `CXX14 = g++ -m$(WIN) -std=c++1y` directly to this file and save.


Also, make sure that Rtools and MinGW are defined in `Sys.getenv("PATH")` and that `Sys.getenv("BINPREF")` is set to MinGW. The code below will allow you to check for and add these settings to your current R session.

```{r, warning=FALSE}
# Check PATH for Rtools and mingw
Sys.getenv("PATH")

# Check for BINPREF
Sys.getenv("BINPREF")

# Add Rtools and mingw to the current session PATH
Sys.setenv(PATH = paste("C:/Rtools/mingw_64/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))

# Set BINPREF to mingw
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
```

Note: `C:/` in the code above is a default location of Rtools and may need to be changed to reflect the location of Rtools on your machine.


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

