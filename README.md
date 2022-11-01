# BES Palaeo R Workshop GAMs

### Slides

* [HTML slidedeck](https://gavinsimpson.github.io/bes-paleo-r-workshop/slides.html)

## Installation instructions for the R part

It is essential that you are running a recent version of R. {gratia} needs at least R 4.1.0. The version of {gratia} on CRAN works in earlier versions of R, but the changes in the about-to-be-released version of {gratia} need 4.1.0 or later.

Please install the following packages. You may already have most of them from Steve's session on Day 1:

```r
pkgs <- c("readr", "ggplot2", "mgcv", "dplyr", "readxl", "patchwork",
          "tibble", "tidyr", "rio", "purrr", "vegan")
install.packages(pkgs, Ncpus = 3) # change to 2, 4 etc if you have more cores
```

Finally, you need to install the development version of {gratia}. You have two options for installing this:

```r
# option 1
options(repos = c(
  gavinsimpson = "https://gavinsimpson.r-universe.dev",
  CRAN = "https://cloud.r-project.org"))
# Download and install gratia in R
install.packages("gratia")

# or option 2
remotes::install_github("gavinsimpson/gratia")
```

Use Option 1 preferably, especially if you had problems install {riojaPlot} using `remotes::install_github()` on Day 1.