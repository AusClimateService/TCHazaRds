# TCHazaRds 

[![CRAN
status](https://www.r-pkg.org/badges/version/TCHazaRds)](https://cran.r-project.org/package=TCHazaRds)

<img align="right" width="250" src="TCHazaRd_logo.png">

`TCHazaRds` is an *R* package for Tropical Cyclone (Hurricane, Typhoon) Spatial Hazard Modelling. There is a tutorial in [TCHazaRds/vignettes](https://htmlpreview.github.io/?https://github.com/AusClimateService/TCHazaRds/blob/main/vignettes/Introduction_to_TCHazaRds.html). The code includes (but is not limited to) parametrisations such as wind reduction factors [Harper et. al. 2001](https://data.longpaddock.qld.gov.au/static/publications/vulnerability-to-cyclones/stage1.pdf) to estimate surface mean 10 min wind speeds from Cyclone tracks. Further modelling details can be found in [O'Grady et. al. 2024](https://journals.ametsoc.org/view/journals/mwre/152/1/MWR-D-23-0063.1.xml).

## Installation

### For Windows
`TCHazaRds` is now on [CRAN](https://cran.r-project.org/package=TCHazaRds) and can downloaded with

```
install.packages("TCHazaRds",dependencies = TRUE)
```

### From source-code

For windows you need to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to get a C++ compiler that R can use. You need a recent version of Rtools42 (rtools42-5355-5357).

Then, in R, install the package.

```
require(remotes)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github("AusClimateService/TCHazaRds", dependencies = TRUE)
```
