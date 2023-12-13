# TCHazaRds 

<img align="right" width="250" src="TCHazaRd_logo.png">

`TCHazaRds` is an *R* package for Tropical Cyclone (Hurricane, Typhoon) Spatial Hazard Modelling. There is a tutorial in [TCHazaRds/vignettes](https://htmlpreview.github.io/?https://github.com/AusClimateService/TCHazaRds/blob/main/vignettes/Introduction_to_TCHazaRds.html). 

## Installation

For windows you need to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to get a C++ compiler that R can use. You need a recent version of Rtools42 (rtools42-5355-5357).

Then, in R, install the package.

```
require(remotes)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github("AusClimateService/TCHazaRds", dependencies = TRUE)
```
