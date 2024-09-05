## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages(require(TCHazaRds))   # this package :)
suppressPackageStartupMessages(require(terra))       # spatial analysis
suppressPackageStartupMessages(require(rasterVis))   # enhanced raster visualization https://oscarperpinan.github.io/rastervis/
suppressPackageStartupMessages(require(sp))          # spatial methods and plotting
suppressPackageStartupMessages(require(knitr))       # formatted table
suppressPackageStartupMessages(require(raster))       # convert for raster plots


## -----------------------------------------------------------------------------
TCi = vect(cbind(c(154,154),c(-26.1,-26)),"lines",crs="epsg:4283") #track line segment
TCi$PRES = 950 #central pressure in hPa
#TCi$RMW = 40 #radius of maximum winds in km
TCi$ISO_TIME = "2022-10-04 20:00:00" #"%Y-%m-%d %H:%M:%S", tz = "UTC"
TCi$LON = geom(TCi)[1,3] #longitude
TCi$LAT = geom(TCi)[1,4] #latitude
TCi$STORM_SPD = perim(TCi)/(3*3600) #speed of the forward motion of the TC m/s
TCi$thetaFm = 90-returnBearing(TCi) #direction of the heading of the TC (Cartesian, clockwise from x axis)


## -----------------------------------------------------------------------------
TC <- vect(system.file("extdata/YASI/YASI.shp", package="TCHazaRds"))
TC$PRES <- TC$BOM_PRES #different agencies each provide a PRES, you need to chose one. 
TC$STORM_SPD = TC$STORM_SPD/1.94 #provided as knots, convert to m/s
TC$thetaFm = 90-returnBearing(TC) #direction of the heading of the TC (Cartesian, clockwise from x axis)
TCi = TC[46]

## -----------------------------------------------------------------------------
paramsTable = read.csv(system.file("extdata/tuningParams/defult_params.csv",package = "TCHazaRds"))
#paramsTable$value[6:7] = c(0,0)
paramsTable#knitr::kable(paramsTable,caption = "Parameter file")

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
r = rast(xmin = 145,xmax=149,ymin = -19,ymax = -16.5,resolution=.01)
values(r) = 0
#GEO_land = land_geometry(r,r)

#
land_v <- vect(system.file("extdata/OSM_500m_QLD/OSM_500m_QLD.shp", package="TCHazaRds"))
land_r = rasterize(land_v,r,touches=TRUE,background=0)
inland_proximity = terra::costDist(land_r,target = 0,scale=1)
GEO_land = land_geometry(land_r,inland_proximity)

#plot(inland_proximity,main = "Inland Distance (m)")
#plot(TC,add=TRUE)


