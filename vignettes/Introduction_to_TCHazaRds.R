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
TCi = TC[47]

## -----------------------------------------------------------------------------
paramsTable = read.csv(system.file("extdata/tuningParams/defult_params.csv",package = "TCHazaRds"))
#paramsTable$value[6:7] = c(0,0)
knitr::kable(paramsTable,caption = "Parameter file")

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
r = rast(xmin = 145,xmax=149,ymin = -19,ymax = -16.5,resolution=.01)
values(r) = 0
GEO_land = land_geometry(r,r)

#
#land_v <- vect(system.file("extdata/OSM_500m_QLD/OSM_500m_QLD.shp", package="TCHazaRds"))
#land_r = rasterize(land_v,r,touches=TRUE,background=0)
#inland_proximity = costDistance(land_r,target = 0,scale=1)
#GEO_land = land_geometry(land_r,inland_proximity)

#plot(inland_proximity,main = "Inland Distance (m)")
#plot(TC,add=TRUE)


## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
ats = seq(0, 65, length=14)
HAZi = TCHazaRdsWindField(GEO_land = GEO_land,TC = TCi,paramsTable=paramsTable)
library(raster)       # convert for raster plots
dummy = raster::raster() 
TC_sp = list("sp.lines",as(TC,"Spatial"),col="black")
sp::spplot(HAZi,"Sw",at=ats,sp.layout = TC_sp)

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
if (.Platform$OS.type == "windows"){
  UV = as(c(HAZi["Uw"],HAZi["Vw"]),"Raster") #need to convert back to raster
  rasterVis::vectorplot(UV, isField='dXY', col.arrows='white', aspX=0.002,aspY=0.002,at=ats ,
  colorkey=list(at=ats), par.settings=viridisTheme)+latticeExtra::layer(sp.lines(as(TC,"Spatial"),col="red"))
}

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
HAZ = TCHazaRdsWindFields(GEO_land=GEO_land,TC=TC,paramsTable=paramsTable)
sp::spplot(max(HAZ$Sw),at=ats,sp.layout = TC_sp)

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
outdate = seq(strptime(TC$ISO_TIME[1],"%Y-%m-%d %H:%M:%S",tz="UTC"),
              strptime(rev(TC$ISO_TIME)[1],"%Y-%m-%d %H:%M:%S",tz="UTC"),
              3600)
HAZI = TCHazaRdsWindFields(outdate=outdate,GEO_land=GEO_land,TC=TC,paramsTable=paramsTable)
sp::spplot(max(HAZI$Sw),at=ats,sp.layout = TC_sp)

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------

outdate = seq(strptime(TC$ISO_TIME[1],"%Y-%m-%d %H:%M:%S",tz="UTC"),
              strptime(rev(TC$ISO_TIME)[1],"%Y-%m-%d %H:%M:%S",tz="UTC"),
              600)

GEO_landp = data.frame(dem=0,lons = 147,lats=-18,f=-4e-4,inlandD = 0)
HAZts = TCHazaRdsWindTimeSereies(GEO_land=GEO_landp,TC=TC,paramsTable = paramsTable)
HAZtsi = TCHazaRdsWindTimeSereies(outdate = outdate,GEO_land=GEO_landp,TC=TC,paramsTable = paramsTable)

main =  paste(TCi$NAME[1],TCi$SEASON[1],"at",GEO_landp$lons,GEO_landp$lats)
if (.Platform$OS.type == "windows"){ 
  suppressWarnings(with(HAZts,plot(date,Sw,format = "%b-%d %HZ",type="l",main = main,ylab = "Wind speed [m/s]")))
  with(HAZtsi,lines(date,Sw,col=2))
  legend("topleft",c("6 hrly","10 min interpolated"),col = c(1,2),lty=1)
}

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
TCi$thetaFm = 90-returnBearing(TCi)
pp <- TCProfilePts(TC_line = TCi,bear=TCi$thetaFm+90,length =150,step=1)
#extract the GEO_land
GEO_land_v = extract(GEO_land,pp,bind=TRUE,method = "bilinear")
HAZp = TCHazaRdsWindProfile(GEO_land_v,TCi,paramsTable)

HAZie = extract(HAZi,pp,bind=TRUE)#,method = "bilinear")

wcol = colorRampPalette(c("white","lightblue","blue","violet","purple"))
#see ?terra::plot
plot(HAZi,"Sw",levels=ats,col = wcol(13),range = range(ats),type="continuous",all_levels=TRUE)
#plot(HAZp,add=TRUE,cex=1.2)
plot(HAZp,"Sw",levels=ats,col = wcol(13),range = range(ats),type="continuous",border="grey")#,all_levels=TRUE)
lines(TC)

## ----out.width = '80%',fig.height=4,fig.width=6, fig.align = "center"---------
plot(HAZp$radialdist,HAZp$Sw,type="l",xlab = "Radial distance [km]",ylab = "Wind speed [m/s]");grid()


