#' Update Parameter List to Calibrated Values
#'
#' @param paramsTable Global parameters to compute TC Hazards.
#' @param infile File containing tuning parameters in a .csv. Default for QLD calibration.
#'
#' @return list of params with updated tuning wind parameters.
#' @export
#'
#' @examples
#' paramsTable <- read.csv(system.file("extdata/tuningParams/defult_params.csv",package = "TCHazaRds"))
#'
#' tunedParams(paramsTable)
tunedParams = function(paramsTable,infile = system.file("extdata/tuningParams/QLD_modelSummaryTable.csv",package = "TCHazaRds")){
  params <- array(paramsTable$value,dim = c(1,length(paramsTable$value)))
  colnames(params) <- paramsTable$param
  params <- data.frame(params)
  if(with(params,rMaxModel != vMaxModel | rMaxModel !=betaModel)) stop("Along track models (beta, rMax, vMax) must have the same value")
  tunedvalues <- utils::read.csv(infile)
  tunedvalues <- tunedvalues[tunedvalues$simulation == "DecayCorrect",]
  vmods = c("Kepert","Hubbert","McConochie")
  tunedvalues = tunedvalues[tunedvalues$Dataset == vmods[params$windVortexModel+1] & tunedvalues$param_mod == params$rMaxModel,]
  paramsTable[which(paramsTable[,1] == "Decay_a1"),]$value <- tunedvalues$a1
  paramsTable[which(paramsTable[,1] == "Decay_a2"),]$value <- tunedvalues$a2
  paramsTable[which(paramsTable[,1] == "Decay_a3"),]$value <- tunedvalues$a3
  return(paramsTable)
}

#' @title Calculate the Geometric Parameters for Terrestrial Wind
#' @description Returns geometric data to compute wind fields.
#' @param dem SpatRaster  object, digital elevation model
#' @param inland_proximity SpatRaster  object, distance from the coast inland
#' @param returnpoints Return SpatVector of points or SpatRaster
#' @return SpatVector with attributes or SpatRaster
#'
#' | Abbreviated attribute       | description  | units |
#' | ------------- | ------------- | ------------- |
#' | dem      | Digital Elevation Model | m |
#' | lat      | Latitude  | degs  |
#' | lon      | Longitude | degs    |
#' | slope      | slope of terrain | - |
#' | aspect      | DEM aspect | - |
#' | inlandD      | distance inland from coast | m |
#' | f        | Coriolis parameter | hz |
#'
#' @md
#' @export
#'
#' @examples
#' require(terra)
#' dem <- rast(system.file("extdata/DEMs/YASI_dem.tif", package="TCHazaRds"))
#' land <- dem; land[land > 0] = 0
#' inland_proximity = distance(land,target = 0)
#' GEO_land = land_geometry(dem,inland_proximity)
#' plot(GEO_land)
land_geometry = function(dem,inland_proximity,returnpoints=FALSE){
    dem[dem < 0] = NA #don't include underwater slopes
    land_slope = terra::terrain(dem,"slope",neighbors = 4)
    land_aspect = terra::terrain(dem,"aspect",neighbors = 4)
    msk = dem
    msk=msk/msk

    rex = terra::ext(dem)
    rre = terra::res(dem)
    latlon = expand.grid(seq(rex[1]+rre[1]/2,rex[2]-rre[1]/2,rre[1]),rev(seq(rex[3]+rre[2]/2,rex[4]-rre[2]/2,rre[2])))

    lons = dem*0;terra::values(lons) <- latlon[,1]
    lats = dem*0;terra::values(lats) <- latlon[,2]
    wearth <- pi*(1.0 / 24.0) / 1800.0
    f = 2*wearth*sin(lats*pi / 180.0)

    if(returnpoints){
      GEO = terra::as.points(dem*msk, na.rm= FALSE)
      names(GEO) = "dem"
      g = terra::geom(GEO)
      GEO$lats = g[,4] #return latitude
      GEO$lons = g[,3] #return longitude
      GEO$slope = terra::extract(land_slope*msk,GEO)[,2]
      GEO$aspect = terra::extract(land_aspect*msk,GEO)[,2]
      GEO$inlandD = terra::extract(inland_proximity*msk,GEO)[,2]
    }
    if(!returnpoints){
      GEO = terra::rast(list(dem,lons,lats,land_slope*msk,land_aspect*msk,inland_proximity*msk,f))
      names(GEO) = c("dem","lons","lats","slope","aspect","inlandD","f")
    }
    return(GEO)
}



#' Calculate Additional TC Parameters, and temporally Interpolate Along a Tropical Cyclone Track
#'
#' @param outdate POSIX times to be interpolated to
#' @param indate POSIX input times
#' @param TClons input central TC longitude
#' @param TClats input central TC latitude
#' @param vFms input forward velocity of TC
#' @param thetaFms input forward direction
#' @param cPs central pressure
#' @param rMaxModel empirical model for radius of maximum wind calculation (rMax in km)
#' @param vMaxModel empirical model for maximum wind velocity calculation (vMax in m/s)
#' @param betaModel empirical model for TC shape parameter beta (dimensionless Beta)
#' @param eP background environmental pressure (hPa)
#' @param rho air density
#'
#' @return list of track data inclining the rMax vMax and Beta.
#' @export
#'
#' @examples
#' paramsTable <- read.csv(system.file("extdata/tuningParams/defult_params.csv",package = "TCHazaRds"))
#' params <- array(paramsTable$value,dim = c(1,length(paramsTable$value)))
#' colnames(params) <- paramsTable$param
#' params <- data.frame(params)
#' require(terra)
#' TCi <- vect(system.file("extdata/YASI/YASI.shp", package="TCHazaRds"))
#' TCi$PRES <- TCi$BOM_PRES
#' t1 <- strptime("2011-02-01 09:00:00","%Y-%m-%d %H:%M:%S", tz = "UTC") #first date in POSIX format
#' t2 <- strptime(rev(TCi$ISO_TIME)[1],"%Y-%m-%d %H:%M:%S", tz = "UTC") #last date in POSIX format
#' outdate <- seq(t1,t2,"hour") #array sequence from t1 to t2 stepping by “hour”

#' TCil = update_Track(outdate = outdate,
#'                    indate = strptime(TCi$ISO_TIME,"%Y-%m-%d %H:%M:%S", tz = "UTC"),
#'                    TClons = TCi$LON,
#'                    TClats = TCi$LAT,
#'                    vFms=TCi$STORM_SPD,
#'                   thetaFms=TCi$thetaFm,
#'                    cPs=TCi$PRES,
#'                   rMaxModel=params$rMaxModel,
#'                    vMaxModel=params$vMaxModel,
#'                    betaModel=params$betaModel,
#'                    eP = params$eP,
#'                    rho = params$rhoa)
#'
update_Track <- function(outdate = NULL, indate, TClons, TClats, vFms, thetaFms, cPs,
                         rMaxModel, vMaxModel, betaModel, eP, rho = NULL) {
  TRACK <- list()
  TRACK$eP <- eP
  TRACK$rho <- rho

  odatei <- as.numeric(indate)

  if (!is.null(outdate[1])) { # interpolate track to out time steps with approx fun
    indatei <- as.numeric(indate)
    outdatei <- as.numeric(outdate)
    odatei <- stats::approx(indatei, indatei, outdatei)$y
    TClons <- stats::approx(indatei, TClons, outdatei)$y
    TClats <- stats::approx(indatei, TClats, outdatei)$y
    vFms <- stats::approx(indatei, vFms, outdatei)$y
    thetaFms <- stats::approx(indatei, thetaFms, outdatei)$y
    cPs <- stats::approx(indatei, cPs, outdatei)$y
  }

  s <- !is.na(cPs) # which values to keep
  if (length(s) < 1) {
    stop("A non-NA Central Pressure value is required")
  }

  TRACK$TClons <- TClons[s]
  TRACK$TClats <- TClats[s]
  TRACK$vFms <- vFms[s]
  TRACK$thetaFms <- thetaFms[s]
  TRACK$cPs <- cPs[s]
  TRACK$odatei <- odatei[s]

  TRACK$dPs <- (eP - TRACK$cPs) # pressure deficit assumed to be 1013 hPa

  TRACK$dt <- c(diff(TRACK$odatei), 0)
  TRACK$dPdt <- c(diff(TRACK$cPs), 0) / (TRACK$dt / 3600) # hPa per hour from Holland 2008
  if (length(s) == 1) {
    TRACK$dt <- 3600
    TRACK$dPdt <- 0.1
  }
  TRACK$dPdt[is.nan(TRACK$dPdt)] <- 0

  # Radius of maximum winds model selection
  if (length(rMaxModel) > 1) {
    TRACK$rMax <- rMaxModel # use the input values
  }
  if (length(rMaxModel) == 1) {
    TRACK$rMax <- rMax_modelsR(rMaxModel = rMaxModel, TClats = TRACK$TClats, cPs = TRACK$cPs, eP = eP,
                               dPdt = TRACK$dPdt, vFms = TRACK$vFms, rho = rho)
  }

  # Compute the Coriolis parameter
  TClatrad <- TRACK$TClats * pi / 180.0
  wearth <- pi * (1.0 / 24.0) / 1800.0
  TRACK$f <- 2.0 * wearth * sin(TClatrad)

  # Maximum wind speed models selection
  if (length(vMaxModel) > 1) {
    TRACK$vMax <- vMaxModel
  }
  if (length(vMaxModel) == 1) {
    TRACK$vMax <- vMax_modelsR(vMaxModel = vMaxModel, cPs = TRACK$cPs, eP = eP, vFms = TRACK$vFms,
                               TClats = TRACK$TClats, dPdt = TRACK$dPdt, beta = 1.3, rho = rho) # Holland 80 beta = 1.3
  }

  # Beta parameter model selection
  if (length(betaModel) > 1) {
    beta <- betaModel
  }
  TRACK$beta <- beta_modelsR(betaModel = betaModel, vMax = TRACK$vMax, rMax = TRACK$rMax, cPs = TRACK$cPs,
                             eP = TRACK$eP, vFms = TRACK$vFms, TClats = TRACK$TClats, dPdt = TRACK$dPdt)

  return(TRACK)
}



#' Temporally Interpolate Along a Tropical Cyclone Track And Compute Along-Track Parameters
#'
#' @param outdate POSIX times to be interpolated to. The output date in "YYYY-MM-DD" format. Default is NULL.
#' @param TC SpatVector of Tropical cyclone track parameters
#' @param paramsTable Global parameters to compute TC Hazards.
#'
#' @return SpatVector of Tropical cyclone track parameters
#' @export
#'
#' @examples
#' require(terra)
#' TCi <- vect(system.file("extdata/YASI/YASI.shp", package="TCHazaRds"))
#' TCi$PRES <- TCi$BOM_PRES
#' TCi$PRES[is.na(TCi$PRES)] = 1010
#' outdate = seq(strptime(TCi$ISO_TIME[1],"%Y-%m-%d %H:%M:%S",tz="UTC"),
#' strptime(rev(TCi$ISO_TIME)[1],"%Y-%m-%d %H:%M:%S",tz="UTC"),3600)
#' paramsTable = read.csv(system.file("extdata/tuningParams/defult_params.csv",package = "TCHazaRds"))
#' TCii = TCvectInterp(outdate = outdate,TC=TCi,paramsTable = paramsTable)
#'
TCvectInterp = function(outdate = NULL, TC, paramsTable) {
  # Extract and organize the parameters from the paramsTable
  params <- array(paramsTable$value, dim = c(1, length(paramsTable$value)))
  colnames(params) <- paramsTable$param
  params <- data.frame(params)

  # Convert the date in ISO_TIME to POSIX class
  indate <- strptime(TC$ISO_TIME, "%Y-%m-%d %H:%M:%S", tz = "UTC")

  # Interpolate and reformat the track data if outdate is provided
  TRACK <- update_Track(outdate = outdate, indate = indate, TClons = TC$LON, TClats = TC$LAT,
                        vFms = TC$STORM_SPD, thetaFms = TC$thetaFm, cPs = TC$PRES,
                        rMaxModel = params$rMaxModel, vMaxModel = params$vMaxModel,
                        betaModel = params$betaModel, eP = params$eP, rho = params$rhoa)

  # Create a spatial vector from the interpolated track data
  v <- terra::vect(cbind(TRACK$TClons, TRACK$TClats))
  v$NAME <- TC$NAME
  v$date <- paste(TRACK$odatei + strptime("1970-01-01 00:00", "%Y-%m-%d %H:%M", tz = "UTC"))
  v$rMax <- TRACK$rMax
  v$vMax <- TRACK$vMax
  v$beta <- TRACK$beta
  v$dPdt <- TRACK$dPdt
  v$cP <- TRACK$cPs
  v$TClat <- TRACK$TClats
  v$vFm <- TRACK$vFm

  return(v)
}



#' Compute the Hazards Associated Over the Period of a TCs Event at one Given Location
#'
#' @param outdate array of POSITx date times to linearly interpolate TC track,optional.
#' @param GEO_land SpatVector or dataframe hazard geometry generated with land_geometry
#' @param TC SpatVector of Tropical cyclone track parameters
#' @param paramsTable Global parameters to compute TC Hazards.
#'
#' @return list() containing a timeseries
#'
#'
#' |  abbreviated attribute       |  description           |  units  |
#' | ------------- | ------------- | ------------- |
#' | date      | POSIX data time object of TC or outdate if provided  |  as.POSIX  |
#' | P      | Atmospheric pressure |  hPa  |
#' | Uw      | Meridional  wind speed  |  m/s  |
#' | Vw      | Zonal wind speed  |  m/s  |
#' | Sw      | Wind speed  |  m/s  |
#' | R      | distance to TC centre  |  m  |
#' | rMax      | radius of maximum wind  |  km  |
#' | vMax      | TC maximum velocity  |  m/s  |
#' | b      | TC wind profile exponent  |  -  |
#' | CP      | TC central Pressure  |  hPa  |
#' | dPdt      | change in TC CP per hour  |  hPa/hr  |
#' | vFm      | velocity of TC forward motion  |  m/s  |
#'
#' @details The function calculates wind speed and direction time series from a tropical cyclone track using various wind profile models.
#' @md
#'
#' @export
#'
#' @examples
#' GEO_land = data.frame(dem=0,lons = 147,lats=-18,f=-4e-4,inlandD = 0)
#'
#' require(terra)
#' TCi <- vect(system.file("extdata/YASI/YASI.shp", package="TCHazaRds"))
#' TCi$PRES <- TCi$BOM_PRES
#'
#' paramsTable = read.csv(system.file("extdata/tuningParams/defult_params.csv",package = "TCHazaRds"))
#' HAZts = TCHazaRdsWindTimeSereies(GEO_land=GEO_land,TC=TCi,paramsTable = paramsTable)
#' main =  paste(TCi$NAME[1],TCi$SEASON[1],"at",GEO_land$lons,GEO_land$lats)
#' with(HAZts,plot(date,Sw,format = "%b-%d %H",type="l",main = main,ylab = "Wind speed [m/s]"))
TCHazaRdsWindTimeSereies <- function(outdate = NULL, GEO_land = NULL, TC, paramsTable) {
  # Extract parameters from the paramsTable and convert them into a data frame.
  params <- array(paramsTable$value, dim = c(1, length(paramsTable$value)))
  colnames(params) <- paramsTable$param
  params <- data.frame(params)

  # Convert ISO_TIME to POSIX format
  indate <- strptime(TC$ISO_TIME, "%Y-%m-%d %H:%M:%S", tz = "UTC")

  # Reformat and interpolate the track if outdate is provided
  TRACK <- update_Track(outdate = outdate, indate = indate, TClons = TC$LON, TClats = TC$LAT,
                        vFms = TC$STORM_SPD, thetaFms = TC$thetaFm, cPs = TC$PRES,
                        rMaxModel = params$rMaxModel, vMaxModel = params$vMaxModel,
                        betaModel = params$betaModel, eP = params$eP, rho = params$rhoa)

  # Extract geographical land information
  lon <- GEO_land$lons
  lat <- GEO_land$lats

  # Calculate distance from the center to each grid point
  Rlam <- with(TRACK, RdistPi(Gridlon = lon, Gridlat = lat, TClon = TClons, TClat = TClats))
  R <- Rlam[, 1]

  # Calculate pressure based on pressure profile model
  if (params$pressureProfileModel == 0)
    P <- with(TRACK, HollandPressureProfilePi(rMax = rMax, dP = dPs, cP = cPs, beta = beta, R = R))
  if (params$pressureProfileModel == 2)
    P <- with(TRACK, DoubleHollandPressureProfilePi(rMax = rMax, dP = dPs, cP = cPs, beta = beta, R = R))

  # Calculate wind speed and direction based on wind profile model
  fs <- rep(GEO_land$f, length(R))
  if (params$windProfileModel == 0)
    VZ <- with(TRACK, HollandWindProfilePi(f = fs, vMax = vMax, rMax = rMax, dP = dPs, rho = rho, beta = beta, R = R))
  if (params$windProfileModel == 1)
    VZ <- with(TRACK, NewHollandWindProfilePi(f = fs, vMax = vMax, rMax = rMax, dP = dPs, rho = rho, beta = beta, R = R))
  if (params$windProfileModel == 2)
    VZ <- with(TRACK, DoubleHollandWindProfilePi(f = fs, vMax = vMax, rMax = rMax, dP = dPs, rho = rho, beta = beta, R = R, cP = cPs))
  if (params$windProfileModel == 4)
    VZ <- with(TRACK, JelesnianskiWindProfilePi(f = fs, vMax = vMax, rMax = rMax, R = R))

  V <- VZ[, 1]

  # Calculate wind vortex model
  if (params$windVortexModel == 0)
    UV <- with(TRACK, KepertWindFieldPi(rMax = rMax, vMax = vMax, vFm = vFms, thetaFm = thetaFms, f = fs, Rlam = Rlam, VZ = VZ))
  if (params$windVortexModel == 1)
    UV <- with(TRACK, HubbertWindFieldPi(rMax = rMax, vFm = vFms, thetaFm = thetaFms, f = fs, Rlam = Rlam, V = V))
  if (params$windVortexModel == 2)
    UV <- with(TRACK, McConochieWindFieldPi(rMax = rMax, vMax = vMax, vFm = vFms, thetaFm = thetaFms, Rlam = Rlam, V = V, f = fs[1]))

  Uw <- UV[, 1]
  Vw <- UV[, 2]

  # Calculate wind direction from wind components
  Va <- atan2(Vw, Uw)
  Dw <- 90 - (Va - pi) * 180 / pi # Convert to clockwise from true north
  Dw[Dw < 0 & !is.na(Dw)] <- 360 + Dw[Dw < 0 & !is.na(Dw)]
  Dw[Dw > 360 & !is.na(Dw)] <- Dw[Dw > 360 & !is.na(Dw)] - 360

  # Calculate wind decay based on parameters
  decay <- 1
  a <- c(params$Decay_a1, params$Decay_a2, params$Decay_a3)
  if (a[1] != 0 & params$surface == 1) {
    decay <- inlandWindDecay(GEO_land$inlandD / 1000, a)
    decay[is.na(decay)] <- 1
  }

  # Bias-correct winds and correct for surface roughness
  Sw <- sqrt(Uw^2 + Vw^2) * decay

  Uw <- cos(Va) * Sw
  Vw <- sin(Va) * Sw

  # Create and populate the output list
  v <- list() # time series list
  v$date <- TRACK$odatei + strptime("1970-01-01 00:00", "%Y-%m-%d %H:%M", tz = "UTC")
  v$Uw <- cos(Va) * Sw
  v$Vw <- sin(Va) * Sw
  v$Sw <- Sw
  v$Dw <- Dw
  v$P <- P
  v$R <- R
  v$rMax <- TRACK$rMax
  v$vMax <- TRACK$vMax
  v$beta <- TRACK$beta
  v$dPdt <- TRACK$dPdt
  v$cP <- TRACK$cPs
  v$TClat <- TRACK$TClats
  v$vFm <- TRACK$vFm

  return(v)
}

#' Compute the Wind and Pressure Spatial Hazards Field Associated with TCs Single Time Step.
#'
#' @param GEO_land SpatVector or dataframe hazard geometry generated with land_geometry
#' @param TC SpatVector or data.frame of Tropical cyclone track parameters for a single time step.
#' @param paramsTable Global parameters to compute TC Hazards.
#'
#' @return SpatRaster with the following attributes
#'
#'
#'
#' | abbreviated attribute       | description     | units |
#' | ------------- | -------------|  -------------|
#' | P      | Atmospheric pressure | hPa  |
#' | Uw      | Meridional  wind speed | m/s |
#' | Vw      | Zonal wind speed | m/s  |
#' | Sw     | Wind speed | m/s  |
#' | Dw     | Wind direction | deg clockwise from true north  |
#'
#' @md
#'
#' @export
#'
#' @examples
#' require(terra)
#' dem <- rast(system.file("extdata/DEMs/YASI_dem.tif", package="TCHazaRds"))
#' land <- dem; land[land > 0] = 0
#' inland_proximity = distance(land,target = 0)
#' GEO_land = land_geometry(dem,inland_proximity)
#'
#' TCi = vect(cbind(c(154,154),c(-26.1,-26)),"lines",crs="epsg:4283") #track line segment
#' TCi$PRES = 950
#' TCi$RMW = 40
#' TCi$ISO_TIME = "2022-10-04 20:00:00"
#' TCi$LON = geom(TCi)[1,3]
#' TCi$LAT = geom(TCi)[1,4]
#' TCi$STORM_SPD = perim(TCi)/(3*3600) #m/s
#' TCi$thetaFm = 90-returnBearing(TCi)
#' #OR
#' TC <- vect(system.file("extdata/YASI/YASI.shp", package="TCHazaRds"))
#' TC$PRES <- TC$BOM_PRES
#' TCi = TC[47]
#' plot(dem);lines(TCi,lwd = 4,col=2)
#'
#' paramsTable = read.csv(system.file("extdata/tuningParams/defult_params.csv",package = "TCHazaRds"))
#' #calculate the wind hazard
#' HAZ = TCHazaRdsWindField(GEO_land,TCi,paramsTable)
#' plot(HAZ)
#'
#' #require(rasterVis) #pretty spatial vector plot
#' #ats = seq(0, 80, length=9)
#' #UV = as(c(HAZ["Uw"],HAZ["Vw"]),"Raster") #need to convert back to raster
#' #vectorplot(UV, isField='dXY', col.arrows='white', aspX=0.002,aspY=0.002,at=ats ,
#' #colorkey=list( at=ats), par.settings=viridisTheme)
#'
TCHazaRdsWindField <- function(GEO_land, TC, paramsTable) {
  # Extract parameters from paramsTable
  params <- array(paramsTable$value, dim = c(1, length(paramsTable$value)))
  colnames(params) <- paramsTable$param
  params <- data.frame(params)

  # Convert TC dates to POSIX class if necessary
  if (methods::is(TC, "SpatVector")) {
    indate <- strptime(TC$ISO_TIME, "%Y-%m-%d %H:%M:%S", tz = "UTC")
    # Reformat track data
    TRACK <- update_Track(outdate = NULL, indate = indate, TClons = TC$LON, TClats = TC$LAT, vFms = TC$STORM_SPD, thetaFms = TC$thetaFm,
                          cPs = TC$PRES, rMaxModel = params$rMaxModel, vMaxModel = params$vMaxModel, betaModel = params$betaModel,
                          eP = params$eP, rho = params$rhoa)
  } else if (methods::is(TC, "data.frame")) {
    TRACK <- TC
    indate <- strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S", tz = "UTC") + TRACK$odatei
  }

  # Check if the central pressure value is provided
  if (is.na(TRACK$cPs)) {
    stop("Central pressure value must be provided")
  }

  lon <- as.numeric(values(GEO_land$lons))
  lat <- as.numeric(values(GEO_land$lats))
  Rlam <- with(TRACK, Rdist(Gridlon = lon, Gridlat = lat, TClon = TClons, TClat = TClats))
  R <- Rlam[, 1]

  # Calculate pressure profile based on the selected model
  if (params$pressureProfileModel == 0) {
    P <- with(TRACK, HollandPressureProfile(rMax = rMax, dP = dPs, cP = cPs, beta = beta, R = R))
  } else if (params$pressureProfileModel == 2) {
    P <- with(TRACK, DoubleHollandPressureProfile(rMax = rMax, dP = dPs, cP = cPs, beta = beta, R = R))
  }

  # Calculate wind profile based on the selected model
  fs <- terra::values(GEO_land$f)
  if (params$windProfileModel == 0) {
    VZ <- with(TRACK, HollandWindProfile(f = f, vMax = vMax, rMax = rMax, dP = dPs, rho = rho, beta = beta, R = R))
  } else if (params$windProfileModel == 1) {
    VZ <- with(TRACK, NewHollandWindProfile(f = f, vMax = vMax, rMax = rMax, dP = dPs, rho = rho, beta = beta, R = R))
  } else if (params$windProfileModel == 2) {
    VZ <- with(TRACK, DoubleHollandWindProfile(f = f, vMax = vMax, rMax = rMax, dP = dPs, rho = rho, beta = beta, R = R, cP = cPs))
  } else if (params$windProfileModel == 4) {
    VZ <- with(TRACK, JelesnianskiWindProfile(f = f, vMax = vMax, rMax = rMax, R = R))
  }

  V <- VZ[, 1]

  # Calculate wind vortex field based on the selected model
  if (params$windVortexModel == 0) {
    UV <- with(TRACK, KepertWindField(rMax = rMax, vMax = vMax, vFm = vFms, thetaFm = thetaFms, f = f, Rlam = Rlam, VZ = VZ, surface = params$surface))
  } else if (params$windVortexModel == 1) {
    UV <- with(TRACK, HubbertWindField(rMax = rMax, vFm = vFms, thetaFm = thetaFms, f = f, Rlam = Rlam, V = V, surface = params$surface))
  } else if (params$windVortexModel == 2) {
    UV <- with(TRACK, McConochieWindField(rMax = rMax, vMax = vMax, vFm = vFms, thetaFm = thetaFms, Rlam = Rlam, V = V, f = f, surface = params$surface))
  }

  # Extract wind components
  ept <- GEO_land["lons"] / GEO_land["lons"]
  Pr <- ept
  terra::values(Pr) <- P
  Uw <- ept
  terra::values(Uw) <- UV[, 1]
  Vw <- ept
  terra::values(Vw) <- UV[, 2]

  # Calculate wind speed and wind direction
  Sw <- sqrt(Uw^2 + Vw^2)
  Va <- terra::atan2(Vw, Uw)
  Dw <- 90 - (Va - pi / 2) * 180 / pi
  Dwv <- terra::values(Dw)

  # Adjust wind direction to be within the range [0, 360]
  s <- which(Dwv < 0 & !is.na(Dwv))
  if (length(s) > 0) {
    Dw[s] <- 360 + Dw[s]
  }
  s <- which(Dwv > 360 & !is.na(Dwv))
  if (length(s) > 0) {
    Dw[s] <- Dw[s] - 360
  }

  # Calculate wind decay based on surface roughness
  decay <- 1
  a <- c(params$Decay_a1, params$Decay_a2, params$Decay_a3)
  if (a[1] != 0 & params$surface == 1) {
    decay <- inlandWindDecay(GEO_land$inlandD / 1000, a)
    decay[is.na(decay)] <- 1
  }

  # Bias-correct winds and correct for surface roughness
  Sw <- Sw * decay
  Uw <- cos(Va) * Sw
  Vw <- sin(Va) * Sw

  # Create a raster stack with wind field components
  terra::time(Pr) <- indate
  terra::units(Pr) <- "hPa"
  terra::longnames(Pr) <- "air_pressure_at_sea_level"

  terra::time(Uw) <- indate
  terra::units(Uw) <- "m/s"
  terra::longnames(Uw) <- "eastward_wind"

  terra::time(Vw) <- indate
  terra::units(Vw) <- "m/s"
  terra::longnames(Vw) <- "northward_wind"

  terra::time(Sw) <- indate
  terra::units(Sw) <- "m/s"
  terra::longnames(Sw) <- "wind_speed"

  terra::time(Dw) <- indate
  terra::units(Dw) <- "Deg"
  terra::longnames(Dw) <- "wind_direction"

  rl <- list(Pr, Uw, Vw, Sw, Dw)
  rout <- terra::rast(rl)
  names(rout) <- c("Pr", "Uw", "Vw", "Sw", "Dw")

  return(rout)
}

#' Compute the Wind and Pressure Spatial Hazards Field Associated with TC track.
#'
#' @param outdate array of POSITx date times to linearly interpolate TC track
#' @param GEO_land SpatVector or dataframe hazard geometry generated with land_geometry
#' @param TC SpatVector of Tropical cyclone track parameters for a single time step
#' @param paramsTable Global parameters to compute TC Hazards
#' @param outfile character. Output netcdf filename
#' @param overwrite TRUE/FALSE, option to overwrite outfile
#' @return SpatRasterDataset with the following attributes.
#'
#'
#'
#' | abbreviated attribute       | description     | units |
#' | ------------- | -------------|  -------------|
#' | P      | Atmospheric pressure | hPa  |
#' | Uw      | Meridional  wind speed | m/s |
#' | Vw      | Zonal wind speed | m/s  |
#' | Sw     | Wind speed | m/s  |
#' | Dw     | Wind direction | deg clockwise from true north  |
#'
#' @md
#'
#' @export
#'
#' @examples
#' require(terra)
#' dem <- rast(system.file("extdata/DEMs/YASI_dem.tif", package="TCHazaRds"))
#' land <- dem; land[land > 0] = 0
#' inland_proximity = distance(land,target = 0)
#' GEO_land = land_geometry(dem,inland_proximity)
#'
#' TCi = vect(cbind(c(154,154),c(-26.1,-26)),"lines",crs="epsg:4283") #track line segment
#' TCi$PRES = 950
#' TCi$RMW = 40
#' TCi$ISO_TIME = "2022-10-04 20:00:00"
#' TCi$LON = geom(TCi)[1,3]
#' TCi$LAT = geom(TCi)[1,4]
#' TCi$STORM_SPD = perim(TCi)/(3*3600) #m/s
#' TCi$thetaFm = 90-returnBearing(TCi)
#' #OR
#' TC <- vect(system.file("extdata/YASI/YASI.shp", package="TCHazaRds"))
#' TC$PRES <- TC$BOM_PRES
#' plot(dem);lines(TC,lwd = 4,col=2)
#'
#' paramsTable = read.csv(system.file("extdata/tuningParams/defult_params.csv",package = "TCHazaRds"))
#' #calculate the wind hazard
#' HAZ = TCHazaRdsWindFields(GEO_land=GEO_land,TC=TC,paramsTable=paramsTable)
#' plot(min(HAZ$Pr))
#'
#' outdate = seq(strptime(TC$ISO_TIME[1],"%Y-%m-%d %H:%M:%S",tz="UTC"),
#'               strptime(rev(TC$ISO_TIME)[1],"%Y-%m-%d %H:%M:%S",tz="UTC"),
#'               3600)
#' HAZi = TCHazaRdsWindFields(outdate=outdate,GEO_land=GEO_land,TC=TC,paramsTable=paramsTable)
#' plot(min(HAZi$Pr))
#'
TCHazaRdsWindFields <- function(outdate = NULL, GEO_land, TC, paramsTable, outfile = NULL, overwrite = FALSE) {
  # Extract and format the parameters
  params <- data.frame(t(paramsTable$value))
  colnames(params) <- paramsTable$param

  # Convert time string to datetime format
  indate <- as.POSIXct(TC$ISO_TIME, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

  # Reformat and interpolate track if outdate is provided
  TRACK <- as.data.frame(update_Track(outdate = outdate, indate = indate, TClons = TC$LON, TClats = TC$LAT,
                                      vFms = TC$STORM_SPD, thetaFms = TC$thetaFm, cPs = TC$PRES,
                                      rMaxModel = params$rMaxModel, vMaxModel = params$vMaxModel,
                                      betaModel = params$betaModel, eP = params$eP, rho = params$rhoa))

  nt <- nrow(TRACK)
  TRACK$PRES <- TRACK$cPs
  TRACK$STORM_SPD <- TRACK$vFms
  TRACK$LON <- TRACK$TClons
  TRACK$LAT <- TRACK$TClons
  TRACK$thetaFm <- TRACK$thetaFms

  # Compute wind fields for each time step
  s <- which(!is.na(TRACK$cPs))
  HAZ_l <- lapply(s, function(x) TCHazaRdsWindField(GEO_land = GEO_land, TC = TRACK[x,], paramsTable = paramsTable))

  # Extract wind field components and create spatial data sets
  Pr <- terra::rast(lapply(HAZ_l, function(x) x$Pr))
  Uw <- terra::rast(lapply(HAZ_l, function(x) x$Uw))
  Vw <- terra::rast(lapply(HAZ_l, function(x) x$Vw))
  Sw <- terra::rast(lapply(HAZ_l, function(x) x$Sw))
  Dw <- terra::rast(lapply(HAZ_l, function(x) x$Dw))

  # Combine spatial data sets into a single sds object
  HAZs <- terra::sds(list(Pr = Pr, Uw = Uw, Vw = Vw, Sw = Sw, Dw = Dw))
  terra::varnames(HAZs) <- c("Pr", "Uw", "Vw", "Sw", "Dw")
  terra::longnames(HAZs) <- c("air_pressure_at_sea_level", "eastward_wind", "northward_wind", "wind_speed", "wind_direction")
  terra::units(HAZs) <- c("hPa", "m/s", "m/s", "m/s", "m/s")

  # Write output to a file if provided
  if (!is.null(outfile)) {
    terra::writeCDF(HAZs, filename = outfile, overwrite = overwrite)
  }

  return(HAZs)
}


#' Return the Bearing for Line Segments
#'
#' @param x spatial vector with line segments (two connected points)
#'
#' @return array of bearings see geosphere::bearing, i.e the Forward direction of the storm geographic bearing, positive clockwise from true north
#' @export
#' @import geosphere, terra
#'
#' @examples
#' ### IBTRACS HAS the WRONG BEARING!!
#' require(terra)
#' northwardTC <- vect(cbind(c(154,154),c(-26.1,-26)),"lines",crs="epsg:4283") #track line segment
#' easthwardTC <- vect(cbind(c(154,154.1),c(-26,-26)),"lines",crs="epsg:4283") #track line segment
#' southhwardTC <- vect(cbind(c(154,154),c(-26,-26.1)),"lines",crs="epsg:4283") #track line segment
#' westwardTC <- vect(cbind(c(154.1,154),c(-26,-26)),"lines",crs="epsg:4283") #track line segment
#' returnBearing(northwardTC)
#' returnBearing(easthwardTC)
#' returnBearing(southhwardTC)
#' returnBearing(westwardTC)
returnBearing <- function(x){
  xy = terra::geom(x)
  np = max(xy[,1])
  bear = sapply(1:np, function(i) {sxy = xy[xy[,1] == i,]; geosphere::bearing(sxy[1,3:4],sxy[2,3:4])})
  bear[bear < 0] = bear[bear < 0]+360
  return(bear)
}





#' Title Read All netcdf Variables Into a list()

#'
#' @param infile Name of the existing netCDF file to be opened.
#'
#' @return a list of lists with netCDF variables read by the ncdf4 package.
#' @export
#'
readncs = function(infile){
  print(paste("reading",infile))
  nc = ncdf4::nc_open(infile)
  nam = names(nc$var)
  out = list()
  out$name = infile
  for(nami in 1:length(nam)) out[[nam[nami]]] = ncdf4::ncvar_get(nc,nam[nami])
  out$dim = nc$dim
  ncdf4::nc_close(nc)
  print("read")
  return(out)
}


#' Compute the  Tropical Cyclone Radius of Maximum Winds
#'
#' @param rMaxModel 0=Powell et.al.(2005),1=McInnes et.al.(2014),2=Willoughby & Rahn (2004),  3=Vickery & Wadhera (2008), 4=Takagi & Wu (2016), 5 = Chavas & Knaff (2022)
#' @param TClats Tropical cyclone central latitude (nautical degrees)
#' @param cPs Tropical cyclone central pressure (hPa)
#' @param eP Background environmental pressure (hPa)
#' @param R175ms radius of 17.5m/s wind speeds (km)
#' @param dPdt rate of change in central pressure over time, hPa per hour from Holland 2008
#' @param vFms Forward speed of the storm m/s
#' @param rho  density of air
#'
#' @return radius of maximum winds (km)
#' @export
#'
#' @examples rMax_modelsR(0,-14,950,1013,200,0,0,1.15)
rMax_modelsR <- function(rMaxModel, TClats, cPs, eP, R175ms = 150, dPdt = NULL, vFms = NULL, rho = 1.15) {
  #not in TCRM
  #0: Powell Soukup et al (2005) updated to Arthur 2021
  #1: McInnes et al 2014 (Kossin, pers. comm. October, 2010).
  #2: Willoughby & Rahn(2004), eq 7
  #3: Vickery & Wadhera (2008) eq 11
  #4: Takagi & Wu (2016)
  #5: Chavas & Knaff (2022)
  dP <- (eP - cPs)
  dP[dP < 1] <- 1

  if (rMaxModel == 0) {
    # Powell Soukup et al (2005) updated to Arthur 2021
    rMaxs <- exp(3.543 - 0.00378 * dP + 0.813 * exp(-0.0022 * (dP)^2) + 0.00157 * (TClats^2))
  } else if (rMaxModel == 1) {
    # McInnes et al 2014 (Kossin, pers. comm. October, 2010)
    rMaxs <- exp(-3.5115 + 0.0068 * cPs + 0.0264 * abs(TClats))
  } else if (rMaxModel == 2) {
    # Willoughby & Rahn(2004), eq 7
    a <- 0.0173 * sqrt(dP * 100 / (exp(1) * rho)) + 0.0313 * (-0.0223 * sqrt(dP * 100 / (exp(1) * rho)) + 0.0281 * abs(TClats))
    b <- 1.0036 - 0.0313 * log(51.6) + 0.0087 * abs(TClats)
    beta <- (1/2) * (a^2 + sqrt(a^4 + 4 * a^2 * b) + 2 * b)
    beta[beta < 0.8] <- 0.8
    beta[beta > 1.9] <- 1.9
    vMax <- sqrt(beta * dP * 100 / (exp(1) * rho))
    rMaxs <- 51.6 * exp(-0.0223 * vMax + 0.0281 * abs(TClats))
  } else if (rMaxModel == 3) {
    # Vickery Wadhera 2008 eq 11
    rMaxs <- exp(3.015 - 6.291e-5 * dP^2 + 0.0337 * abs(TClats))
  } else if (rMaxModel == 4) {
    # Takagi Wu 2016 Figure 3
    rMaxs <- 0.676 * cPs - 578
  } else if (rMaxModel == 5) {
    x <- 0.6 * (1 - dP / 215)
    bs <- -4.4e-5 * dP^2 + 0.01 * dP + 0.03 * dPdt - 0.014 * abs(TClats) + 0.15 * vFms^x + 1
    vMax <- 0.6252 * sqrt(dP * 100)
    f <- 2 * 7.292e-5 * sin(abs(TClats * pi / 180))
    ruser_m_vec <- R175ms * 1000
    Muser_vec <- ruser_m_vec * 17.491 + 0.5 * abs(f) * ruser_m_vec^2
    halffcorruser_m_vec <- 0.5 * f * ruser_m_vec
    coefs.b <- 0.699
    coefs.c <- -0.00618
    coefs.f <- -0.00210
    MmaxMuser_predict <- coefs.b * exp(coefs.c * (vMax - 17.491) + coefs.f * (vMax - 17.491) * halffcorruser_m_vec)
    Mmax_predict <- MmaxMuser_predict * Muser_vec
    rmax_predict <- (vMax / f) * (sqrt(1 + (2 * f * Mmax_predict / (vMax^2))) - 1)
    rMaxs <- rmax_predict / 1000
  }

  return(rMaxs)
}

#' Compute the Tropical Cyclone Maximum Wind Speeds
#'
#' @param vMaxModel 0=Arthur (1980),1=Holland (2008),2=Willoughby & Rahn (2004).3=Vickery & Wadhera (2008),4=Atkinson and Holliday (1977)
#' @param cPs Tropical cyclone central pressure (hPa)
#' @param eP Background environmental pressure (hPa)
#' @param vFms Forward speed of the storm m/s
#' @param TClats Tropical cyclone central latitude
#' @param dPdt rate of change in central pressure over time, hPa per hour from Holland 2008
#' @param beta exponential term for Holland vortex
#' @param rho density of air
#'
#' @return maximum wind speed m/s.
#' @export
#'
#' @examples
#' vMax_modelsR(vMaxModel=1,cPs=950,eP=1010,vFms = 1,TClats = -14,dPdt = .1)
vMax_modelsR <- function(vMaxModel, cPs, eP, vFms = NULL, TClats = NULL, dPdt = NULL, beta = 1.3, rho = 1.15) {
  # Calculate max wind speed in the model. The beta parameter is only used in Holland 08, where it equals 1.3.
  # vMaxModel is the type of model used:
  # 0: Arthur 2021 Willoughby & Rahn (2004) default.
  # 1: Holland (2008)
  # 2: Willoughby & Rahn (2004) via the beta parameter
  # 3: Vickery & Wadhera (2008)
  # 4: Atkinson and Holliday (1977)

  dP = eP - cPs
  dP[dP < 1] = 1 # Needs to be at least lower than the background

  if (vMaxModel == 0) {
    # see https://github.com/GeoscienceAustralia/tcha/blob/6e6de8df2b3fd5c9287d060f1f7f1d4ff63e87a7/wind-radii/rmax_fit.py#L205
    vMax = 0.6252 * sqrt(dP * 100)
  }

  if (vMaxModel == 1) {
    x = 0.6 * (1 - dP / 215)
    bs = -4.4e-5 * dP^2 + 0.01 * dP + 0.03 * dPdt - 0.014 * abs(TClats) + 0.15 * vFms^x + 1
    vMax = sqrt(abs(100 * bs / (rho * exp(1)) * dP))  # Holland 2008 model "A Revised Hurricane Pressure – Wind Model" Equation 11.
  }

  if (vMaxModel == 2) {
    # The WR04 Vmax and Rmax equations needed to be refactored to solve for beta
    a <- 0.0173 * sqrt(dP * 100 / (exp(1) * rho)) + 0.0313 * (-0.0223 * sqrt(dP * 100 / (exp(1) * rho)) + 0.0281 * abs(TClats))
    b <- 1.0036 - 0.0313 * log(51.6) + 0.0087 * abs(TClats)
    # x - a * sqrt(x) = b
    # https://www.wolframalpha.com/input?i=x-+a*sqrt%28x%29%3Db
    beta <- (1 / 2) * (a^2 + sqrt(a^4 + 4 * a^2 * b) + 2 * b)
    beta[beta < 0.8] = 0.8
    beta[beta > 1.9] = 1.9
    vMax <- sqrt(beta * dP * 100 / (exp(1) * rho))
  }

  if (vMaxModel == 3) { # Vickery and Wadhera (2008)
    rMaxs <- exp(3.015 - 6.291e-5 * dP^2 + 0.0337 * abs(TClats)) # Equation 11
    TClatrad <- TClats * pi / 180.0
    wearth <- pi * (1.0 / 24.0) / 1800.0
    f <- 2.0 * wearth * sin(TClatrad)
    beta = 1.833 - 0.326 * sqrt(abs(f) * rMaxs * 1000) # Equation 23 rMax is in units [m] not km
    vMax = sqrt(beta * dP * 100 / (exp(1) * rho)) # Equation 4  dP in Pa, not hPa, so * 100
  }

  if (vMaxModel == 4) {
    vMax = (3.04 * ((1010.0 - cPs) / 100.0)^0.644) # Atkinson and Holliday (1977), *Tropical Cyclone Minimum SeaLevel Pressure Maximum Sustained Wind Relationship for Maximum 10m, 1 - minute wind speed. Uses ``pEnv`` as 1010 hPa.
  }

  return(vMax)
}

#' Compute the Exponential TC beta Profile-Curvature Parameter
#'
#' @param betaModel 0=Holland (2008),1=Powell (2005),2=Willoughby & Rahn (2004),3=Vickery & Wadhera (2008),4=Hubbert (1991)
#' @param vMax maximum wind speed m/s. see \code{vMax_modelsR}
#' @param rMax radius of maximum winds (km). see \code{rMax_modelsR}
#' @param cPs Tropical cyclone central pressure (hPa)
#' @param eP Background environmental pressure (hPa)
#' @param vFms Forward speed of the storm m/s
#' @param TClats Tropical cyclone central latitude
#' @param dPdt rate of change in central pressure over time, hPa per hour from Holland 2008
#' @param rho density of air
#'
#' @return exponential beta parameter
#' @export
#'
#' @examples beta_modelsR(0,10,10,960,1013,3,-15,1)
beta_modelsR <- function(betaModel, vMax, rMax, cPs, eP, vFms, TClats, dPdt, rho = 1.15) {
  # 0: Powell et al (2005) see Arthur 2021
  # 1: Holland (2008)
  # 2: Willoughby & Rahn(2004) default? https://github.com/GeoscienceAustralia/tcrm/blob/1916233f7dfdecf6a1b5f5b0d89f3eb1d164bd3e/wind/vmax.py#L96
  # 3: Vickery & Wadhera (2008)
  # 4: Hubbert (1991)

  dP <- (eP - cPs)
  dP[dP < 1] <- 1

  if (betaModel == 0) {
    beta <- 1.881093 - 0.010917 * abs(TClats) - 0.005567 * rMax # Powell (2005)
    beta[beta < 0.8] <- 0.8
    beta[beta > 1.9] <- 1.9
  } else if (betaModel == 1) {
    x <- 0.6 * (1 - dP / 215)
    bs <- -4.4e-5 * dP^2 + 0.01 * dP + 0.03 * dPdt - 0.014 * abs(TClats) + 0.15 * vFms^x + 1
    beta <- bs # Holland (2008) # cyclone08v2.for beta = 1.6bs
  } else if (betaModel == 2) {
    a <- 0.0173 * sqrt(dP * 100 / (exp(1) * rho)) + 0.0313 * (-0.0223 * sqrt(dP * 100 / (exp(1) * rho)) + 0.0281 * abs(TClats))
    b <- 1.0036 - 0.0313 * log(51.6) + 0.0087 * abs(TClats)
    beta <- (1/2) * (a^2 + sqrt(a^4 + 4 * a^2 * b) + 2 * b)
    beta[beta < 0.8] <- 0.8
    beta[beta > 1.9] <- 1.9
  } else if (betaModel == 3) {
    TClatrad <- TClats * pi / 180.0
    wearth <- pi * (1.0 / 24.0) / 1800.0
    f <- 2.0 * wearth * sin(TClatrad)
    beta <- 1.833 - 0.326 * sqrt(abs(f) * rMax * 1000) # Vickery and Wadhera (2008) equ 23 rMax is in units [m] not km
  } else if (betaModel == 4) {
    beta <- 1.5 + (980 - cPs) / 120 # Hubbert (1991)
  }

  return(pmin(pmax(beta, 0.8), 1.9)) # Clip beta values to be between 0.8 and 1.9
}



#' Reduce Winds Overland
#'
#' @param d inland distance in km
#' @param a three parameter of decay model a1,a2,a3
#'
#' @return a reduction factor Km
#' @export
#'
#' @examples
#' inlandWindDecay(10)
inlandWindDecay = function(d,a = c(0.66,1,0.4)){
  d[d < 0] = 0
  f = a[1]+a[2]*(1-a[1]/a[2])*exp(-d/a[3])
  return(f)
}



#' Transect points from a origin through a point or with a bearing and to the opposite side.
#'
#' @param TC_line origin of the transect
#' @param Through_point a point to pass through
#' @param bear the bearing
#' @param length the length of the transect in Km
#' @param step the spacing of the transect in Km
#'
#' @return spatial vector of transect profile points with distances in Km (negative for left hand side)
#' @export
#'
#' @examples
#' require(terra)
#' TCi <- vect(cbind(c(154.1,154),c(-26.1,-26)),"lines",crs="epsg:4283") #track line segment
#' TCi$PRES <- 950
#' TCi$RMW <- 40
#' TCi$ISO_TIME <- "2022-10-04 20:00:00"
#' TCi$LON <- geom(TCi)[1,3]
#' TCi$LAT <- geom(TCi)[1,4]
#' TCi$STORM_SPD <- perim(TCi)/(3*3600) #m/s
#' TCi$thetaFm <- 90-returnBearing(TCi)
#' #Through_point <- isd[isd$OID==isdsi]
#' pp <- TCProfilePts(TC_line = TCi,Through_point=NULL,bear=TCi$thetaFm+90,length =100,step=10)
#' plot(pp,"radialdist",type="continuous")
#' lines(TCi,col=2)
TCProfilePts = function(TC_line,Through_point=NULL,bear=NULL,length =200,step=2){
  xy0 = terra::geom(TC_line)[1,3:4]
  #XY0 = terra::vect(rbind(xy0),"points")

  if(is.null(Through_point) & is.null(bear)) stop("need a through point or bearing")
  if(!is.null(Through_point)) {
    xy1 = terra::geom(Through_point)[1,3:4]
    #XY1 = terra::vect(rbind(xy1),"points")
  }
  if(is.null(bear)) bear = geosphere::bearing(xy0,xy1)
  dist = seq(0,length,step)*1000
  bear = 90-bear+180
  bear[bear < 0] = bear[bear < 0]+360
  ptsR = geosphere::destPoint(p=xy0,b=bear,d= dist)
  ptsL = geosphere::destPoint(p=xy0,b=bear+180,d= dist)
  np = dim(ptsL)[1]
  out = terra::vect(rbind(ptsL[np:1,],ptsR[-1,]),"points")
  out$radialdist = c(rev(-dist),dist[-1])/1000
  return(out)
}


#' Compute the Wind and Pressure Spatial Hazards Profile Associated with TCs Single Time Step.
#'
#' @param GEO_land SpatVector or dataframe hazard geometry generated with land_geometry
#' @param TC SpatVector or data.frame of Tropical cyclone track parameters for a single time step.
#' @param paramsTable Global parameters to compute TC Hazards.
#'
#' @return SpatRaster with the following attributes
#'
#'
#'
#' | abbreviated attribute       | description     | units |
#' | ------------- | -------------|  -------------|
#' | P      | Atmospheric pressure | hPa  |
#' | Uw      | Meridional  wind speed | m/s |
#' | Vw      | Zonal wind speed | m/s  |
#' | Sw     | Wind speed | m/s  |
#' | Dw     | Wind direction | deg clockwise from true north  |
#'
#' @md
#'
#' @export
#'
#' @examples
#' require(terra)
#' dem <- rast(system.file("extdata/DEMs/YASI_dem.tif", package="TCHazaRds"))
#' land <- dem; land[land > 0] = 0
#' inland_proximity = distance(land,target = 0)
#' GEO_land = land_geometry(dem,inland_proximity)
#'
#' TCi = vect(cbind(c(154,154),c(-26.1,-26)),"lines",crs="epsg:4283") #track line segment
#' TCi$PRES = 950
#' TCi$RMW = 40
#' TCi$ISO_TIME = "2022-10-04 20:00:00"
#' TCi$LON = geom(TCi)[1,3]
#' TCi$LAT = geom(TCi)[1,4]
#' TCi$STORM_SPD = perim(TCi)/(3*3600) #m/s
#' TCi$thetaFm = 90-returnBearing(TCi)
#' #OR
#' TC <- vect(system.file("extdata/YASI/YASI.shp", package="TCHazaRds"))
#' TC$PRES <- TC$BOM_PRES
#' TCi = TC[47]
#' TCi$thetaFm = 90-returnBearing(TCi)
#'
#' #extract a profile/transect at right angles (90 degrees) from the TC heading/bearing direction
#' pp <- TCProfilePts(TC_line = TCi,bear=TCi$thetaFm+90,length =100,step=1)
#' plot(dem);lines(TCi,lwd = 4,col=2)
#' points(pp)
#' GEO_land_v = extract(GEO_land,pp,bind=TRUE,method = "bilinear")
#'
#' paramsTable = read.csv(system.file("extdata/tuningParams/defult_params.csv",package = "TCHazaRds"))
#' #calculate the wind hazard
#' HAZ = TCHazaRdsWindProfile(GEO_land_v,TCi,paramsTable)
#' plot(HAZ$radialdist,HAZ$Sw,type="l",xlab = "Radial distance [km]",ylab = "Wind speed [m/s]");grid()
#' plot(HAZ,"Sw",type="continuous")
#'
TCHazaRdsWindProfile = function(GEO_land,TC,paramsTable){
  params <- array(paramsTable$value,dim = c(1,length(paramsTable$value)))
  colnames(params) <- paramsTable$param
  params <- data.frame(params)
  #SpatVector can't "hold" POSIX class.
  if(methods::is(TC,"SpatVector")) {
    indate=strptime(TC$ISO_TIME,"%Y-%m-%d %H:%M:%S",tz = "UTC")
    #reformat
    TRACK = update_Track(outdate=NULL,indate=indate,TClons=TC$LON,TClats=TC$LAT,vFms=TC$STORM_SPD,thetaFms=TC$thetaFm,cPs=TC$PRES,
                         rMaxModel=params$rMaxModel,vMaxModel=params$vMaxModel,betaModel=params$betaModel,eP = params$eP,rho = params$rhoa)
  }
  if(methods::is(TC,"data.frame")){
    TRACK = TC
    indate = strptime("1970-01-01 00:00:00","%Y-%m-%d %H:%M:%S",tz = "UTC")+TRACK$odatei
  }

  lon = as.numeric(GEO_land$lons)
  lat = as.numeric(GEO_land$lats)
  Rlam <- with(TRACK,Rdist(Gridlon =lon,Gridlat=lat,TClon=TClons,TClat=TClats))
  R <- Rlam[,1]

  #Rcpp::sourceCpp("src/TCHazaRds.cpp")

  #0: Holland (1980)
  #2: McConochie et al. (2004) "Double-Holland"
  if(params$pressureProfileModel == 0) P <- with(TRACK,HollandPressureProfile(rMax=rMax,dP = dPs, cP=cPs,beta=beta,R=R))
  if(params$pressureProfileModel == 2) P <- with(TRACK,DoubleHollandPressureProfile(rMax=rMax,dP=dPs,cP=cPs,beta=beta,R=R))

  #windProfileModel #https://geoscienceaustralia.github.io/tcrm/docs/setup.html?highlight=jelesnianski#windProfileinterface
  #0: Holland (1980)
  #2: McConochie et al. (2004) "Double-Holland"
  #4: Jelesnianski (1966)
  fs=GEO_land$f
  if(params$windProfileModel == 0) VZ <- with(TRACK,HollandWindProfile(      f=f, vMax=vMax, rMax=rMax, dP=dPs, rho=rho, beta=beta, R=R))
  if(params$windProfileModel == 1) VZ <- with(TRACK,NewHollandWindProfile(   f=f, vMax=vMax, rMax=rMax, dP=dPs, rho=rho, beta=beta, R=R))
  if(params$windProfileModel == 2) VZ <- with(TRACK,DoubleHollandWindProfile(f=f, vMax=vMax, rMax=rMax, dP=dPs, rho=rho, beta=beta, R=R, cP=cPs))
  if(params$windProfileModel == 4) VZ <- with(TRACK,JelesnianskiWindProfile( f=f, vMax=vMax, rMax=rMax,                             R=R))

  V=VZ[,1]
  #0 Kepert & Wang (2001)
  #1 Hubbert et al (1991)
  #2 McConochie et al (2004) "Double-Holland"
  if(params$windVortexModel == 0) UV <- with(TRACK,KepertWindField(    rMax=rMax, vMax=vMax, vFm=vFms, thetaFm=thetaFms, f=f , Rlam=Rlam, VZ=VZ,  surface=params$surface))
  if(params$windVortexModel == 1) UV <- with(TRACK,HubbertWindField(   rMax=rMax,            vFm=vFms, thetaFm=thetaFms, f=f , Rlam=Rlam, V=V,    surface=params$surface))
  if(params$windVortexModel == 2) UV <- with(TRACK,McConochieWindField(rMax=rMax, vMax=vMax, vFm=vFms, thetaFm=thetaFms,       Rlam=Rlam, V=V,f=f,surface=params$surface))

  GEO_land$Pr  <- P
  GEO_land$Uw  <- UV[,1]
  GEO_land$Vw  <- UV[,2]
  GEO_land$Va  <- atan2(GEO_land$Vw,GEO_land$Uw)
  GEO_land$Dw  <- 90-(GEO_land$Va-pi/2)*180/pi
  GEO_land$Dwv <- GEO_land$Dw
  s <- which(GEO_land$Dwv < 0 & !is.na(GEO_land$Dwv))
  if( length(s > 0) ) GEO_land$Dw[s] <- 360+GEO_land$Dw[s]
  s = which(GEO_land$Dwv > 360 & !is.na(GEO_land$Dwv))
  if( length(s > 0) ) GEO_land$Dw[s] <- GEO_land$Dw[s]-360

  GEO_land$decay <- 1
  a = c(params$Decay_a1,params$Decay_a2,params$Decay_a3)
  if(a[1] != 0 & params$surface == 1){
    GEO_land$decay = inlandWindDecay(GEO_land$inlandD/1000,a)
    GEO_land$decay[is.na(GEO_land$decay)] <- 1
  }

  #bias correct winds and correct for surface roughness
  GEO_land$Sw = sqrt(GEO_land$Uw^2+GEO_land$Vw^2)*GEO_land$decay

  GEO_land$Uw = cos(GEO_land$Va) * GEO_land$Sw
  GEO_land$Vw = sin(GEO_land$Va) * GEO_land$Sw

  return(GEO_land)
}



