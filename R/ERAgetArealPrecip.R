#' Get ERA areal precipitation 
#' @description Extracts the areal precipitation for all time intervals from the ERA NetCDF file for a region. If the data are 3-hour values accumulated over 12-hour periods, then they will be deaccumulated, and negative precipitation will be set to zero
#' @param ncdfFile  Required. Name of the NetCDF file containing ERA data.
#' @param timezone Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}
#'
#' @return If unsuccessful, returns \code{FALSE}. If successful, returns a list containing the following: 
#'\tabular{ll}{
#'\bold{Name}\tab \bold{meaning}\cr
#'\code{precip}\tab 3d array (lon x lat x time) of precip (mm)\cr
#'\code{lonres}\tab longitude resolution (\eqn{^\circ}{degrees})\cr
#'\code{latres}\tab latitude resolution (\eqn{^\circ}{degrees})\cr
#'\code{minLon}\tab mininum longitude of data(\eqn{^\circ}{degrees}W)\cr
#'\code{maxLon}\tab maximum longitude of data(\eqn{^\circ}{degrees}W)\cr
#'\code{minLat}\tab minimum latitude of data(\eqn{^\circ}{degrees})\cr
#'\code{minLat}\tab maximum latitude of data(\eqn{^\circ}{degrees})\cr
#'\code{datetime}\tab the datetime of each time layer\cr
#' }
#' @export
#'
#' @examples \dontrun{threehour_precip <- ERAgetArealPrecip('_grib2netcdf-atls17-95e2cf679cd58ee9b4db4dd119a05a8d-szSVRK.nc', timezone='CST')
#' }
ERAgetArealPrecip <- function(ncdfFile, timezone='', quiet=TRUE){
  # check parameters
  if (ncdfFile == ''){
    cat('Error: must specify NetCDF file name\n')
    return(FALSE)
  }
  varName <- 'tp'
  
  # calculate hour offset between timezone and GMT
  datetime <- '2000-01-01 01:00'
  s1 <- as.POSIXct(datetime, format='%Y-%m-%d %H:%M', tz='UTC')
  s2 <- as.POSIXct(datetime, format='%Y-%m-%d %H:%M', tz=timezone)
  houroffset <- as.numeric(difftime(s2, s1, units='hours'))
  
  nc <-  RNetCDF::open.nc(ncdfFile, write=FALSE)
  nctimes <- RNetCDF::var.get.nc(nc, variable='time', unpack=TRUE)
  precip <- RNetCDF::var.get.nc(nc, variable=varName, unpack=TRUE)
  lons <- as.vector(RNetCDF::var.get.nc(nc, variable='longitude', unpack=TRUE))
  lats <- as.vector(RNetCDF::var.get.nc(nc, variable='latitude', unpack=TRUE))
  RNetCDF::close.nc(nc)
  
  # get extent and find resolution
  maxLon <- max(lons) - 360
  minLon <- min(lons) - 360
  maxLat <- max(lats)
  minLat <- min(lats)
  
  lonRes <- abs(lons[2]-lons[1])
  latRes <- abs(lats[2]-lats[1])
  
  # create dataframe of lats and lons
  num.lats <- length(lats)
  num.lons <- length(lons)
  latnums <- seq(1, num.lats)
  lonnums <- seq(1, num.lons)
  all.lons <-  rep(lons, each=num.lats)
  all.lonnums <- rep(lonnums, each=num.lats)
  all.lats <- rep(lats, times=num.lons)
  all.latnums <- rep(latnums, times=num.lons)
  all.locs <- data.frame(all.lonnums, all.lons, all.latnums, all.lats)
  names(all.locs)[c(2,4)] <- c('longitude', 'latitude')
  
  # add times and values
  xy <- as.matrix(all.locs[,c('longitude', 'latitude')], ncol=2)
  
  # convert times to R times
  secs <- nctimes * 3600.0
  datetimes.utc <- as.POSIXct(secs, origin='1900-01-01', tz='UTC')
  datetimes <- datetimes.utc + (houroffset * 3600)
  datetimes <- lubridate::force_tz(datetimes, tzone=timezone)
  
  # convert precip from m to mm
  precip <- precip * 1000
  
  # now deaccumulate
  # find time interval
  t2 <- datetimes.utc[2]
  t1 <- datetimes.utc[1]
  hourdiff <- as.numeric(difftime(t2, t1, units='hours'))
  
  if ((hourdiff != 3) & !quiet)
    cat("Couldn't disaggregate precips\n")
  else{
    # generate sequence of layers to subtract
    num.times <- length(datetimes.utc)
    layernum <- c(1,2,3,4)
    all_layernum <- rep_len(layernum, length.out=num.times)
    
    # now disaggregate
    disagg <- array(data=0, dim=c(num.lats, num.lons, num.times))
    disagg[,,all_layernum==1] <- precip[,,all_layernum==1]
    disagg[,,all_layernum==2] <- precip[,,all_layernum==2] - precip[,,all_layernum==1]
    disagg[,,all_layernum==3] <- precip[,,all_layernum==3] - precip[,,all_layernum==2]
    disagg[,,all_layernum==4] <- precip[,,all_layernum==4] - precip[,,all_layernum==3]
    
    # check for negative precip
    neglocs <- which(disagg < 0)
    if(length(neglocs) >0 ){
      if (!quiet)
        cat('Found negative precips, values are set to zero\n')
      disagg[neglocs] <- 0
    }
    precip <- disagg
  }

  # create list
  ret <- list(precip=precip, lonRes=lonRes, latRes=latRes, minLon=minLon, maxLon=maxLon, minLat=minLat, maxLat=maxLat, datetime=datetimes)
  return(ret)
  
}