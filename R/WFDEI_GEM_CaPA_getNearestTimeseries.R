#' Finds timeseries of WFDEI-GEM-CaPA data nearest to specified location
#' @name WFDEI-GEM-WFDEI_GEM_CaPA_getNearestTimeseries
#' @description Reads a NetCDF file containing WFDEI-GEM-CaPA reanalysis data and extracts the timeseries of the specified variable. Some of the the commonly-used variables are:
#' \tabular{llr}{
#'\bold{Parameter} \tab \bold{Units} \tab \bold{Variable}\cr
#'10 m eastward wind component \tab m/s \tab u10\cr
#'10 m northward wind component \tab m/s \tab v10\cr
#'2 metre temperature	\tab K \tab t2m \cr
#'2 metre dewpoint \tab K \tab d2m \cr
#'Downward surface solar radiation* \tab J/m\eqn{^2}{^2} \tab ssrd\cr
#'Downward surface thermal radiation* \tab J/m\eqn{^2}{^2} \tab strd\cr
#'Surface net solar radiation* \tab J/m\eqn{^2}{^2} \tab ssr \cr
#'Surface net thermal radiation* \tab J/m\eqn{^2}{^2} \tab str \cr
#'Total precipitation* \tab m of water \tab tp \cr
#' }
#' Parameters marked with an asterisk are \emph{cumulative} values, and must be deaccumulated using the deaccumERA function.
#' @param ncdfFile Required. Name of the NetCDF file containing ERA data.
#' @param varName Required. Name of the NetCDF variable to extract.
#' @param pointLon Required. Decimal longitude of desired location. Note that the NetCDF longitude is 0-360\eqn{^\circ}{}, so add 360 to negative longitudes.
#' @param pointLat Required. Decimal latitude of desired location.
#' @param projection Optional. Projection to be used to convert latitudes and longitudes to locations. Used for finding the nearest ERA gridpoint. The default, \option{+proj=utm +zone=13 +ellps=WGS84}, is \emph{only} valid for Western Canada. If you are processing data for the whole world, you can use the August Epicycloidal Projection which is \option{+proj=august +lon_0=90w} .
#' @param timezone Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}, as it will list the actual latitude and longitude of the ERA tile.
#'@param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If unsuccessful, returns \code{FALSE}. If successful, returns a standard \pkg{CRHMr} dataframe containing the datetime and the extracted data, which are unpacked (i.e. the NetCDF multiplier and offset have been applied).
#' @author Kevin Shook
#' @export
#' @seealso  \code{\link{ERAdeaccum}}
#' @examples \dontrun{
#' ta <- WFDEI-GEM-CaPAgetNearestTimeseries('1979.nc', 'ta', 69, 51.69, timezone='Etc/GMT+7')
#' p <- WFDEI-GEM-CaPAgetNearestTimeseries('1979.nc', 'pr', 243.5, 51.69, timezone='Etc/GMT+7')
#' u <- WFDEI-GEM-CaPAgetNearestTimeseries('1979.nc', 'wind_speed', 243.5, 51.69, timezone='Etc/GMT+7')
#' ps <- WFDEI-GEM-CaPAgetNearestTimeseries('1979.nc', 'd2m', 243.5, 51.69, timezone='Etc/GMT+7')}

WFDEI_GEM_CaPA_getNearestTimeseries <- function(ncdfFile, varName, pointLon, pointLat, projection='+proj=utm +zone=13 +ellps=WGS84', timezone='', quiet=TRUE, logfile=''){

  # check parameters
  if (ncdfFile == ''){
    cat('Error: must specify NetCDF file name\n')
    return(FALSE)
  }

  if (varName == ''){
    cat('Error: must specify variable name\n')
    return(FALSE)
  }

  if (pointLon == ''){
    cat('Error: must specify longitude\n')
    return(FALSE)
  }

  if (pointLat == ''){
    cat('Error: must specify latitude\n')
    return(FALSE)
  }

  if (timezone == ''){
    cat('Error: must specify time zone\n')
    return(FALSE)
  }

  # calculate hour offset between timezone and GMT
  houroffset <- CRHMr::GMToffset(timezone)

  # read in values netcdf file
  nc <-  RNetCDF::open.nc(ncdfFile, write=FALSE)
  nctimes <- RNetCDF::var.get.nc(nc, variable='time', unpack=TRUE)
  all.values <- RNetCDF::var.get.nc(nc, variable=varName, unpack=TRUE)
  lons <- as.vector(RNetCDF::var.get.nc(nc, variable='rlon', unpack=TRUE))
  lats <- as.vector(RNetCDF::var.get.nc(nc, variable='rlat', unpack=TRUE))
  RNetCDF::close.nc(nc)

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
  xy <- as.matrix(all.locs[,c('all.lons', 'all.lats')], ncol=2)

  # now project values
  projected <- proj4::project(xy, proj=projection)
  projected <- as.data.frame(projected)
  names(projected) <- c('x', 'y')
  all.locs <- cbind(all.locs, projected)

  # project point
  xy2 <- list(pointLon,pointLat)
  projected2 <- proj4::project(xy2, proj=projection)
  point_x=projected2$x
  point_y=projected2$y

  # now get distance
  all.locs$xdistance <- (all.locs$x - point_x) ^ 2
  all.locs$ydistance <- (all.locs$y - point_y) ^ 2
  all.locs$distance <- (all.locs$xdistance + all.locs$ydistance) ^ 0.5

  # get minimum
  minpoint <- which.min(all.locs$distance)
  closest.lat <- all.locs$all.lats[minpoint]
  closest.lon <- all.locs$all.lons[minpoint]
  closest.latpoint <- all.locs$all.latnums[minpoint]
  closest.lonpoint <- all.locs$all.lonnums[minpoint]

  # tell user data location
  if (!quiet){
    cat('                Lon      Lat\n')
    cat('Specified: ', format(pointLon, width=8), ' ', format(pointLat, width=8), '\n', sep='' )
    cat('Nearest  : ', format(closest.lon, width=8), ' ', format(closest.lat, width=8), '\n', sep='' )
  }

  # extract time series values
  values <- as.vector(all.values[closest.lonpoint, closest.latpoint,])

  # convert times to R times
  secs <- nctimes * 3600.0
  datetime.utc <- as.POSIXct(secs, origin='1979-01-01', tz='UTC')
  datetime = datetime.utc + (houroffset * 3600)

  # force timezone
  datetime <- lubridate::force_tz(datetime, tzone=timezone)

  obs <- data.frame(datetime, values)
  names(obs) <- c('datetime', varName)

  # sort by datetime, just in case
  obs <- obs[order(obs$datetime),]


  # output info to screen (if req'd) and write to log file
  file.info <- CRHMr::CRHM_summary(obs)
  if (!quiet)
    print(file.info)

  comment <- paste('ERAgetNearestTimeseries ncdfFile:',ncdfFile,
                   ' varName:', varName,
                   ' pointLon:', pointLon,
                   ' pointLat:', pointLat,
                   ' timezone:', timezone, sep='')
  result <- CRHMr::logAction(comment, logfile)

  if(result)
    return(obs)
  else
    return(result)
}
