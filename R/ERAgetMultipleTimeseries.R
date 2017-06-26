#' Extracts timeseries of ERA data nearest to multiple specified locations
#' @description Reads a NetCDF file containing ERA reanalysis data and extracts the timeseries of the specified variable for each location. This is faster than calling the function \code{ERAgetNearestTimseries} for each location, as the reanalysis data are only read in once. Some of the the commonly-used variables are:
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
#' @param outDir Optional. Directory to hold output files. If not specified, the current directory will be used.
#' @param siteNames Required. A vector containing the names of the sites. The names will be used for the output obs files.
#' @param pointLons Required. A vector containing decimal longitudes of desired locations. Note that the NetCDF longitude is 0-360\eqn{^\circ}{}, so add 360 to negative longitudes.
#' @param pointLats Required. A vector containing decimal latitudes of desired location.
#' @param projection Optional. Projection to be used to convert latitudes and longitudes to locations. Used for finding the nearest ERA gridpoint. The default, \option{+proj=utm +zone=13 +ellps=WGS84}, is \emph{only} valid for Western Canada. If you are processing data for the whole world, you can use the August Epicycloidal Projection which is \option{+proj=august +lon_0=90w} .
#' @param timezones Required. A vector containing the name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time. If there are fewer timezone values than variable names, the timezone will be recycled.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}, as it will list the actual latitude and longitude of the ERA tile.
#'@param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If unsuccessful, returns \code{FALSE}. If successful, returns \code{TRUE} and all variables for all locations are written to \code{obs} files in the specified directory.
#' @author Kevin Shook
#' @export
#' @seealso \code{\link{ERAgetNearestTimeseries}} \code{\link{ERAdeaccum}}
#' @examples \dontrun{
#' ERAgetMultipleLocationTimeseries('ssrd.nc', 'ssrd', s$`Station Names`,  
#' pointLons = s$Longitude, pointLats = s$Latitude, projection='+proj=august +lon_0=90w', 
#' timezones=s$timezone, quiet=FALSE)}

ERAgetMultipleLocationTimeseries <- function(ncdfFile, varName, siteNames, outDir='', pointLons, pointLats, projection='+proj=utm +zone=13 +ellps=WGS84', timezones='', quiet=TRUE, logfile=''){

  # check parameters
  if (ncdfFile == ''){
    cat('Error: must specify NetCDF file name\n')
    return(FALSE)
  }

  if (varName == ''){
    cat('Error: must specify variable name\n')
    return(FALSE)
  }

 if(outDir == ''){
   if (!quiet)
     cat('No directory specified, current directory will be used')
   # check dirwctory name for trailing slash - add if needed
   
 }
  
  numSites <- length(siteNames)
  numLons <- length(pointLons)
  numLats <- length(pointLats)
  numZones <- length(timezones)
  
  if (numLons == 0){
    cat('Error: must specify longitude\n')
    return(FALSE)
  }

  if (numLats == 0){
    cat('Error: must specify latitude\n')
    return(FALSE)
  }

  if (numZones == 0){
    cat('Error: must specify time zone\n')
    return(FALSE)
  }
  
  # check to make sure same number of values in all vectors

  if(numSites != numLons){
    cat('Error: differing number of site names and longitudes\n')
    return(FALSE)
  }
  
  if(numSites != numLats){
    cat('Error: differing number of site names and latitudes\n')
    return(FALSE)
  }
  
  if (numZones != numSites){
    if (numZones > numSites){
      if (!quiet)
        cat(numZones, ' timezones to be trimmed to ', numSites  ,' sites\n', sep='')
      timezones <- timezones[1:numSites]
    } else{
      cat(numZones, ' timezones to be recycled to ', numSites  ,' sites\n', sep='')
      timezones <- timezones[1:numSites]
    }
    
    reps <- ceiling(numSites / numZones)
    timezones <- rep.int(timezones, reps)
    timezones <- timezones[1:numSites]
  }
  
  # change poinLons if needed
  pointLons[pointLons < 0] <- pointLons[pointLons < 0] + 360

  
  nc <-  RNetCDF::open.nc(ncdfFile, write=FALSE)
  # read in values netcdf file
  for (i in 1:numSites){
    siteName <- siteNames[i]
    pointLat <- pointLats[i]
    pointLon <- pointLons[i]
    timezone <- timezones[i]
    
    nctimes <- RNetCDF::var.get.nc(nc, variable='time', unpack=TRUE)
    all.values <- RNetCDF::var.get.nc(nc, variable=varName, unpack=TRUE)
    lons <- as.vector(RNetCDF::var.get.nc(nc, variable='longitude', unpack=TRUE))
    lats <- as.vector(RNetCDF::var.get.nc(nc, variable='latitude', unpack=TRUE))
    
  
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
      cat(siteName, '\n', sep='')
      cat('                Lon      Lat\n')
      cat('Specified: ', format(pointLon, width=8), ' ', format(pointLat, width=8), '\n', sep='' )
      cat('Nearest  : ', format(closest.lon, width=8), ' ', format(closest.lat, width=8), '\n', sep='' )
    }
  
    # extract time series values
    values <- as.vector(all.values[closest.lonpoint, closest.latpoint,])
  
    # convert times to R times
    houroffset <- CRHMr::GMToffset(timezone[i])
    secs <- nctimes * 3600.0
    datetime.utc <- as.POSIXct(secs, origin='1900-01-01', tz='UTC')
    datetime = datetime.utc - (houroffset * 3600)
  
    # force timezone
    datetime <- lubridate::force_tz(datetime, tzone=timezone)
  
    obs <- data.frame(datetime, values)
    names(obs) <- c('datetime', varName)
  
    # sort by datetime, just in case
    obs <- obs[order(obs$datetime),]
    
    
    # write to file
    if (outDir != '')
      obsFilename <- paste(outDir, siteNames[i],"_", varName, ".obs", sep='')
    else
      obsFilename <- paste(siteNames[i],"_", varName, ".obs", sep='')
    
    CRHMr::writeObsFile(obs, obsfile = obsFilename, comment='Values from ERAgetMultipleTimeSeries')
    
    # output info to screen (if req'd) and write to log file
    file.info <- CRHMr::CRHM_summary(obs)
    if (!quiet)
      print(file.info)
    
  }
  RNetCDF::close.nc(nc)
  
  


  comment <- paste('ERAgetNearestTimeseries ncdfFile:',ncdfFile,
                   ' varName:', varName, sep='')
  result <- CRHMr::logAction(comment, logfile)

  if(result)
    return(TRUE)
  else
    return(result)
}
