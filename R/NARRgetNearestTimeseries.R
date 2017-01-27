#' Finds timeseries of NARR data nearest to specified location
#'
#' @description Reads a NARR NetCDF file, which contains all of the 3-hour values of a single variable for a single year, and extracts the values for a single location. Note that the values are \emph{NOT} quality controlled - negative precipitation values are possible.
#' @param ncdfFile Name of the NetCDF file containing ERA data.
#' @param varName Required. Name of the NetCDF variable to extract.
#' @param pointLon Required. Decimal longitude of desired location. Note that the NetCDF longitude is \emph{East}, your value should be negative.
#' @param pointLat  Required. Decimal latitude of desired location.
#' @param timezone Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}, as it will list the actual latitude and longitude of the ERA tile.
#'@param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If unsuccessful, returns \code{FALSE}. If successful, returns a standard \pkg{CRHMr} dataframe containing the datetime and the extracted data, which are unpacked (i.e. the NetCDF multiplier and offset have been applied). If \code{quiet=FALSE}, then the NARR location used and its distance to the specified location will be displayed. 
#'
#' @return If unsuccessful, returns \code{FALSE}. If successful, returns a \pkg{CRHMr} data frame of the specified variable. Note that no unit conversions are done, so the variable will be in the NARR units.
#' @export
#'
#' @examples \dontrun{
#' p1979 <- NARRgetNearestTimeseries('acpcp.1979.nc', varName='acpcp', pointLon = -113, 
#' pointLat = 52, timezone = 'CST', quiet=FALSE)}
NARRgetNearestTimeseries <- function(ncdfFile, varName, pointLon, pointLat, 
                               timezone='', quiet=TRUE, logfile=''){
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
  
  # project lat and lon to x and y
  projection = "+proj=lcc +lat_1=50 +lat_2=50 +lat_0=50 +lon_0=-107.0 +x_0=5632642.22547 +y_0=4612545.65137"
  xy <- list(pointLon,pointLat)
  projected <- proj4::project(xy, proj=projection)
  point_x <- projected$x
  point_y <- projected$y
  
  # calculate hour offset between timezone and GMT
  datetime <- '2000-01-01 01:00'
  s1 <- as.POSIXct(datetime, format='%Y-%m-%d %H:%M', tz='UTC')
  s2 <- as.POSIXct(datetime, format='%Y-%m-%d %H:%M', tz=timezone)
  houroffset <- as.numeric(difftime(s2, s1, units='hours'))
  
  # read all data from file
  nc <-  RNetCDF::open.nc(ncdfFile, write=FALSE)
  nctimes <- RNetCDF::var.get.nc(nc, variable='time', unpack=TRUE)
  data <- RNetCDF::var.get.nc(nc, variable=varName, unpack=TRUE)
  lons <- RNetCDF::var.get.nc(nc, variable='lon', unpack=TRUE)
  lats <- RNetCDF::var.get.nc(nc, variable='lat', unpack=TRUE)
  RNetCDF::close.nc(nc)
  
  # convert times to R times
  secs <- nctimes * 3600.0
  datetimes.utc <- as.POSIXct(secs, origin='1800-01-01', tz='UTC')
  datetimes <- datetimes.utc + (houroffset * 3600)
  datetimes <- lubridate::force_tz(datetimes, tzone=timezone)
  
  # try to project lat and lon
  lonlat <- list(lons, lats)
  projected <- proj4::project(lonlat, proj=projection)
  x_vals <- projected$x
  y_vals <- projected$y
  
  delta_x <- x_vals - point_x
  delta_y <- y_vals - point_y
  
  distance <- (delta_x^2 + delta_y^2) ^ 0.5
  min_dist <- min(distance)
  min_loc <- which(distance <= min_dist)
  
  lon <- lons[min_loc]
  lat <- lats[min_loc]  
  
  if (!quiet){
    cat('Closest location at: lon = ', lon, ' lat = ', lat, ', distance = ', 
        formatC(min_dist/1000, digits=1, format='f'), 'km\n', sep='')
  }
  
  # find location in data matrix
  col_num <- ceiling(min_loc /nrow(lons))
  row_num <- min_loc - ((col_num-1) * nrow(lons))

  # test locations
  # lon_test <- lons[row_num, col_num]
  # lat_test <- lats[row_num, col_num]
  
  # extract values
  values <- data[row_num, col_num ,]
  
  # create dataframe
  df <- data.frame(datetimes, values)
  names(df) <- c('datetime', varName)
  
  # output info to screen (if req'd) and write to log file
  file.info <- CRHMr::CRHM_summary(df)
  if (!quiet)
    print(file.info)
  
  comment <- paste('NARRgetNearestTimeseries ncdfFile:',ncdfFile,
                   ' varName:', varName,
                   ' pointLon:', pointLon,
                   ' pointLat:', pointLat,
                   ' timezone:', timezone, sep='')
  result <- CRHMr::logAction(comment, logfile)
  
  if(result)
    return(df)
  else
    return(result)
  
}