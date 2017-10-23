#' Reads a time series for the nearest point from undadjusted CanRCM4
#' @description This function extracts a time series of hourly values from a NetCDF file of CanRCM4 data, which has not been bias-corrected. The values are stored at a spatial resolution of 0.125 degees, so the closest point to the specified location will be used. Each NetCDF file contains a single variable. Typically the first letters of the file name designate the variable. The variables are converted to values, and variable names, appropriate for CRHM when they are extracted.
#' \tabular{lllr}{
#'\bold{NetCDF parameter} \tab \bold{netCDF units} \tab \bold{CRHM Variable} \tab \bold{CRHM units}\cr
#'pr - hourly precip\tab mm/s \tab p \tab mm\cr
#'tas - surface air temp. \tab K \tab t \tab  \eqn{^\circ}{}C\cr
#'huss - specific humidity \tab dimensionless \tab qair \tab dimensionless \cr
#'sfcWind - surface (10m) wind speed \tab m/s \tab u10 \tab m/s \cr
#'ps - surface pressure \tab Pa \tab ps \tab Pa \cr
#'rsds - incoming SW radiation \tab W/\eqn{^2}{^2}  \tab Qsi \tab W/\eqn{^2}{^2} \cr
#'rlds - incoming LW radiation \tab W/\eqn{^2}{^2}  \tab Qli \tab W/\eqn{^2}{^2} \cr
#' }
#' @param netCDFfile Required. The name of a NetCDF file containing a single variable.
#' @param longitude Required. The longitude of the point being sought. Valid values appear to be between -90 and -142, but the input value is \emph{not} checked for validity, in case the model extent changes.
#' @param latitude Required. The latitude of the point being sought. Valid values appear to be between 45 and 75, but the input value is \emph{not} checked for validity, in case the model extent changes.
#' @param startDate Optional. Beginning date of data to be extracted. A string formatted as "yyyy-mm-dd". The default value of \option{1979-01-01} is the beginning of the data.
#' @param endDate Optional. Beginning date of data to be extracted. A string formatted as "yyyy-mm-dd". The default value of \option{2100-12-01} is the end of the data.
#' @param timezone Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.
#'@param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#'
#' @return If successful, returns a standard \pkg{CRHMr} data frame containing the datetime and the variable. If unsuccessful, returns the value \code{FALSE}.
#' @export
#' @seealso \code{\link{qa2ea}}
#'
#' @examples \dontrun{ f <- "CanRCM4_r8i2p1r1_hist+fut_1979_2100_huss.nc4"
#' r <- RCM4UnadjustedGetNearestTimeseries(f, longitude=-115.27, 
#' latitude=52.03,  startDate = "1980-01-01", endDate="1980-12-31")
#' }
CanRCM4unadjustedGetNearestTimeseries <- function(netCDFfile = "", 
                      longitude = 0, latitude = 0, 
                      startDate = "1979-01-01",
                      endDate = "2100-12-01",
                      timezone = "Etc/GMT+7", 
                      logfile = "") {
  
  # check parameters
  if (netCDFfile == '') {
    cat('Error: missing input netCDF file name\n')
    return(FALSE)
  }
  
  if (latitude == 0) {
    cat('Error: missing latitude\n')
    return(FALSE)
  }
  
  if (longitude == 0) {
    cat('Error: missing latitude\n')
    return(FALSE)
  }
  
  # calculate hour offset between timezone and GMT
  houroffset <- CRHMr::GMToffset(timezone)
  
  # read in netCDF data
  nc <- ncdf4::nc_open(netCDFfile, write = FALSE)
  
  # get variables and locations
  lons <- ncdf4::ncvar_get( nc, attributes(nc$dim)$names[1])
  lats <- ncdf4::ncvar_get( nc, attributes(nc$dim)$names[2])
  varName <- attributes(nc$var)$names[1]
  
  # times are hours since 1979-1-1 00:00:00
  nc_times <- ncdf4::ncvar_get(nc, attributes(nc$dim)$names[3])
  nc_times <- as.numeric(nc_times)

  # find location - closest 1/8 degree
  int_lat <- floor(latitude)
  decimal_lat <- latitude - int_lat 
  lat_frac <- round(decimal_lat / 0.125)
  closest_lat <- int_lat + lat_frac
  
  int_lon <- floor(longitude)
  decimal_lon <- longitude - int_lon
  lon_frac <- round(decimal_lon / 0.125)
  closest_lon <- int_lon + lon_frac
  
  latNum <- which(lats == closest_lat)
  lonNum <- which(lons == closest_lon)
  

  # create empty obs dataframe
  startDate <- as.Date(startDate, format = "%Y-%m-%d")
  obsStartDate <- format(startDate - 1, format = "%Y-%m-%d")
  
  # convert units and set NA values
  if (varName == "pr") {                     # precip
    CRHMvariable <- "p"
  } else if (varName == "tas") {             # surface air temp
    CRHMvariable <- "t"
  } else if (varName == "huss") {            # specific humidity (dimensionless)
    CRHMvariable <- "qair"
  } else if (varName == "sfcWind") {         # 10m wind speed (m/s)
    CRHMvariable <- "u10"
  } else if (varName == "ps") {              # surface pressure in Pa
    CRHMvariable <- "ps"
  }  else if (varName == "rsds") {           # incoming SW (W/m2)
    CRHMvariable <- "Qli"
  } else if (varName == "rlds") {            # incoming LW (W/m2)
    CRHMvariable <- "Qsi"
  } else { 
    CRHMvariable <- "unknown"
  }
  
  
  empty <- CRHMr::createObsDataframe(obsStartDate, endDate, timestep = 1,
                                     variables = CRHMvariable, timezone = "UTC",
                                     logfile = logfile)
  
  numVals <- nrow(empty)
  
  # get starting point of data
  firstDatetime <- empty$datetime[1]
  ncdf_first_datetime <- as.POSIXct("1979-1-1 00:00", format = "%Y-%m-%d %H:%M", tz = "UTC")
  startPoint <- as.numeric(difftime(firstDatetime, ncdf_first_datetime, units = "hours"))
  
  # extract values 
  data <-  ncdf4::ncvar_get(nc, varName, start = c(lonNum, latNum, startPoint), 
                       count = c(1, 1, numVals)) 
  ncdf4::nc_close(nc)
  
  # convert units and set NA values
  if (varName == "pr") {                     # precip
    data[data > 1e19] <- NA_real_
    data <- data * 3600                      # convert mm/s to mm
  } else if (varName == "tas") {             # surface air temp
    data[data > 1e19] <- NA_real_
    data <- data - 273.15                    # convert K to C
  } else if (varName == "huss") {            # specific humidity (dimensionless)
    data[data > 1e19] <- NA_real_
  } else if (varName == "sfcWind") {         # 10m wind speed (m/s)
    data[data > 1e19] <- NA_real_
  } else if (varName == "ps") {              # surface pressure in Pa
    data[data > 1e19] <- NA_real_
  }  else if (varName == "rsds") {           # incoming SW (W/m2)
    data[data > 1e19] <- NA_real_
  } else if (varName == "rlds") {            # incoming LW (W/m2)
    data[data > 1e19] <- NA_real_
  } else { 
    data[data > 1e19] <- NA_real_
  }
  
  empty[,2] <- data
  filled <- empty
  rm(empty)

  # force timezone
  filled$datetime <- filled$datetime - (houroffset * 3600)
  filled$datetime <- lubridate::force_tz(filled$datetime, tzone = timezone)
  return(filled)
}