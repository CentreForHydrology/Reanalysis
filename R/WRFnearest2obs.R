#' Extracts CRHM obs variables from a WRF NetCDF closest to specified point
#' @description This function is intended to work with NetCDF created by the research group of Dr. Yanping Li of Global Water Futures \url{gwf.usask.ca}. This function creates a \pkg{CRHMr} data frame from a WRF NetCDF specified by their longitude and latitude, extracting values from the nearest location.
#' @param NetCDF Required. The name of the netCDF file containing the WRF data.
#' @param obsfile Optional. If specified the values are written to the \code{obs} file.
#' @param longitude Required. The longitude of the point being sought. The input value is \emph{not} checked for validity, in case the model extent changes.
#' @param latitude Required. The latitude of the point being sought. The input value is \emph{not} checked for validity, in case the model extent changes.
#' @param airPressure Optional. The mean air pressure in Pa. If set to zero (the default) then the air pressure within the NetCDf file will be used for converting the specific humidity to relative humidity. Because the air pressure may be incorrect or missing, you may wish to specify a mean air pressure to be used for the conversion. Note that if there is no air pressure within the NetCDF and you do not specify an air pressure value, this will trigger an error message, and the RH cannot be calculated. 
#' @param startDatetime Optional. Beginning date of data to be extracted. A string formatted as "yyyy-mm-dd hh:mm". The default value is \option{2000-10-01 00:00}.
#' @param endDatetime Optional. Ending date of data to be extracted. A string formatted as "yyyy-mm-dd hh:mm". The default value is \option{2100-12-01 23:00}.
#' @param timezone Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.
#' @param trim Optional. If \code{TRUE} (the default) then the final output will be trimmed to CRHM day boundaries.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If successful, returns the value \code{TRUE} and (optionaly) writes the specified .obs file. If unsuccessful, returns the value \code{FALSE}.
#' @author Kevin Shook
#' @seealso \code{\link{WRFbyloc2obs}}
#' @export
#'
#' @examples \dontrun{ f <- "wrf2d_tquvr.nc"
#' r <- WRFnearest2obs(f, longitude=-115.27, 
#' latitude=52.03, startDatetime = "1980-01-01 00:00", endDatetime="1980-12-31 23:00")
#' }
WRFnearest2obs <- function(NetCDF = "", obsfile = "", 
                           longitude = 0,
                           latitude = 0,
                           airPressure = 0,
                           startDatetime = "2000-10-01 00:00",
                           endDatetime = "2100-12-01 23:00",
                           timezone = "", 
                           trim = TRUE, 
                           quiet = TRUE, 
                           logfile = "") {
  # create global variables
  
  var_names <- c(0)

  
  # check parameters
  if (NetCDF == '') {
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
  
  if (timezone == "") {
    cat('Error: missing time zone\n')
    return(FALSE)
  }
  
  # read in netCDF data
  nc <- ncdf4::nc_open(NetCDF, write = FALSE)
  
  # get longs and lats - if they exist
  num_vars <- as.numeric(nc$nvars)
  for (i in 1:num_vars) {
    var_names[i] <- nc$var[[i]]$name
  }
  if (length(grep("XLONG", var_names)) == 0) {
    cat("NetCDF does not contain longitudes or latitudes \n")
    return(FALSE)
  }
  lons <- ncdf4::ncvar_get(nc, varid = "XLONG")
  lats <- ncdf4::ncvar_get(nc, varid = "XLAT")
  
  ncdf4::nc_close(nc)
  
  dist <- ((lons - longitude) ^ 2 + (lats - latitude) ^ 2) ^ 0.5
  mindist <- min(dist)
  minloc <- which.min(dist)
  minlat <- lats[minloc]
  minlon <- lons[minloc]

  # get X and Y of closest point
  locs <- which(dist == mindist, arr.ind = TRUE)
  
  # transpose locations, because locations are retrieved transposed
  west_east_loc <- locs[1]
  south_north_loc <- locs[2]
  
  # extract data
  values <- WRFbyloc2obs(NetCDF = NetCDF,
                           obsfile = "", 
                           startDatetime = startDatetime,
                           endDatetime = endDatetime,
                           west_east_loc = west_east_loc,
                           south_north_loc = south_north_loc,
                           airPressure = airPressure,
                           timezone = timezone, 
                           trim = trim, 
                           quiet = quiet, 
                           logfile = logfile) 
  
  comment <- paste("WRF hourly data. Specified lon, lat:", longitude,",", 
                   latitude ,"Extracted lon,lat", 
                   minlon, ",", minlat, sep = "")
  
  if (obsfile != "") {
    result <- CRHMr::writeObsFile(values, obsfile, comment = comment, quiet, 
                                  logfile)
    if (!result)
      return(result)
    else
      return(values)
  }
   else {
     return(values)
   }
}
