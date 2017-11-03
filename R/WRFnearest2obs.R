#' Gets CRHM obs veriables closest to specified point
#' @description Creates a \pkg{CRHMr} data frame of all values from a WRF netCDF file for a specified location.
#' @param netCDFfile Required. The name of the netCDF file containing the WRF data.
#' @param obsfile Optional. If specified the values are written to the \code{obs} file.
#' @param longitude Required. The longitude of the point being sought. The input value is \emph{not} checked for validity, in case the model extent changes.
#' @param latitude Required. The latitude of the point being sought. The input value is \emph{not} checked for validity, in case the model extent changes.
#' @param startDatetime Optional. Beginning date of data to be extracted. A string formatted as "yyyy-mm-dd hh:mm". The default value is \option{2000-10-01 00:00}.
#' @param endDatetime Optional. Ending date of data to be extracted. A string formatted as "yyyy-mm-dd hh:mm". The default value is \option{2100-12-01 23:00}.
#' @param timezone Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.
#' @param trim Optional. If \code{TRUE} (the default) then the final output will be trimmed to CRHM day boundaries.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If successful, returns the value \code{TRUE} and (optionaly) writes the specified .obs file. If unsuccessful, returns the value \code{FALSE}.
#'
#' @return If successfl, returns either \code{TRUE} (if an obs file is specified) or a \pkg{CRHMr} obs data frame (if no obs file is specified). If unsucessful, retuns \code{FALSE}. 
#' @export
#'
#' @examples \dontrun{ f <- "wrf2d_tquvr.nc"
#' r <- WRFnearest2obs(f, longitude=-115.27, 
#' latitude=52.03,  startDatetime = "1980-01-01 00:00", endDate="1980-12-31 23:00")
#' }
WRFnearest2obs <- function(netCDFfile = "", obsfile = "", 
                           longitude = 0,
                           latitude = 0,
                           startDatetime = "2000-10-01 00:00",
                           endDatetime = "2100-12-01 23:00",
                           timezone = "", 
                           trim = TRUE, 
                           quiet = TRUE, 
                           logfile = "") {
  # create global variables
  
  var_names <- c(0)

  
  # check parameters
  if (netCDFfile == '') {
    cat('Error: missing input netCDF file name\n')
    return(FALSE)
  }


  # calculate hour offset between timezone and GMT
  houroffset <- CRHMr::GMToffset(timezone)
  
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
  nc <- ncdf4::nc_open(netCDFfile, write = FALSE)
  
  
  # create dataframe of x and y values
  x_count <- nc$dim$west_east$len
  y_count <- nc$dim$south_north$len
  
  # extract projection info
  projType <- ncdf4::ncatt_get(nc, varid = 0, attname = "MAP_PROJ_CHAR")$value
  if (projType != "Lambert Conformal") {
    cat("Error: projection is ", projType, ". Don't know how to reproject\n", sep = "")
    return(FALSE)
  }
  
  dx <- ncdf4::ncatt_get(nc, varid = 0, attname = "DX")$value
  dy <- ncdf4::ncatt_get(nc, varid = 0, attname = "DY")$value
  cen_lat <- ncdf4::ncatt_get(nc, varid = 0, attname = "CEN_LAT")$value
  cen_lon <- ncdf4::ncatt_get(nc, varid = 0, attname = "CEN_LON")$value 
  truelat1 <- ncdf4::ncatt_get(nc, varid = 0, attname = "TRUELAT1")$value
  truelat2 <- ncdf4::ncatt_get(nc, varid = 0, attname = "TRUELAT2")$value
  stand_lon <- ncdf4::ncatt_get(nc, varid = 0, attname = "STAND_LON")$value
  
  # assemble proj4 string
  proj4string <- paste("+proj=lcc +lon 0=",
                       stand_lon,"e +lat 1=",truelat1,
                       "n +lat 2=",truelat2,"n", sep = "")
  
  # get locations of center point and specified point
  xy1 <- list(cen_lon, cen_lat)
  projected1 <- proj4::project(xy1, proj = proj4string)
  center_x <- projected1$x
  center_y <- projected1$y
  
  xy2 <- list(longitude, latitude)
  projected2 <- proj4::project(xy2, proj = proj4string)
  specified_x <- projected2$x
  specified_y <- projected2$y
  
  
  # assemble array of locations
  maxx <- floor(x_count / 2)
  minx <- -maxx
  maxy <- floor(y_count / 2)
  miny <- -maxy
  
  xseq <- seq(minx, maxx)
  yseq <- seq(miny, maxy)
  
  xlocs <- center_x + (xseq * dx)
  ylocs <- center_y + (yseq * dy)
  
  all_x <-  rep(xlocs, each = y_count)
  all_x_locs <- rep(seq(1:x_count), each = y_count)
  
  all_y <-  rep(ylocs, each = x_count)
  all_y_locs <- rep(seq(1:y_count), each = x_count)
  
  all_locs <- data.frame(all_x_locs, all_x, all_y_locs, all_y)
  
  # now get distance
  all_locs$x_distance <- (all_locs$all_x - specified_x) ^ 2
  all_locs$y_distance <- (all_locs$all_y - specified_y) ^ 2
  all_locs$distance <- (all_locs$x_distance + all_locs$y_distance) ^ 0.5
  
  # get minimum
  minpoint <- which.min(all_locs$distance)
  min_distance <- min(all_locs$distance) / 1000
  closest_x <- all_locs$all_x[minpoint]
  closest_y <- all_locs$all_y[minpoint]
  closest_x_loc <- all_locs$all_x_locs[minpoint]
  closest_y_loc <- all_locs$all_y_locs[minpoint]
  
  # figure out time
  data_start_date <- ncdf4::ncatt_get(nc, varid = 0, 
                                      attname = "START_DATE")$value
  data_start_datetime <- as.POSIXct(data_start_date, 
                                    format = "%Y-%m-%d_%H:%M:%S",
                                    tz = "UTC")
  
  selected_start_datetime <- as.POSIXct(startDatetime, 
                                    format = "%Y-%m-%d %H:%M",
                                    TZ = "UTC")
  selected_end_datetime <- as.POSIXct(endDatetime, 
                                        format = "%Y-%m-%d %H:%M",
                                        TZ = "UTC")
  
  
  
  # now extract values
  nctimes <- as.numeric(ncdf4::ncvar_get(nc, varid = 'XTIME', raw_datavals = FALSE))
  dt <- nctimes[2] - nctimes[1]   # minutes

  timeLocs <- seq(1, length(nctimes))
  data_datetime <- timeLocs * (dt * 60) + data_start_datetime  
  timeDF <- data.frame(timeLocs, data_datetime)

  selected_times <- timeDF[(timeDF$data_datetime >= selected_start_datetime) &
                          (timeDF$data_datetime <= selected_end_datetime),] 
  
  startPoint <- selected_times[1,1]
  numVals <- nrow(selected_times)
  
  # get variable names
  num_vars <- as.numeric(nc$nvars)
   for (i in 1:num_vars) {
     var_names[i] <- nc$var[[i]]$name
   }
  
  # get values, if they exist
  if (length(grep("Q2", var_names)) > 0) {
    Q2_exists <- TRUE
  } else {
    Q2_exists <- FALSE
  }
  
  if (length(grep("T2", var_names)) > 0) {
    T2_exists <- TRUE
  } else {
    T2_exists <- FALSE
  }
  
  if (length(grep("U10", var_names)) > 0) {
    U10_exists <- TRUE
  } else {
    U10_exists <- FALSE
  }
  
  if (length(grep("V10", var_names)) > 0) {
    V10_exists <- TRUE
  } else {
    V10_exists <- FALSE
  }
  
  if ( U10_exists & V10_exists) {
    wind_exists <- TRUE
  } else {
    wind_exists <- FALSE
  }
    
  
  if (length(grep("I_RAINNC", var_names)) > 0) {
    I_RAINNC_exists <- TRUE
  } else {
    I_RAINNC_exists <- FALSE
  }
  
  if (length(grep("RAINNC", var_names)) > 1) {
    RAINNC_exists <- TRUE
  } else {
    RAINNC_exists <- FALSE
  }
  
  if (RAINNC_exists & I_RAINNC_exists) {
    precip_exists <- TRUE
  } else {
    precip_exists <- FALSE
  }
  
  
  if (length(grep("PSFC", var_names)) > 0) {
    PSFC_exists <- TRUE
  } else {
    PSFC_exists <- FALSE
  }
  
  if (length(grep("SWDOWN", var_names)) > 0) {
    SWDOWN_exists <- TRUE
  } else {
    SWDOWN_exists <- FALSE
  }
  
  if (Q2_exists) {
    Q2 <- as.numeric(ncdf4::ncvar_get(nc, varid = "Q2", 
                                      start = c(closest_x_loc , closest_y_loc , 
                                                startPoint), 
                                      count = c(1, 1, numVals)))
  }
  
  if (T2_exists) {
    T2 <- as.numeric(ncdf4::ncvar_get(nc, varid = "T2", 
                                      start = c(closest_x_loc , closest_y_loc , 
                                                startPoint), 
                                      count = c(1, 1, numVals)))
  }
  
  if (wind_exists) {
    U10 <- as.numeric(ncdf4::ncvar_get(nc, varid = "U10", 
                                       start = c(closest_x_loc , closest_y_loc , 
                                                 startPoint), 
                                       count = c(1, 1, numVals)))
    V10 <- as.numeric(ncdf4::ncvar_get(nc, varid = "V10",  
                                       start = c(closest_x_loc , closest_y_loc , 
                                                 startPoint), 
                                       count = c(1, 1, numVals)))
  }
  
  if (precip_exists) {
    I_RAINNC <- as.numeric(ncdf4::ncvar_get(nc, varid = "I_RAINNC",
                                    start = c(closest_x_loc , closest_y_loc , 
                                              startPoint), 
                                    count = c(1, 1, numVals)))

    RAINNC <- as.numeric(ncdf4::ncvar_get(nc, varid = "RAINNC",
                                      start = c(closest_x_loc , closest_y_loc , 
                                                startPoint), 
                                      count = c(1, 1, numVals)))
  }

  if (SWDOWN_exists) {
    QSI <- as.numeric(ncdf4::ncvar_get(nc, varid = "SWDOWN", 
                                      start = c(closest_x_loc , closest_y_loc , 
                                                startPoint), 
                                      count = c(1, 1, numVals)))
  }
  
  if (PSFC_exists) {
    PSFC <- as.numeric(ncdf4::ncvar_get(nc, varid = "PSFC", 
                                       start = c(closest_x_loc , closest_y_loc , 
                                                 startPoint), 
                                       count = c(1, 1, numVals)))
  }
  
  
  ncdf4::nc_close(nc)
  
  # create data frame
  
  df1 <- data.frame(selected_times$data_datetime)
  names(df1) <- "datetime"
  
  
  # add data, doing unit conversions
  
  if (wind_exists) {
    u10 <- (U10 ^ 2 + V10 ^ 2) ^ 0.5
    df1$u <- u10
  }
  
  if (precip_exists) {
    p <-  I_RAINNC * 100 + RAINNC

    # remove resets
    accum <- data.frame(df1$datetime, cummax(p))
    names(accum) <- c("datetime", "p")
    deaccum <- CRHMr::weighingGaugeInterval(accum)
    df1$p <- deaccum$p_interval
    
  }
  
  if (T2_exists) {
    df1$t <- T2 - 273.15
  }

  if (Q2_exists & PSFC_exists) {
    # convert Q to RH 
    df1$rh <- CRHMr::qair2rh(Q2, df1$t, PSFC * 0.01) * 100
  }

  if (SWDOWN_exists) {
    df1$qsi <- QSI
  }  
  
  var_names <- names(df1)


  # figure out the start and end dates
  startDate <- as.Date(startDatetime)
  endDate <- as.Date(df1$datetime[nrow(df1)])
  startDate <- format(startDate, format = "%Y-%m-%d")
  endDate <- format(endDate, format = "%Y-%m-%d")
  
  # create the obs data frame
  df2 <- CRHMr::createObsDataframe(start.date = startDate, end.date = endDate,
                                   timezone = "UTC", logfile = logfile)
  df2 <- data.frame(df2$datetime)
  names(df2) <- "datetime"
  
  # merge the WRF data into the obs dataframe
  merged <- merge(df2, df1, by = "datetime", all = TRUE)
  
  
  # convert to required time zone
  datetime = merged$datetime + (houroffset * 3600)
  
  # force timezone
  datetime <- lubridate::force_tz(datetime, tzone = timezone)
  merged$datetime <- datetime
  
  
  # trim if required
  if (trim)
    df <- CRHMr::trimObs(merged, quiet, logfile )
  else
    df <- merged
  
  # write comment showing actual location
  # project closest grid point to long, lat
  
  closestloc <- list(closest_x_loc, closest_y_loc)
  closest_unprojected <- proj4::project(closestloc, proj = proj4string,
                                        inverse = TRUE)
  comment <- paste("WRF hourly data. Extracted lon,lat", closest_unprojected[[1]], 
        ",", closest_unprojected[[2]], " distance from specified: ",
        min_distance, " km", sep = "")
  
  if (obsfile != "") {
    result <- CRHMr::writeObsFile(df, obsfile, comment = comment, quiet, logfile)
    if (!result)
      return(result)
    else
      return(df)
  }
   else {
     return(df)
     
   }

}
