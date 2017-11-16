#' Extracts CRHM obs variables from a WRF NetCDF at a single specified location
#' @description This function is intended to work with NetCDF created by the research group of Dr. Yanping Li of Global Water Futures \url{gwf.usask.ca}. This function creates a \pkg{CRHMr} data frame from a WRF NetCDF specified by their location within the file. This function will work even if there are only values for a single location inside the NetCDF.
#' @param NetCDF Required. The name of the NetCDF containing the WRF data.
#' @param obsfile Optional. If specified the values are written to the \code{obs} file.
#' @param startDatetime Optional. Beginning datetime of data to be extracted. A string formatted as "yyyy-mm-dd hh:mm". The default value is \option{2000-10-01 00:00}.
#' @param endDatetime Optional. Endning datetime of data to be extracted. A string formatted as "yyyy-mm-dd hh:mm". The default value is \option{2100-12-01 23:00}.
#' @param west_east_loc Optional. The NetCDF west_east location as an integer. The default value of \code{1} specifies the westernmost location.
#' @param south_north_loc Optional. The NetCDF south_north location as an integer. The default value of \code{1} specifies the southernernmost location.
#' @param airPressure Optional. The mean air pressure in Pa. If set to zero (the default) then the air pressure within the NetCDF file will be used for converting the specific humidity to relative humidity. Because the air pressure may be incorrect or missing, you may wish to specify a mean air pressure to be used for the conversion. Note that if there is no air pressure within the NetCDF and you do not specify an air pressure value, this will trigger an error message, and the RH cannot be calculated.
#' @param timezone Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.
#' @param trim Optional. If \code{TRUE} (the default) then the final output will be trimmed to CRHM day boundaries.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If successful, returns either \code{TRUE} (if an obs file is specified) or a \pkg{CRHMr} obs data frame (if no obs file is specified). The data frame will contain the standard CRHM variables: t, p, u, either rh or ea, qsi, and qli, depending on their presence in the NetCDF. Note that the precipitation is the accumulated value, because it can contain negative spikes and resets. To deaccumulate the precipitation, you will have to use the \pkg{CRHMr} weighing gauge functions, such as \code{weighingGauge1}, to fill gaps, \code{weighingGauge2} or \code{weighingGauge5} to remove resets. You can use \code{weighingGaugeInterval} to deaccumulate the values. Note that because these functions use thresholds to define missing values, spike and resets, they need to be used interatively. You can use the \code{weighingGaugePlot} to see how well the data are being processed. You may also need to remove too-large precipitation events. If unsucessful, returns \code{FALSE}.
#' @seealso \code{\link{WRFnearest2obs}}
#' @author Kevin Shook
#' @export
#'
#' @examples \dontrun{f <- "wrf2d_tquvr.nc"
#' r <- WRFbyloc2Obs(f, startDatetime = "1980-01-01 00:00", endDatetime = "1980-12-31 23:00")}
#' 
WRFbyloc2obs <- function(NetCDF = "", obsfile = "", 
                           startDatetime = "2000-10-01 00:00",
                           endDatetime = "2100-12-01 23:00",
                           west_east_loc = 1,
                           south_north_loc = 1,
                           airPressure = 0,
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


  # calculate hour offset between timezone and GMT
  houroffset <- CRHMr::GMToffset(timezone)
  
  if (timezone == "") {
    cat('Error: missing time zone\n')
    return(FALSE)
  }
  
  # read in netCDF data
  nc <- ncdf4::nc_open(NetCDF, write = FALSE)
  
  # create dataframe of x and y values
  x_count <- nc$dim$west_east$len
  y_count <- nc$dim$south_north$len
  
  if (west_east_loc > x_count) {
    cat('Error: value of west_east_loc, ', west_east_loc, ', is greater than\n', sep = "")
    cat("number of locations, ", x_count, ", in file ", NetCDF, "\n", sep = "")
    return(FALSE)
  }
  
  if (south_north_loc > y_count) {
    cat('Error: value of south_north_loc, ', south_north_loc, ', is greater than\n', sep = "")
    cat("number of locations, ", y_count, ", in file ", NetCDF, "\n", sep = "")
    return(FALSE)
  }
  
  
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

  # get variable names
  num_vars <- as.numeric(nc$nvars)
   for (i in 1:num_vars) {
     var_names[i] <- nc$var[[i]]$name
   }
  
  nctimes <- nc$dim$Time$vals
  dt <- nctimes[2] - nctimes[1]   # minutes
  
  timeLocs <- seq(1, length(nctimes))
  data_datetime <- timeLocs * (dt * 3600) + data_start_datetime  
  timeDF <- data.frame(timeLocs, data_datetime)
  
  selected_times <- timeDF[(timeDF$data_datetime >= selected_start_datetime) &
                             (timeDF$data_datetime <= selected_end_datetime),] 
  
  startPoint <- selected_times[1,1]
  numVals <- nrow(selected_times)
  
  if (is.na(startPoint)) {
    cat('Error: specified start or end datetimes not found in NetCDF\n')
    return(FALSE)
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
  
  if (length(grep("GLW", var_names)) > 0) {
    GLW_exists <- TRUE
  } else {
    GLW_exists <- FALSE
  }
  
  if (Q2_exists) {
    Q2 <- as.numeric(ncdf4::ncvar_get(nc, varid = "Q2", 
                                      start = c(west_east_loc, south_north_loc, 
                                                startPoint), 
                                      count = c(1, 1, numVals)))
  }
  
  if (T2_exists) {
    T2 <- as.numeric(ncdf4::ncvar_get(nc, varid = "T2", 
                                      start = c(west_east_loc, south_north_loc,
                                                startPoint), 
                                      count = c(1, 1, numVals)))
  }
  
  if (wind_exists) {
    U10 <- as.numeric(ncdf4::ncvar_get(nc, varid = "U10", 
                                       start = c(west_east_loc, south_north_loc,
                                                 startPoint), 
                                       count = c(1, 1, numVals)))
    
    V10 <- as.numeric(ncdf4::ncvar_get(nc, varid = "V10",  
                                       start = c(west_east_loc, south_north_loc,
                                                 startPoint), 
                                       count = c(1, 1, numVals)))
  }
  
  if (precip_exists) {
    I_RAINNC <- as.numeric(ncdf4::ncvar_get(nc, varid = "I_RAINNC",
                                        start = c(west_east_loc, south_north_loc,
                                                  startPoint), 
                                    count = c(1, 1, numVals)))

    RAINNC <- as.numeric(ncdf4::ncvar_get(nc, varid = "RAINNC",
                                       start = c(west_east_loc, south_north_loc,
                                                startPoint), 
                                      count = c(1, 1, numVals)))
  }

  if (SWDOWN_exists) {
    QSI <- as.numeric(ncdf4::ncvar_get(nc, varid = "SWDOWN", 
                                       start = c(west_east_loc, south_north_loc,
                                                startPoint), 
                                      count = c(1, 1, numVals)))
  }
  
  if (GLW_exists) {
    QLI <- as.numeric(ncdf4::ncvar_get(nc, varid = "GLW", 
                                       start = c(west_east_loc, south_north_loc,
                                                 startPoint), 
                                       count = c(1, 1, numVals)))
  }
  
  if (PSFC_exists) {
    PSFC <- as.numeric(ncdf4::ncvar_get(nc, varid = "PSFC", 
                                       start = c(west_east_loc, south_north_loc,
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

    # remove resets and drop-outs
    # accum <- data.frame(df1$datetime, cummax(p))
    # names(accum) <- c("datetime", "p")
    # deaccum <- CRHMr::weighingGaugeInterval(accum)
    df1$p <- p
  }
  
  if (T2_exists) {
    df1$t <- T2 - 273.15
  }

  if (Q2_exists & T2_exists) {
    if (!quiet) {
      cat('Air temperature is present, so will calculate RH\n')
    }
    if ((airPressure <= 0) & PSFC_exists) {
      # convert Q to RH 
      df1$rh <- CRHMr::qair2rh(Q2, df1$t, PSFC) * 100
    } else if (airPressure > 0 ) {
     df1$rh <- CRHMr::qair2rh(Q2, df1$t, airPressure) * 100
   }
   else (
     cat("Error - cannot calculate RH without air pressure\n")
   )
    
  }
  
  if (Q2_exists & !T2_exists) {
    if (!quiet) {
      cat('No air temperature, so cannot calculate RH, so will calculate ea\n')
    }
    if ((airPressure <= 0) & PSFC_exists) {
      # convert Q to ea 
      df1$ea <- CRHMr::qair2ea(Q2, PSFC )
    } else if (airPressure > 0 ) {
      df1$ea <- CRHMr::qair2ea(Q2, airPressure)
    }
    else {
      cat("Error - cannot calculate ea without air pressure\n")
    }
  }

  if (SWDOWN_exists) {
    df1$qsi <- QSI
  }  
  
  if (GLW_exists) {
    df1$qli <- QLI
  }  
  

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
                                      
  comment <- "WRF hourly data"
  
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
