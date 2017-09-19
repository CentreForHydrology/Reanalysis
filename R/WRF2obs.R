#' Converts a WRF netcdf file to obs
#'
#' @description Creates a \pkg{CRHMr} data frame of all values from a WRF netCDF file containing t, RH, p, and wind vectors for a \emph{single} location.
#' @param netCDFfile Required. The name of the netCDF file containing the WRF data.
#' @param obsfile Optional. If specified the values are written to the \code{obs} file.
#' @param startDateTime Optional. Date and time of first interval. Format is "yyyy-mm-dd hh:mm".
#' @param timezone Optional. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. You can use \option{Etc/GMT+6} or \option{Etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings.
#' @param trim Optional. If \code{TRUE} (the default) then the final output will be trimmed to CRHM day boundaries.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If successful, returns the value \code{TRUE} and (optionaly) writes the specified .obs file. If unsuccessful, returns the value \code{FALSE}.
#'
#' @export
#'
#' @examples \dontrun{w <- WRF2Obs("wrf2d_tquvr.nc")}
WRF2Obs <- function(netCDFfile = "", obsfile = "", startDatetime = "2000-10-01 00:00",
                    timezone="Etc/GMT+7", trim = TRUE, quiet = TRUE, logfile = "") {
  # check parameters
  if (netCDFfile == '') {
    cat('Error: missing input netCDF file name\n')
    return(FALSE)
  }
  
  # calculate hour offset between timezone and GMT
  houroffset <- CRHMr::GMToffset(timezone)

  # read in netCDF data
  nc <- RNetCDF::open.nc(netCDFfile, write = FALSE)
  
  # get variables and locations
  nctimes <- RNetCDF::var.get.nc(nc, variable = 'XTIME', unpack = TRUE)
  Q2 <- RNetCDF::var.get.nc(nc, variable = "Q2", unpack = TRUE)
  T2 <- RNetCDF::var.get.nc(nc, variable = "T2", unpack = TRUE)
  U10 <- RNetCDF::var.get.nc(nc, variable = "U10", unpack = TRUE)
  V10 <- RNetCDF::var.get.nc(nc, variable = "V10", unpack = TRUE)
  I_RAINNC <- RNetCDF::var.get.nc(nc, variable = "I_RAINNC", unpack = TRUE)
  RAINNC <- RNetCDF::var.get.nc(nc, variable = "RAINNC", unpack = TRUE)  
  
  RNetCDF::close.nc(nc)
  
  # creat data frame and convert variables
  u10 <- (U10 ^ 2 + V10 ^ 2) ^ 0.5
  p <-  I_RAINNC * 100 + RAINNC
  
  df1 <- data.frame(nctimes)
  names(df1) <- "datetime"
  
  startTime <- as.POSIXct(startDatetime, format = "%Y-%m-%d %H:%M", 
                          tz = "UTC")
  df1$datetime <- (df1$datetime * 60) + startTime
  
  # add data, doing unit conversions
  df1$t <- T2 - 273.15
  
  # convert Q to RH
  df1$rh <- CRHMr::qair2rh(Q2, df1$t) * 100
  df1$u <- u10
  df1$p <- p
  
  # remove resets
  accum <- data.frame(df1$datetime, cummax(df1[,5]))
  names(accum) <- c("datetime", "p")
  deaccum <- CRHMr::weighingGaugeInterval(accum)
  df1$p <- deaccum$p_interval
  
  # figure out the start and end dates
  startDate <- as.Date(startDatetime)
  endDate <- as.Date(df1$datetime[nrow(df1)])
  startDate <- format(startDate, format = "%Y-%m-%d")
  endDate <- format(endDate, format = "%Y-%m-%d")
  
  # create the obs data frame
  df2 <- CRHMr::createObsDataframe(start.date = startDate, end.date = endDate,
                                   timezone = "UTC", logfile = logfile)
  
  # merge the WRF data into the obs dataframe
  merged <- merge(df1, df2, by = "datetime", all = TRUE)
  
  # delete extra columns
  merged <- merged[,c(1,2,3,4,5)]
  names(merged) <- c("datetime", "t.1", "rh.1", "u.1", "p.1")
  
  # convert to required time zone
  datetime = merged$datetime + (houroffset * 3600)
  
  # force timezone
  datetime <- lubridate::force_tz(datetime, tzone = timezone)
  merged$datetime <- datetime
  
  
  # trim if required
  if (trim)
    df <- CRHMr::trimObs(merged, quiet, logfile)
  
  if (obsfile != "") {
    result <- CRHMr::writeObsFile(df, obsfile, comment = "WRF hourly data", quiet, logfile)
    if (!result)
      return(result)
    else
      return(df)
  } else
    return(df)
  
}