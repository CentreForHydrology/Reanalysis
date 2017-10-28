#' Extracts all CRHM obs variables for a specified location and duration and constructs an obs file of hourly values.
#' @description This function extracts the time series of 3-hourly values from a set of NetCDF files of CanRCM4 data, which has been bias-corrected using the WFDEI 3-hour reanalysis values. The values are stored at a spatial resolution of 0.125 degees, so the closest point to the specified location will be used. Note that because the reanalysis data omit February 29, the values returned by this function will hbe interplolate from February 28 and March 1 for leap days, except for precipitation values, which are set to zero. The variables are converted to values, and variable names, appropriate for CRHM when they are extracted. The three-hourly values of temeperature, wind speed
#' @param startDate Optional. Beginning date of data to be extracted. A string formatted as "yyyy-mm-dd". The default value of \option{1979-01-01} is the beginning of the data.
#' @param endDate Optional. Beginning date of data to be extracted. A string formatted as "yyyy-mm-dd". The default value of \option{2100-12-01} is the end of the data.
#' @param longitude Required. The longitude of the point being sought. Valid values appear to be between -90 and -142, but the input value is \emph{not} checked for validity, in case the model extent changes.
#' @param latitude Required. The latitude of the point being sought. Valid values appear to be between 45 and 75, but the input value is \emph{not} checked for validity, in case the model extent changes.
#' @param sunTimeOffset Optional. The offset (in hours) is added to the solar time to convert it to local time. The default value \code{2} shifts the daily peak to occur at 2pm. This is required to downscale the Qsi values.
#' @param locationName Optional. Name for the location. This value is used to construct the name of the output file.
#' @param fileStr Optional. The name for the NetCDF .nc files containing the data, exclusive of their path and variable name.
#' @param inDir Optional. The path to the NetCDF .nc files. Default is the current directory.
#' @param outDir Optional. The path output file(s). Default is the current directory.
#' @param timezone Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.
#' @param write3hour Optional. Should the three-hour values be written to a .obs file? Default is \code{FALSE}.
#' @param quiet Optional. If \code{TRUE} then comments will be written to the screen as the input data is processed. Note that this function may take a long time to execute, so the default is the \emph{opposite} value used by most functions.
#'
#' @return If successful, writes \pkg{CRHMr} obs file(s) and returns \code{TRUE}. If unsuccessful, returns \code{FALSE}.
#' @export
#' @seealso \code{\link{CanRCM4adjustedGetNearestTimeseries}} \code{\link{qa2ea}}
#'
#' @examples \dontrun{CanRCM4AdjustedCreateHourlyObs(startDate = "1980-01-01", 
#' endDate = "1980-12-31", longitude = -101.704683, latitude = 50.845585, 
#' locationName = "SmithCreek")
#' }
CanRCM4AdjustedCreateHourlyObs <- function(startDate = "1979-01-01", 
                                           endDate = "2100-12-01",
                                           longitude = 0,
                                           latitude = 0, 
                                           sunTimeOffset = 2,
                                           locationName = "",
                                           fileStr = "_CanRCM4_hist+fut_1979_2100",
                                           inDir ="./",
                                           outDir = "./",
                                           timezone = "etc/GMT+7",
                                           write3hour = FALSE, 
                                           quiet = FALSE) {
  
  # check parameters
  
  if (latitude == 0) {
    cat('Error: missing latitude\n')
    return(FALSE)
  }
  
  if (longitude == 0) {
    cat('Error: missing longitude\n')
    return(FALSE)
  }
  
  if (startDate == "") {
    cat('Error: missing start date\n')
    return(FALSE)
  }
  
  if (endDate == "") {
    cat('Error: missing end date\n')
    return(FALSE)
  }
  
  
  startYear <- as.numeric(format(as.POSIXct(startDate, format = "%Y-%m-%d"), 
                                 format = "%Y"))
  endYear <- as.numeric(format(as.POSIXct(endDate, format = "%Y-%m-%d"), 
                               format = "%Y"))
  
  vars <- c("tas", "pr", "huss", "sfcWind", "ps", "rsds", "rlds")
  CRHMvars <- c("t", "p", "specHum", "u10", "press", "qsi", "qli")
  numVars <- length(vars)

  for (i in 1:numVars) {
    
    varName <- vars[i]
    if (!quiet) {
      cat("variable: ", varName, "\n", sep = "")
    }
    

    # extract data and write obs file
    netcdf <- paste(inDir, varName, fileStr, ".nc4", sep = "")

    
    obs <- CanRCM4adjustedGetNearestTimeseries(netCDFfile = netcdf, 
                                             longitude = longitude,
                                             latitude = latitude, 
                                             startDate = startDate, 
                                             endDate = endDate,
                                             timezone = timezone)
    
    if (i == 1)
      obs3hr <- obs
    else
      obs3hr <- cbind(obs3hr, obs[,2])
  }
  rm(obs)
  
  names(obs3hr) <- c("datetime", CRHMvars)
  
  # remove leading missing lines
  firstDays <- obs3hr[1:16,2]
  missingFirstDays <- which(is.na(firstDays))
  obs3hr <- obs3hr[-(missingFirstDays),]

  
  # infill gaps in leap years
  missingLocs <- which(is.na(obs3hr[,2]))
  
  # get values before and after gaps
  prev <- missingLocs - 8
  after <- missingLocs + 8
  
  # get infill vals
  numcols <- ncol(obs3hr)
  
  for (col in 2:(numcols)) {
    obs3hr[missingLocs, col] <- (obs3hr[prev, col] + obs3hr[after, col]) / 2
  }  
  
  # convert specific humidity to ea
  obs3hr$ea <- qa2ea(obs3hr$specH, obs3hr$press)
  
  # select variables and output
  obs3hr <- obs3hr[,c("datetime", "t", 
                "p",  "u10", "ea", "qsi", "qli")]
  
  if (write3hour) {
    obs3hrFile <- paste(outDir, locationName,"_",startYear,"-",endYear,
                        "_adjusted_3hr_CanRCM4.obs", sep = "")
    if (!quiet) {
      cat("writing 3-hour obs file to ", obs3hrFile, "\n", sep = "")
    }
    CRHMr::writeObsFile(obs3hr, obs3hrFile, quiet = quiet)    
  }

  
  # convert to 1-hour values
  
  if (!quiet) {
    cat("converting 3-hourly to hourly\n")
  }
  
  firstDate <- format(obs3hr$datetime[1], format = "%Y-%m-%d")
  numObs <- nrow(obs3hr)
  lastDate <- format(obs3hr$datetime[numObs], format = "%Y-%m-%d")
  
  hourlyObs <- CRHMr::createObsDataframe(start.date = firstDate, 
                                         end.date = lastDate,
                                  timestep = 1, 
                                  variables = c("t", "p",  "u10", "ea", "qsi", "qli"), 
                                  timezone = timezone)
  
  # create reps of each variable
  colnums <- c(2:7)
  for (col in colnums) {
    vals <- obs3hr[,col]
    reps <- rep(vals, each = 3)
    hourlyObs[, col] <- reps
  }
  names(hourlyObs) <- c("datetime", c("t", "p",  "u10", "ea", "qsi", "qli"))

  # convert precip by dividing by 3
  hourlyObs$p <- hourlyObs$p / 3
  
  # distribute Qsi to hourly
  qsi <- CRHMr::distributeQsi(obs3hr, QsiColnum = 5, latitude = latitude,
                              sunTimeOffset = 2, timeStep = 1, 
                              solarMethod = "PotSolarInst")
  hourlyObs$qsi <- qsi[,2]
  rm(obs3hr)
  
  # trim
  hourlyObs <- CRHMr::trimObs(hourlyObs)

  
  # write obs file
  obs1hrFile <- paste(outDir, locationName,"_",startYear,"-",endYear,
                      "_adjusted_1hr_CanRCM4.obs", sep = "")
  result <- CRHMr::writeObsFile(hourlyObs, obsfile = obs1hrFile, 
                      comment = "CanRCM4 adjusted downscaled to hourly")
  return(result)
}