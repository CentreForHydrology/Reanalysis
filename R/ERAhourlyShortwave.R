#' Estimates hourly incoming shortwave radiation from ERA-Interim 3-hourly values or ERA-40 6-hourly values
#' @name ERAhourlyShortwave
#' @description This function is called \emph{after} using the function \code{getNearestERAtimeseries}. This function removes negative values, deaccumulates the ERA-Interim 12-hour cumulative values, using the function \code{ERAdeaccum} and interpolates the 3-hour values to hourly values, based on the extra-terrestrial hourly radiation. ERA-40 6-hourly instantaneous values are interpolated directly to hourly values.
#' @param ERAssrd Required. The \pkg{CRHMr} obs dataframe of ERA ssrd values. The values must not be deaccumulated, as the \code{deaccumERA} function is called by this function.
#' @param ssrdColnum Optional. The column number containing the ssrd values, not including the datetime. Default is column 1.
#' @param latitude Required. The latitude of the point at which the ssrd data is extracted. This value is used to calculate the extra-terrestrial incoming solar radiation. Not used for ERA-40 data.
#' @param sunTimeOffset Optional. Number of hours that local noon is offset from solar noon. The default is 2 hours. Not used for ERA-40 data.
#' @param solarMethod Optional. The method to be used for calculating the extra-terrestrial radiation for 3-hourly ERA-Interim data. The default is \option{simpleMaxSolar}. The other supported method is \option{PotSolarInst}, which requires the package \pkg{EcoHydRology} to be installed.
#' @param interpolationMethod Optional. This is the method to be usied for interpolating instantaneous values from ERA-40. The default is \option{linear}, but \option{spline} can also be specified.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an \R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used
#' @return If successful, returns an obs dataframe containing the interpolated hourly ERA incoming shortwave radiation (Qsi). If unsuccessful, returns the value \code{FALSE}.
#' @author Kevin Shook
#' @export
#' @seealso \code{\link{ERAgetNearestTimeseries}}  \code{\link{ERAhourlyLongwave}} \code{\link[CRHMr]{distributeQsi}}
#'
#' @examples \dontrun{
#' hourlyQsi <- ERAhourlyShortwave(ssrd, latitude = 51.69)}
#'
ERAhourlyShortwave <- function(ERAssrd, ssrdColnum=1, latitude, 
                               sunTimeOffset=2, solarMethod='simpleMaxSolar', interpolationMethod = 'linear',
                               quiet=TRUE, logfile=''){

  # check parameters

  obsName <- deparse(substitute(ERAssrd))
  if (obsName == ''){
    cat('Error: must specify ssrd dataframe\n')
    return(FALSE)
  }
  if (nrow(ERAssrd) == 0){
    cat('Error: missing values\n')
    return(FALSE)
  }
  if (latitude == ''){
    cat('Error: must specify latitude\n')
    return(FALSE)
  }

  # remove negative values
  ERAssrd[,ssrdColnum+1] <- pmax(ERAssrd[,ssrdColnum+1], 0)
  
  # figure out interval from first 2 values
  first.datetime <- ERAssrd$datetime[1]
  second.datetime <- ERAssrd$datetime[2]
  interval.hours <- abs(as.numeric(difftime(second.datetime, first.datetime, units='hours')))
  
  # check time interval to decide whether data is -Interim or -40
  if (interval.hours == 3){
    # deaccumulate ERA 12-hour cumulative values to 3-hourly
    deaccum <- ERAdeaccum(ERAssrd, ssrdColnum, quiet=quiet, logfile=logfile)
    deaccum <- as.vector(deaccum[,ssrdColnum+1])
    ssrdColnum <- ssrdColnum + 1
  
    # replace negative values with original values
    deaccum[deaccum < 0] <- ERAssrd[deaccum < 0,ssrdColnum]
  
    # make sure the replacement values are > 0 as well
    deaccum < max(deaccum, 0)
  
    # now have 3-hour values in J/m2 - convert to W/m2
    deaccum <- deaccum / (3 * 3600)
  
    # deaccumulate 3-hour means to hourly values
    # create hourly dataframe
    datetime <- ERAssrd[,1]
    
    # create data frame 
    meanH <- data.frame(datetime, deaccum)
    
    # now interpolate
    hourlyQsi <-  CRHMr::distributeQsi(QsiObs=meanH, QsiColnum=1, 
                                       latitude=latitude, sunTimeOffset=sunTimeOffset, 
                                       timeStep=1, solarMethod=solarMethod, 
                                       quiet=TRUE, logfile='')
    
    names(hourlyQsi)[2] <- 'Qsi'
    
    
  } else{
    deaccum <- ERAdeaccum(ERAssrd, ssrdColnum, quiet=quiet, logfile=logfile)
    deaccum$ssrd <- deaccum$ssrd / (interval.hours * 3600)
    
    maxInterval <- interval.hours + 1
    # downscale longwave
    hourlyQsi <- CRHMr::distributeInst(deaccum, obsCols = 1, timeStep = 1, 
                                      interpolationMethod = interpolationMethod, maxLength = maxInterval,
                                      quiet = quiet, logfile = logfile)
    
    names(hourlyQsi)[2] <- 'Qsi'
  }


  # output info to screen (if req'd) and write to log file
  file.info <-  CRHMr::CRHM_summary(hourlyQsi)
  if (!quiet)
    print(file.info)

  comment <- paste('ERAhourlyShortwave ERAssrd:', obsName,
                   ' latitude:', latitude,
                   ' offset:', sunTimeOffset, sep='')
  result <-  CRHMr::logAction(comment, logfile)

  if(result)
    return(hourlyQsi)
  else
    return(result)

}
