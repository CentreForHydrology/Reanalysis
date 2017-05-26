#' Calculates the hourly vapour pressure from the 3-hour ERA-Interim or 6 hour ERA-40 dew point temperature
#' @name ERAhourlyVP
#' @description Interpolates ERA instantaneous dew point temperatures to hourly values. The ERA dew point temperatures are first converted from K to \eqn{^\circ}{}C, if required.
#' @param ERAd2m Required. The \pkg{CRHMr} obs dataframe of ERA d2m values. The values must not be deaccumulated, as the \code{deaccumERA} function is called by this function.
#' @param d2mColnum Optional. The column number containing the d2m values, not including the datetime. Default is column 1.
#' @param method Optional. The methods to be used for interpolation of the dew point temperature. Currently supported methods are \option{linear} and \option{spline}. The default is to use linear interpolation.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @author Kevin Shook
#' @return If successful, returns an obs dataframe containing the interpolated hourly 2m vapour pressures (ea) in kPA. If unsuccessful, returns the value \code{FALSE}.
#' @seealso \code{\link{ERAgetNearestTimeseries}} \code{\link[CRHMr]{interpolate}} \code{\link{ERAhourlyAirtemp}}
#' @export
#'
#' @examples \dontrun{
#' hourlyVP <- ERAhourlyVP(d2m)}
ERAhourlyVP <- function(ERAd2m, d2mColnum=1, method='linear', quiet=TRUE, logfile=''){
  obsName <- deparse(substitute(ERAd2m))
  if (obsName == ''){
    cat('Error: must specify ERA d2m dataframe\n')
    return(FALSE)
  }
  if (nrow(ERAd2m) == 0){
    cat('Error: missing values\n')
    return(FALSE)
  }

  # check units - convert from K to C if required
  ERAd2m <- ERAd2m[, c(1, 1+d2mColnum)]
  minTemp <- max(ERAd2m[,2])

  # check to see if values have already been converted from K to C
  if(minTemp > 150){
    if(!quiet)
      cat('Values appear to be in K, will be converted to C')
    ERAd2m[,2] <- ERAd2m[,2] - 273.15
  }
  
  # get timestep
  dt  <- difftime(ERAd2m$datetime[1], ERAd2m$datetime[2], units='hours')
  dt <- abs(as.numeric(dt))
  
  # get interpolation max length
  interpolationLength <- dt + 1

  # interpolate to hourly
  # first generate hourly time series
  firstDatetime <- ERAd2m$datetime[1]
  lastDatetime <- ERAd2m$datetime[nrow(ERAd2m)]

  # generate hourly datetimes
  datetime <- seq(from=firstDatetime, to=lastDatetime, by=3600)

  # merge
  hourlyd2m <- data.frame(datetime)
  merged <- merge(hourlyd2m, ERAd2m, by='datetime', all.x=TRUE)

  # now interpolate
  hourlydt2m <-  CRHMr::interpolate(merged, varcols=1, methods=method, maxlength=interpolationLength,
                             quiet, logfile)

  # get vapour pressure
  hourlyVP <-  CRHMr::saturatedVP(hourlydt2m[,2])
  hourlyVP <- data.frame(datetime, hourlyVP)
  names(hourlyVP)[2] <- 'ea'

  # output log files
  obs.info <-  CRHMr::CRHM_summary(hourlyVP)
  if (!quiet)
    print(obs.info)

  comment <- paste('ERAhourlyVP:', obsName,
                   ' method:', method,
                   sep='')

  result <-  CRHMr::logAction(comment, logfile)

  if (result)
    return (hourlyVP)
  else
    return(result)

}
