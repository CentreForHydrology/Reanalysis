#' Estimates hourly air temperatures from ERA 3 hourly values
#' @name ERAhourlyAirtemp
#' @description Interpolates ERA 3-hour instantaneous air temperatures to hourly values. The ERA air temperatures are first converted from K to \eqn{^\circ}{}C, if required.
#' @param ERAt2m Required. The \pkg{CRHMr} obs dataframe of ERA t2m values. The values must not be deaccumulated, as the \code{deaccumERA} function is called by this function.
#' @param t2mColnum Optional. The column number containing the t2m values, not including the datetime. Default is column 1.
#' @param method Optional. The methods to be used for interpolation of the air temperature. Currently supported methods are \option{linear} and \option{spline}. The default is to use linear interpolation.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @author Kevin Shook
#' @return If successful, returns an obs dataframe containing the interpolated hourly 2m air temperatures. If unsuccessful, returns the value \code{FALSE}.
#' @seealso \code{\link{ERAgetNearestTimeseries}} \code{\link[CRHMr]{interpolate}}
#' @export
#'
#' @examples \dontrun{
#' hourlyTemps <- ERAhourlyAirtemp(t2m)}
ERAhourlyAirtemp <- function(ERAt2m, t2mColnum=1, method='linear', quiet=TRUE, logfile=''){
  # check parameters

  obsName <- deparse(substitute(ERAt2m))
  if (obsName == ''){
    cat('Error: must specify ERA t2m dataframe\n')
    return(FALSE)
  }
  if (nrow(ERAt2m) == 0){
    cat('Error: missing values\n')
    return(FALSE)
  }

  # check units - convert from K to C if required
  ERAt2m <- ERAt2m[, c(1, 1+t2mColnum)]
  minTemp <- max(ERAt2m[,2])

  # check to see if values have already been converted from K to C
  if(minTemp > 150){
    if(!quiet)
      cat('Values appear to be in K, will be converted to C')
    ERAt2m[,2] <- ERAt2m[,2] - 273.15
  }

  # interpolate to hourly
  # first generate hourly time series
  firstDatetime <- ERAt2m$datetime[1]
  lastDatetime <- ERAt2m$datetime[nrow(ERAt2m)]

  # generate hourly datetimes
  datetime <- seq(from=firstDatetime, to=lastDatetime, by=3600)

  # merge
  hourly <- data.frame(datetime)
  merged <- merge(hourly, ERAt2m, by='datetime', all.x=TRUE)

  # now interpolate
  hourlyT <- CRHMr::interpolate(merged, varcols=1, methods=method, maxlength=4,
                              quiet, logfile)

  names(hourlyT) <- c('datetime', 't')

  # output log files
  obs.info <-  CRHMr::CRHM_summary(hourlyT)
  if (!quiet)
    print(obs.info)

  comment <- paste('ERAhourlyAirtemp:', obsName,
                   ' method:', method,
                   sep='')

  result <-  CRHMr::logAction(comment, logfile)


  if (result)
    return (hourlyT)
  else
    return(result)

}
