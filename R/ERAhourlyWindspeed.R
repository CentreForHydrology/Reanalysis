#' Estimates hourly wind speeds from ERA 3 hourly wind vectors
#' @name ERAhourlyWindspeed
#' @param ERAu10 Required. The \pkg{CRHMr} obs dataframe of ERA u10 wind vector values. The values must not be deaccumulated, as the \code{deaccumERA} function is called by this function.
#' @param u10Colnum Optional. The column number containing the u10 values, not including the datetime. Default is column 1.
#' @param ERAv10 Required. The \pkg{CRHMr} obs dataframe of ERA v10 wind vector values. The values must not be deaccumulated, as the \code{deaccumERA} function is vcalled by this function.
#' @param v10Colnum Optional. The column number containing the v10 values, not including the datetime. Default is column 1.
#' @param method Optional. The methods to be used for interpolation of the wind speeds. Currently supported methods are \option{linear} and \option{spline}. The default is to use linear interpolation.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @author Kevin Shook
#' @return If successful, returns an obs dataframe containing the interpolated hourly 10m wind speeds (u10) in m/s and directions in degrees. If unsuccessful, returns the value \code{FALSE}.
#' @seealso \code{\link{ERAgetNearestTimeseries}} \code{\link[CRHMr]{interpolate}}
#' @export
#'
#' @examples \dontrun{
#' hourlyU <- ERAhourlyWindspeed(ERAu10=u10, ERAv10=v10)}
ERAhourlyWindspeed <- function(ERAu10, u10Colnum=1, ERAv10, v10Colnum=1,
                            method='linear', quiet=TRUE, logfile=''){

  # check parameters
  uName <- deparse(substitute(ERAu10))
  if (uName == ''){
    cat('Error: must specify ERA u dataframe\n')
    return(FALSE)
  }
  if (nrow(ERAu10) == 0){
    cat('Error: missing values for ERAu \n')
    return(FALSE)
  }

  vName <- deparse(substitute(ERAv10))
  if (vName == ''){
    cat('Error: must specify ERA v dataframe\n')
    return(FALSE)
  }
  if (nrow(ERAv10) == 0){
    cat('Error: missing values for ERAv \n')
    return(FALSE)
  }

  # subset dataframes
  ERAu10 <- ERAu10[,c(1,u10Colnum+1)]
  ERAv10 <- ERAv10[,c(1,v10Colnum+1)]


  # merge dataframes together

  merged <- merge(ERAu10, ERAv10, by='datetime')
  names(merged)[2:3] <- c('u10', 'v10')

  # interpolate to hourly
  # first generate hourly time series
  firstDatetime <- merged$datetime[1]
  lastDatetime <- merged$datetime[nrow(merged)]

  # generate hourly datetimes
  datetime <- seq(from=firstDatetime, to=lastDatetime, by=3600)

  # merge
  hourly <- data.frame(datetime)
  merged2 <- merge(hourly, merged, by='datetime', all.x=TRUE)

  # now interpolate
  hourly <-  CRHMr::interpolate(merged2, varcols=c(1,2), methods=method, maxlength=4,
                             quiet, logfile)

  hourly$windspeed <- (hourly$u10 ^ 2 + hourly$v10 ^ 2) ^ 0.5
  angle <- atan2(-hourly$u10, -hourly$v10)

  # convert to degrees and orient to north
  hourly$direction <- (angle * 180 / pi) + 180

  hourly <- hourly[, c('datetime', 'windspeed', 'direction')]
  names(hourly)[2] <-'u10'

  # output log files
  obs.info <-  CRHMr::CRHM_summary(hourly)
  if (!quiet)
    print(obs.info)

  comment <- paste('ERAhourlyWindspeed u_values:', uName,
                   'v_values:', vName,
                   ' method:', method,
                   sep='')

  result <-  CRHMr::logAction(comment, logfile)


  if (result)
    return (hourly)
  else
    return(result)

}
