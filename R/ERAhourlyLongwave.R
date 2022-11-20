#' Estimates hourly incoming longwave radiation from ERA-Interim 3-hourly values or ERA-40 6-hourly values
#' @name ERAhourlyLongwave
#' @description This function is called \emph{after} using the function \code{getNearestERAtimeseries}. This function deaccumulates the ERA-Interim 12-hour cumulative values,using the function \code{ERAdeaccum} and interpolates the 3-hour values to hourly values, based on the extra-terrestrial hourly radiation. ERA-40 6-hourly instantaneous values are interpolated directly to hourly values.
#' @param ERAstrd Required. The \pkg{CRHMr} obs dataframe of ERA strd values. The values must not be deaccumulated, as the \code{ERAdeaccum} function is called by this function.
#' @param strdColnum Optional. The column number containing the strd values, not including the datetime. Default is column 1.
#' @param ERAt2m Required. The \pkg{CRHMr} obs dataframe of ERA t2m values. The values must not be deaccumulated, as the \code{deaccumERA} function is called by this function.
#' @param t2mColnum Optional. The column number containing the t2m values, not including the datetime. Default is column 1.
#' @param method Optional. The methods to be used for interpolation of the air temperature. Currently supported methods are \option{linear} and \option{spline}. The default is to use linear interpolation.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used
#' @return If successful, returns an obs dataframe containing the interpolated hourly ERA incoming longwave radiation (Qli) in W/m\eqn{^2}{^2}, and the hourly 2m air temperatures (t) in \eqn{^\circ}{}C. If unsuccessful, returns the value \code{FALSE}.
#' @author Kevin Shook
#' @export
#' @seealso \code{\link{ERAgetNearestTimeseries}}  \code{\link{ERAdeaccum}} \code{\link{ERAhourlyShortwave}}
#'
#' @references The value of the Stefan-Boltzmann constant was obtained from \url{http://physics.nist.gov/cgi-bin/cuu/Value?sigma}
#' @examples \dontrun{
#' hourlyQli <- ERAhourlyLongwave(ERAstrd=strd, ERAt2m=t2m)}
ERAhourlyLongwave <- function(ERAstrd, strdColnum=1, ERAt2m, t2mColnum=1,
                              method='linear', quiet=TRUE, logfile=''){

  # check parameters
  obsName1 <- deparse(substitute(ERAstrd))
  if (obsName1 == ''){
    cat('Error: must specify ERA strd dataframe\n')
    return(FALSE)
  }
  if (nrow(ERAstrd) == 0){
    cat('Error: missing values\n')
    return(FALSE)
  }

  obsName2 <- deparse(substitute(ERAt2m))
  if (obsName2 == ''){
    cat('Error: must specify ERA t2m dataframe\n')
    return(FALSE)
  }
  if (nrow(ERAt2m) == 0){
    cat('Error: missing values\n')
    return(FALSE)
  }



  airtemp <- ERAt2m[, c(1,t2mColnum+1)]
  
  # figure out interval from first 2 values
  first.datetime <- ERAstrd$datetime[1]
  second.datetime <- ERAstrd$datetime[2]
  interval.hours <- abs(as.numeric(difftime(second.datetime, first.datetime, units='hours')))
  
  # check time interval to decide whether data is -Interim or -40
  
  # deaccumulate longwave
  if (interval.hours == 3){
    longwaveDeaccum <- ERAdeaccum(ERAstrd, colnum=strdColnum, quiet=quiet, 
                                  logfile=logfile)
    # convert longwave from J/m2 to W/m2
    longwaveDeaccum$strd <- longwaveDeaccum$strd / (interval.hours * 3600)
    
    # downscale longwave
    hourlyLW <- CRHMr::distributeQli(longwaveDeaccum, 1, airtemp, 1)
    
  } else {
    longwaveDeaccum <- ERAstrd
    # convert longwave from J/m2 to W/m2
    longwaveDeaccum$strd <- longwaveDeaccum$strd / (interval.hours * 3600)
    
    maxInterval <- interval.hours + 1
    # downscale longwave
    hourlyLW <- CRHMr::distributeInst(longwaveDeaccum, obsCols = 1, timeStep = 1, 
                                      interpolationMethod = method, maxLength = maxInterval,
                                      quiet = quiet, logfile = logfile)
  }
    
  # output log files
  obs.info <-  CRHMr::CRHM_summary(hourlyLW)
  if (!quiet)
    print(obs.info)

  comment <- paste('ERAhourlyLongwave:', obsName1, sep='')
  result <-  CRHMr::logAction(comment, logfile)


  if (result)
    return (hourlyLW)
  else
    return(result)

}
