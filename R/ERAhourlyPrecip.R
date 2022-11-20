#' Etimates hourly precipitation from ERA cumulatve 3 hourly values
#' @name ERAhourlyPrecip
#' @description This function is called \emph{after} using the function \code{getNearestERAtimeseries}. This function removes negative values, deaccumulates the 12-hour cumulative values to 3-hour values, using the function \code{deaccumERA}, and divides the results by 3 to give to hourly values.
#' @param ERAtp  Required. The \pkg{CRHMr} obs dataframe of ERA tp values. The values must not be deaccumulated, as the \code{deaccumERA} function is called by this function.
#' @param tpColnum Optional. The column number containing the tp values, not including the datetime. Default is column 1.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used
#' @author Kevin Shook
#' @return If successful, returns an obs dataframe containing the interpolated hourly precipitation (p). If unsuccessful, returns the value \code{FALSE}.
#' @seealso \code{\link{ERAgetNearestTimeseries}}
#' @export
#'
#' @examples \dontrun{
#' hourlyP <- ERAhourlyPrecip(tp)}
ERAhourlyPrecip <- function(ERAtp, tpColnum=1, quiet=TRUE, logfile=''){

  obsName <- deparse(substitute(ERAtp))
  if (obsName == ''){
    cat('Error: must specify ERA tp dataframe\n')
    return(FALSE)
  }
  if (nrow(ERAtp) == 0){
    cat('Error: missing values\n')
    return(FALSE)
  }

  # check to see if values have been converted from m to mm
  if (max(ERAtp[,tpColnum+1]) < 1){
    if(!quiet)
      cat('Values appear to be in m, will be converted to mm\n')
    ERAtp[,tpColnum+1] <- ERAtp[,tpColnum+1] * 1000
  }



  # deaccumulate precipitation from 12h to 3h
  # remove negative values
  ERAtp[,tpColnum+1] <- pmax(ERAtp[,tpColnum+1], 0)

  # deaccumulate ERA 12-hour cumulative values to 3-hourly
  deaccum <- ERAdeaccum(ERAtp, tpColnum, quiet=quiet, logfile=logfile)
  deaccum <- as.vector(deaccum[,tpColnum+1])
  tpColnum <- tpColnum + 1

  # replace negative values with original values
  deaccum[deaccum < 0] <- ERAtp[deaccum < 0,tpColnum]

  # make sure the replacement values are > 0 as well
  deaccum < max(deaccum, 0)

  # deaccumulate 3-hour means to hourly values
  # create hourly dataframe
  datetime <- ERAtp$datetime

  # figure out interval from first 2 values
  first.datetime <- datetime[1]
  second.datetime <- datetime[2]

  last.datetime <- datetime[length(datetime)]

  # set back by 2 intervals
  interval.hours <- abs(as.numeric(difftime(second.datetime, first.datetime, units='hours')))
  first.datetime <- first.datetime - (interval.hours - 1) * 3600
  datetime <- seq(from=first.datetime, to=last.datetime, by=3600)



  # create 3 x 1-hourly values for each 3-hour mean
  deaccum <- deaccum / 3
  p <- rep(deaccum, each=interval.hours)

  hourly <- data.frame(datetime, p)

  # output info to screen (if req'd) and write to log file
  file.info <-  CRHMr::CRHM_summary(hourly)
  if (!quiet)
    print(file.info)

  comment <- paste('ERAhourlyPrecip ERAtp:', obsName, sep='')
  result <-  CRHMr::logAction(comment, logfile)

  if(result)
    return(hourly)
  else
    return(result)

}
