#' Deaccumulates ERA5 cumulative precipitation time series
#' @description This function is used to deaccumulate variables stored as 24-hour cumulative values by ERA to 1-hour values. Note that this function only works for a single location, i.e. \emph{NOT} areal values.
#' @param ERAobs Required. A \pkg{CRHMr} obs data frame of ERA 5 data, created by \code{ERAgetnearestTimeSeries}.
#' @param colnum Optional. The column number containing the values to be deaccumulated, not including the datetime. Default is column 1.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#'@param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#'
#' @return If successful, returns an obs data frame containing the deaccumulated value. If unsuccessful, returns the value \code{FALSE}.
#' @author Kevin Shook
#' @export
#'
#' @examples \dontrun{
#' deaccum <- ERA5deaccum(ERAobs)}
#'
ERAdeaccum <- function(ERAobs, colnum=1, quiet=TRUE, logfile=''){
  # check parameters

  obsName <- deparse(substitute(ERAobs))
  if (obsName == ''){
    cat('Error: must specify ERA dataframe\n')
    return(FALSE)
  }
  if (nrow(ERAobs) == 0){
    cat('Error: missing values\n')
    return(FALSE)
  }

  # do deaccumulation
  colnum <- colnum + 1
  deaccum <-  c(ERAobs[1,colnum], diff(ERAobs[,colnum]))

  hours <- as.numeric(format(ERAobs$datetime, "%H"))
  ERAobs[hours != 1, colnum] <- deaccum[hours != 1]


  # output info to screen (if req'd) and write to log file
  file.info <- CRHMr::CRHM_summary(ERAobs)
  if (!quiet)
    print(file.info)

  comment <- paste('ERAdeaccum ERAobs:', obsName, sep='')
  result <- CRHMr::logAction(comment, logfile)

  if(result)
    return(ERAobs)
  else
    return(result)

}
