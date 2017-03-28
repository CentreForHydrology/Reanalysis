#' Interpolates 3 and 6 hour WATCH variables to hourly values
#'
#' @description Interpolates \code{t}, \code{u10}, \code{ea}, \code{qsi}, and \code{qli} values (6-hour in \code{WFD}, 3-hour in \code{WFDEI}) to hourly. The total \code{p} values are evenly divided over the hourly intervals. Outputs hourly \code{t}, \code{u10}, \code{ea}, \code{qsi}, \code{qli}, and \code{p}.
#' @param infile Required. Name of file created by CreateWFDobs or CreateWFDobs.
#' @param outfile Required. Hourly obs file to be created.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not use
#' @return If successful, returns the value \code{TRUE} and writes the specified .obs file of hourly values. If unsuccessful, returns the value \code{FALSE}.
#' @author Kevin Shook
#' @examples
#'\dontrun{
#' WATCHhourlyObs('VermilionWATCH_WFD.obs', 'VermilionWFDhourly.obs')
#' WATCHhourlyObs('VermilionWATCH_WFDEI.obs', 'VermilionWFDEIhourly.obs')}
#' @export
WATCHhourlyObs <-
function(infile='', outfile='', logfile=''){
  # reads in 6 and 3 hour WATCH data and interpolates to hourly values
  # works with WFD and WFDEI data
  eol.val <- CRHMr::win.eol()
  # check parameters
  if (infile == ''){
    cat('Error: missing input .obs file\n')
    return(FALSE)
  }

  if (outfile == 0){
    cat('Error: missing output .obs file\n')
    return(FALSE)
  }
  classes <- c('character', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric' )
  original <- utils::read.table(infile, header=FALSE, skip=8, sep='\t',na.strings = "NA",
                         stringsAsFactors=FALSE, colClasses=classes)
  # read first line from .obs file to write to output
  fileName <- infile
  con <- file(fileName,open="r")
  original.header <- readLines(con, n=1)
  close(con)

  names(original) <- c('datetime', 't','e','u10','p', 'Qsi', 'Qli') #dp, added Qsi and Ql
  # create time series for interpolation
  start.datetime <- original$datetime[1]
  end.datetime <- original$datetime[length(original$datetime)]

  times <- format(timeDate::timeSequence(from = start.datetime, to = end.datetime,
                                     by='hour', format="%Y %m %d %H %M"),
                  format="%Y %m %d %H %M")
  time <- data.frame(times)
  names(time) <- 'datetime'
  filled <- merge(time, original, by='datetime', all=TRUE)

 # interpolate
  t <- zoo::na.approx(filled$t)
  e <- zoo::na.approx(filled$e)
  rh <- mapply(FUN='CRHMr::vp2rh', t, e)
  rh <- formatC(rh, digits=2, mode='real', format='f')
  t <- formatC(t, mode='real', digits=2, format='f')
  u10 <- format(zoo::na.approx(filled$u10), digits=2)
  p <- format(zoo::na.approx(filled$p, method='constant', f=0)/3, digits=2)
  Qsi <- zoo::na.approx(filled$Qsi)  #dp, added the line
  Qli <- zoo::na.approx(filled$Qli) #dp, added the line
  Qsi <- formatC(Qsi, mode='real', digits=2, format='f') #dp, added the line
  Qli <- formatC(Qli, mode='real', digits=2, format='f') #dp, added the liine

 # assemble into a dataframe

 # find end of data to account for trailing NA values
  end.loc <- min(length(t), length(rh), length(u10), length(p))
  datetime <- filled$datetime[1:end.loc]
  t <- t[1:end.loc]
  rh <- rh[1:end.loc]
  u10 <- u10[1:end.loc]
  p <- p[1:end.loc]
  Qsi <- Qsi[1:end.loc] #dp, added the line
  Qli <- Qli[1:end.loc] #dp, added the line
  output <- data.frame(datetime, t, rh, u10, p, Qsi, Qli)

 # find start of data to make sure it starts at 01:00
  first.times <- times[1:24]
  first.hours <- as.numeric(format(timeDate::timeDate(first.times, format = "%Y %m %d %H %M"),
                                   format='%H'))
  start.loc <- which(first.hours == 1)
 # format time and trim data

  output <- output[-(1:(start.loc-1)),]
  cat(original.header,  eol.val,
      't 1 (C)',  eol.val,
      'rh 1 (%)',  eol.val,
      'u10 1 (m/s)', eol.val,
      'p 1 (mm)',  eol.val,
      'Qsi 1 (W/m2)', eol.val, #dp, added the line
      'Qli 1 (W/m2)', eol.val, #dp, added the li
      '$ea ea(t, rh)',eol.val,
      '$u refwind(u10, 10, 2, 0.5)  "convert 10m wind speed to 2m wind speed"',  eol.val,
      '########################################',  eol.val, file=outfile, sep='')

  utils::write.table(output, file=outfile, sep='\t',col.names=FALSE, row.names=FALSE,
              quote=FALSE, eol = eol.val, append=TRUE)

  # record to logfile
  comment <- paste('WATCHhoulyrObs infile:', infile, ' outfile:',
                   outfile, sep='')
  result <- CRHMr::logAction(comment, logfile)

  return(result)
}
