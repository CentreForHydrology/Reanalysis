#' Creates a CRHM .obs file of 3-hourly values from WATCH reanalysis data WFDEI files
#'
#' @description Extracts data from WATCH WDFEI netCDF files and builds a CRHM .obs file of 3-hour data containing \code{t}, \code{ea}, \code{u10}, and \code{p} values. The output values can be interpolated to hourly using the function \code{HourlyWATCHobs}. The windspeeds are at 10m, hence they are denoted as u10. Air temperatures are at 2m. The values for ea are computed from the atmospheric pressure (at 10m) and the absolute humidity (at 2m).
#' @param nc.location Required. A character string of the directory holding the WATCH WFDEI netCDf files. This is a file path WITHOUT a terminal slash, e.g. \option{z:\\WATCH\\WFDEI}
#' @param startyear Optional. Year to begin. Must be in the range \code{1979-2012}. Default is \code{1979}.
#' @param endyear Optional. Year to end. Must be in the range \code{1979-2012}. Default is \code{2012}.
#' @param lon Required. Decimal longitude to extract for.
#' @param lat Required. Decimal latitude to extract for.
#' @param houroffset Required. Number of hours that the local location is offset from UTC (GMT). Must be negative in the western hemisphere. For Mountain Standard Time, the offset is \code{-7} hours.
#' @param obsFileName Required. Name of the .obs file to be created.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If successful, returns the value \code{TRUE} and writes the specified .obs file. Each month's data is written as it is created. If unsuccessful, returns the value \code{FALSE}.
#' @author Kevin Shook.
#' @references R code for conversion of air pressure and absolute humidity was taken from project PEcAn The Predictive Ecosystem Analyzer \url{http://pecanproject.github.io}. The source code is available at \url{https://github.com/PecanProject/pecan/blob/master/modules/data.atmosphere/R/metutils.R.}
#' @examples
#' \dontrun{
#' location <- 'z:\data\WATCH\WFDEI'
#' obsFileName <- 'VermilionWATCH_WFDEI.obs'
#' WATCHcreateWFDEIobs(location, 1979, 2001, -111.9, 53.2, -7, obsFileName)}
#' @export
WATCHcreateWFDEIobs <-
function(nc.location='', startyear=1979, endyear=2012,
                         lon=0, lat=0, houroffset=0, obsFileName='', quiet=TRUE, logfile=''){

  # check parameters
  if (nc.location == ''){
    cat('Error: missing location of netCDF files\n')
    return(FALSE)
  }

  if (lon == 0){
    cat('Error: missing longitude\n')
    return(FALSE)
  }

  if (lat == 0){
    cat('Error: missing latitude\n')
    return(FALSE)
  }

  if (obsFileName == ''){
    cat('Error: missing name for .obs file\n')
    return(FALSE)
  }

  # creates a CRHM .obs file from WATCH WFDEI reanalysis data
  eol.val <- CRHMr::win.eol()
  # get glat and glon
  lat.glat <- glat(lat)
  lon.glon <- glon(lon)

  # write header for .obs file
  cat('Created from WATCH WFDEI reanalysis data. Lat = ',lat, ' Lon = ',lon,
      ' glat = ', lat.glat, 'glon = ', lon.glon, eol.val,
      't 1 (C)',eol.val,
      'ea 1 (kPa)', eol.val,
      'u10 1 (m/s)', eol.val,
      'p 1 (mm)', eol.val,
      'Qsi 1 (W/m2)', eol.val,  #dp added this line
      'Qli 1 (W/m2)',  eol.val, #dp added this line
      '########################################',eol.val, file=obsFileName, sep='')

  vars <- c('Rainf','Snowf','Tair','Qair','Wind','PSurf', 'SWdown', 'LWdown')#dp, added SWdown and LWdown
  if (!quiet)
    cat('Reading data')
  # now get data files
  for (year in startyear:endyear){
    if (!quiet)
      cat('\nYear=', year,'\n', sep='')
    for (month in 1:12){
      if (!quiet)
        cat(month,' ',sep='')
      # get each variable for each month
      month.str <- formatC(month, width = 2, format = "d", flag = "0")
      for (variable in vars){
        # create file name & open .nc file
        if (variable == 'Rainf' | variable == 'Snowf')
          filename <- paste(nc.location, '/',variable,'_WFDEI_CRU/',variable,'_WFDEI_CRU_',
                            year,month.str,'.nc', sep='')
        else
          filename <- paste(nc.location, '/',variable,'_WFDEI/',variable,'_WFDEI_',
                            year,month.str,'.nc', sep='')
        values <- WATCHgetWFDEIvalues(lon.glon, lat.glat, filename, variable, month, year, houroffset)

        if (variable == vars[1])
          all.values <- values
        else
          all.values <- merge(all.values, values, by='DateTime', all=TRUE)
      }

      # having read all values for this month, do data conversions and write .obs file
      # calculate ea using code from https://github.com/PecanProject/pecan/blob/master/modules/data.atmosphere/R/metutils.R
      all.values$t <- format(all.values$Tair - 273.15, digits=1)    # K -> C
      all.values$PSurf <- all.values$PSurf * 0.01    # Pa -> mb
      all.values$e <- all.values$Qair * all.values$PSurf / (0.378 * all.values$Qair + 0.622)
      # convert back to kPa and format
      all.values$ea <- format(all.values$e * 0.1, digits=3)
      all.values$p <- format((all.values$Rainf + all.values$Snowf) * 10800,
                             digits=2) # mm/s -> mm

      all.values$u <- format(all.values$Wind, digits=2)
      all.values$Qsi <- format(all.values$SWdown, digits=2) #dp, added the line
      all.values$Qli <- format(all.values$LWdown, digits=2) #dp, added the line

      ObsData <- all.values[,c('DateTime', 't', 'ea', 'u',
                               'p', 'Qsi', 'Qli')] #dp, added Qsi and Qli

      ObsData$DateTime <- format(ObsData$DateTime, format='%Y %m %d %H %M')
      utils::write.table(ObsData, file=obsFileName, sep='\t',col.names=FALSE, row.names=FALSE,
                  quote=FALSE, eol = eol.val, append=TRUE)
    }
  }
  # record to logfile
  comment <- paste('createWFDEIobs startyear:',startyear,
                   ' endyear: ', endyear, ' lon: ',' lat: ',
                   lat, ' houroffset:',houroffset, ' obsFileName:', obsFileName, sep='')
  result <- CRHMr::logAction(comment, logfile)
  return(result)
}
