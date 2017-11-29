#' Creates a CRHM .obs file of 3-hourly and 6-hourly values from WATCH reanalysis data WFD files
#'
#' @description Extracts data from WATCH WFD netCDF files and builds a CRHM .obs file of 3-hour data containing \code{t}, \code{ea}, \code{u10}, and \code{p} values. The values of \code{t}, \code{ea} and \code{u} are 6-hourly, with \code{NA} values inserted. The data are output at MST. The output values can be interpolated to hourly values using the function \code{HourlyWATCHObs}. The windspeeds are at 10m, so they are denoted as \code{u10}. Air temperatures are at 2m. The values for \code{ea} are computed from the atmospheric pressure (at 10m) and the absolute humidity (at 2m).
#' @param nc.location Required. A character string of the directory holding the WATCH WFD netCDf files. This is a file path WITHOUT a terminal slash, e.g. \option{z:\\WATCH\\WFD}.
#' @param startyear Optional. Year to begin. Must be in the range \code{1901-2001}. Default is \code{1901}.
#' @param endyear Optional. Year to end. Must be in the range \code{1901-2001}. Default is \code{2001}.
#' @param lon Required. Decimal longitude to extract for.
#' @param lat Required. Decimal latitude to extract for.
#' @param houroffset Required. Number of hours that the local location is offset from UTC (GMT). Must be negative in the western hemisphere. For Mountain Standard Time, the offset is \code{-7} hours.
#' @param obsFileName Required. Name of the .obs file to be created.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If successful, returns the value \code{TRUE} and writes the specified .obs file. Each month's data is written as it is created. If unsuccessful, returns the value \code{FALSE}.
#' @author Kevin Shook
#' @seealso \code{\link{WATCHcreateWFDEIobs}}
#' @references R code for conversion of air pressure and absolute humidity was taken from project PEcAn The Predictive Ecosystem Analyzer \url{http://pecanproject.github.io}. The source code is available at \url{https://github.com/PecanProject/pecan/blob/master/modules/data.atmosphere/R/metutils.R.}
#' @examples
#' \dontrun{
#' location <- 'z:\data\WATCH\WFD'
#' obsFileName <- 'VermilionWATCH_WFD.obs'
#' WATCHcreateWFDobs(location, 1979, 2001, -111.9, 53.2, -7, obsFileName)}
#' @export
WATCHcreateWFDobs <-
function(nc.location='', startyear=1901, endyear=2001,
         lon=0, lat=0, houroffset=0, obsFileName='', quiet=TRUE, logfile=''){

  # setup variables
  data("Reanalysis", envir = environment())
  land <- land
  eol.val <- CRHMr::win.eol()

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

  # get closest half degree to specified lat and lon
  lat.halfdegree <- halfdegree(lat)
  lon.halfdegree <- halfdegree(lon)

  # write header for .obs file
  cat('Created from WATCH WFD reanalysis data. Lat = ',lat.halfdegree, ' Lon = ',lon.halfdegree,
       eol.val,'t 1 (C)',eol.val,'ea 1 (kPa)',
       eol.val,'u10 1 (m/s)',
       eol.val,'p 1 (mm)',
       eol.val,'Qsi 1 (W/m2)', #dp, added the line
       eol.val,'Qli 1 (W/m2)', #dp, added the line
       eol.val,
      '########################################',eol.val, file=obsFileName, sep='')

  # find lat and lon in landfile
  land.selected <- land[(land$Longitude == lon.halfdegree) &
                          (land$Latitude == lat.halfdegree),]
  loc <- which(land$Land == land.selected$Land)
  vars <- c('Rainf','Snowf','Tair','Qair','Wind','PSurf', 'SWdown', 'LWdown')
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
          filename <- paste(nc.location, '/',variable,'_WFD/',variable,'_WFD_CRU_',
                            year,month.str,'.nc', sep='')
        else
          filename <- paste(nc.location, '/',variable,'_WFD/',variable,'_WFD_',
                            year,month.str,'.nc', sep='')
        values <- WATCHgetWFDvalues(loc, filename, variable, month, year, houroffset)

        if (variable == vars[1])
          all.values <- values
        else
          all.values <- merge(all.values, values, by='DateTime', all=TRUE)
      }

      # having read all values for this month, do data conversions and write .obs file
      # calculate ea using code from https://github.com/PecanProject/pecan/blob/master/modules/data.atmosphere/R/metutils.R
      all.values$t <- format(all.values$Tair - 273.15, digits=1)    # K -> °C
      all.values$PSurf <- all.values$PSurf * 0.01    # Pa -> mb
      all.values$e <- all.values$Qair * all.values$PSurf / (0.378 * all.values$Qair + 0.622)
      # convert back to kPa and format
      all.values$ea <- format(all.values$e * 0.1, digits=3)
      all.values$p <- format((all.values$Rainf +
                                all.values$Snowf) * 10800, digits=2) # mm/s -> mm

      all.values$u <- format(all.values$Wind, digits=2)
      all.values$Qsi <- format(all.values$SWdown, digits=2) #dp, added the line
      all.values$Qli <- format(all.values$LWdown, digits=2) #dp, added the line

      ObsData <- all.values[,c('DateTime', 't', 'ea', 'u',
                               'p', 'Qsi', 'Qli')] #dp, added Qsi and Qli

      ObsData$DateTime <- format(ObsData$DateTime, format='%Y %m %d %H %M')
      utils::write.table(ObsData, file=obsFileName, sep='\t',col.names=FALSE, row.names=FALSE,
                  quote=FALSE, eol =  eol.val, append=TRUE)
    }
  }
  # record to logfile
  comment <- paste('createWFDobs startyear:',startyear,
                   ' endyear: ', endyear, ' lon: ',' lat: ',
                   lat, ' houroffset:',houroffset,
                   ' obsFileName:', obsFileName, sep='')
  result <- logAction(comment, logfile)

  return(result)
}
