#' Creates a CRHM .obs file of hourly values from WATCH reanalysis data WFDEI files
#'
#' @description Extracts data from WATCH WDFEI netCDF files and builds a CRHM .obs file of hourly values containing \code{t}, \code{ea}, \code{u10}, \code{p} \code{Qsi} and \code{Qli}. All values are interpolated from 3 and 6 hur data. The windspeeds are at 10m, hence they are denoted as u10. Air temperatures are at 2m. The values for ea are computed from the atmospheric pressure (at 10m) and the absolute humidity (at 2m). Because the original NetCDF files are very large, this function typically runs very slowly. Also, because this function assembles and processes all of the data in memory, it can require a lot of RAM to execute.
#' @param nc.location Required. A character string of the directory holding the WATCH WFDEI netCDf files. This is a file path WITHOUT a terminal slash, e.g. \option{z:\\WATCH\\WFDEI}
#' @param startyear Optional. Year to begin. Must be in the range \code{1979-2012}. Default is \code{1979}.
#' @param endyear Optional. Year to end. Must be in the range \code{1979-2012}. Default is \code{2012}.
#' @param lon Required. Decimal longitude to extract for.
#' @param lat Required. Decimal latitude to extract for.
#' @param precipType Optional. The precipitation type used; can be \option{CRU} (the default) or \option{GPCC}.
#' @param sunTimeOffset Optional. Number of hours that local noon is offset from solar noon. The default is 2 hours. This may be incorrect. It seems that WATCH data may not incorporate the temporal offset correctly - zero may work better.
#' @param solarMethod The method to be used for calculating the extra-terrestrial radiation. The default method is \option{simpleMaxSolar}. Note that this method is only valid for latitudes between 49 and 55\eqn{^\circ}{ }N. The other supported method is \option{PotSolarInst}, which requires the package \pkg{EcoHydRology} to be installed
#' @param interpolationMethod Optional. A vector containing the methods to be used for interpolation for each of the variables. Currently supported methods are \option{linear} and \option{spline}. The default is to use linear interpolation. If fewer methods than columns are specified, the methods are recycled.
#' @param obsFileName Required. Name of the .obs file to be created.
#' @param timezone Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. You can use \option{Etc/GMT+6} or \option{Etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If successful, returns the value \code{TRUE} and writes the specified .obs file. If unsuccessful, returns the value \code{FALSE}.
#' @author Kevin Shook.
#' @note Because this function can be slow, and uses a lot of memory, you may wish to run it repeatedly for short intervals. You can then use the function \code{appendObs} in \pkg{CRHMr} to join the files together. The time shifting results in obs files which do not begin and end on day boundaries, so you should use the function \code{trimObs} in \pkg{CRHMr} to trim the obs file - \emph{after} the final file has been assembled.
#' @references \R code for conversion of air pressure and absolute humidity was taken from project PEcAn The Predictive Ecosystem Analyzer \url{http://pecanproject.github.io}. The source code is available at \url{https://github.com/PecanProject/pecan/blob/master/modules/data.atmosphere/R/metutils.R.}
#' @seealso \code{\link{WATCHcreateHourlyWFDobs}} \code{\link[CRHMr]{trimObs}} \code{\link[CRHMr]{appendObs}}
#' @examples
#' \dontrun{
#' location <- 'z:\data\WATCH\WFDEI'
#' obsName <- 'VermilionWATCH_WFDEI.obs'
#' WATCHcreateHourlyWFDEIobs(nc.location=location, 
#' startyear=1979, endyear=2001, lon=-111.9, lat=53.2, timezone='Etc/GMT+7', obsFileName=obsName)
#' # read in file and trim it
#' obs <- CRHMr::readObsFile(obsFile=obsName, timezone='Etc/GMT+7')
#' trimmedObs <- CRHMr::trimObs(obs)
#' CRHMr::writeObsFile(trimmedObs, obsFileName)}
#' 
#' @export
WATCHcreateHourlyWFDEIobs <-
function(nc.location='', startyear=1979, endyear=2012,
                         lon=0, lat=0, precipType='CRU',sunTimeOffset=2, solarMethod='simpleMaxSolar',
         interpolationMethod='linear', obsFileName='',timezone='', quiet=TRUE, logfile=''){

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
  
  if (timezone == ''){
    cat('Error: must specify time zone\n')
    return(FALSE)
  }

  if (obsFileName == ''){
    cat('Error: missing name for .obs file\n')
    return(FALSE)
  }
  
  if ((precipType != 'CRU') & (precipType != 'GPCC')){
    cat('Error: unknown precipType\n')
    return(FALSE)
  }
  
  
  # calculate hour offset between timezone and GMT
  datetime <- '2000-01-01 01:00'
  s1 <- as.POSIXct(datetime, format='%Y-%m-%d %H:%M', tz='UTC')
  s2 <- as.POSIXct(datetime, format='%Y-%m-%d %H:%M', tz=timezone)
  houroffset <- -1*as.numeric(difftime(s2, s1, units='hours'))
  
  # get glat and glon
  lat.glat <- glat(lat)
  lon.glon <- glon(lon)

  # create header for .obs file
  obsHeader <- paste('Created from WATCH WFDEI reanalysis data. Lat = ',lat, ' Lon = ',lon,
      ' glat = ', lat.glat, 'glon = ', lon.glon, sep='')


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
          filename <- paste(nc.location, '/',variable,'_WFDEI_',precipType,'/',variable,
                            '_WFDEI_',precipType,'_',
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

      # having read all values for this month, do data conversions
      all.values$t <- all.values$Tair - 273.15         # K -> degrees C
      all.values$PSurf <- all.values$PSurf * 0.01    # Pa -> mb
      all.values$e <- all.values$Qair * all.values$PSurf / (0.378 * all.values$Qair + 0.622)
      # convert back to kPa and format
      all.values$ea <- all.values$e * 0.1
      all.values$p <- all.values$Rainf + all.values$Snowf * 10800  # mm/s -> mm

      all.values$u <- all.values$Wind
      all.values$Qsi <- all.values$SWdown            #dp, added the line
      all.values$Qli <- all.values$LWdown            #dp, added the line
      
      
      if((year == startyear) & (month == 1))
        obsData <- all.values[,c('DateTime', 't', 'ea', 'u',
                                 'p', 'Qsi', 'Qli')] #dp, added Qsi and Qli        
      else
        obsData <- rbind(obsData, all.values[,c('DateTime', 't', 'ea', 'u',
                                                'p', 'Qsi', 'Qli')])

    }
  }
  names(obsData)[1] <- 'datetime'
  # having read in all data, interpolate to hourly values
  # trim obs
  hourly_p <- CRHMr::distributeP(obs=obsData, p.cols=4, timestep=1)
  hourly_Qsi <- CRHMr::distributeQsi(obsData, QsiColnum = 5, latitude=lat, 
                             sunTimeOffset = sunTimeOffset, timeStep = 1, 
                             solarMethod = solarMethod)
  
  hourly_t_ea_u <- CRHMr::distributeInst(obsData, obsCols=c(1,2,3), timeStep=1, 
                                         interpolationMethod=interpolationMethod,
                                      maxLength=7, 
                                      quiet=quiet, logfile = logfile)
  
  hourly_Qli <- CRHMr::distributeQli(obsData, QliColnum = 6, tObs = hourly_t_ea_u, tColnum = 1,
                                     timeStep = 1, interpolationMethod = interpolationMethod)
  
  # merge all together
  all <- merge(hourly_t_ea_u, hourly_p, by='datetime')
  all <- merge(all, hourly_Qsi, by='datetime')
  all <- merge(all, hourly_Qli, by='datetime')
  names(all) <- c('datetime', 't', 'ea', 'u', 'p', 'Qsi', 'Qli')
  
  
  # record to logfile
  comment <- paste('WATCHcreateHourlyWFDEIobs startyear:',startyear,
                   ' endyear: ', endyear, ' lon: ',' lat: ',
                   lat, ' houroffset:',houroffset, ' obsFileName:', obsFileName, sep='')
  result <- CRHMr::logAction(comment, logfile)
  
  # write to file
  result <- CRHMr::writeObsFile(all, obsFileName, comment=obsHeader, quiet=quiet, logfile = logfile)
  return(result)
  
}
