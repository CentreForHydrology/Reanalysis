#' Create obs files from WATCH WFD data for groups of sites
#'
#' @description Extracts all variables from WFD NetCDF files and builds obs files for all specified locations
#' @param siteFile Required. A .csv file containing all of the variables required to describe the site locations. These are: 
#' \describe{
#'  \item{Name}{Name of site}
#'  \item{Land}{Land ID number}
#'  \item{Longitude}{Site longitude}
#'  \item{Latitude}{Site latitude}
#'  \item{glon}{Site glon value}
#'  \item{glat}{Site glat value}
#'  \item{timezone}{Time zone for site. Must be in Etc format, e.g. Etc/GMT+7}
#'  \item{SolarOffset}{Local offset of solar noon in hours} 
#' }
#' @param ncLocation Required. Location of all of the WATCH WFD files. Must have a trailing '\\' symbol.
#' @param outputLocation  Required. Location for all of the output files. Must have a trailing '\\' symbol.
#' @param startyear Optional. Year to begin extraction. Default is 1901.
#' @param endyear Optional. Year to begin extraction. Default is 2001.
#' @param solarMethod The method to be used for calculating the extra-terrestrial radiation. The default method is \option{PotSolarInst}, which requires the package \pkg{EcoHydRology} to be installed. The other supported method is  \option{simpleMaxSolar}. Note that this method is only valid for latitudes between 49 and 55\eqn{^\circ}{ }N.
#' @param interpolationMethod Optional. A vector containing the methods to be used for interpolation for each of the variables. Currently supported methods are \option{linear} and \option{spline}. The default is to use linear interpolation. If fewer methods than columns are specified, the methods are recycled.
#' @param quiet Optional. Suppresses display of messages, except for errors. If you are calling this function in an \R script, you will usually leave \code{quiet=TRUE} (i.e. the default). If you are working interactively, you will probably want to set \code{quiet=FALSE}.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used.
#' @return If successful, returns the value \code{TRUE} and writes the specified .obs file. Each month's data is written as it is created. If unsuccessful, returns the value \code{FALSE}.
#' @author Kevin Shook
#' @references R code for conversion of air pressure and absolute humidity was taken from project PEcAn The Predictive Ecosystem Analyzer \url{http://pecanproject.github.io}. The source code is available at \url{https://github.com/PecanProject/pecan/blob/master/modules/data.atmosphere/R/metutils.R.}
#' @export
#'
#' @examples \dontrun{
#' points <- './Reanalysis/testdata/WATCH_selected_points.csv'
#' outputLocation <- './Reanalysis/testdata/'
#' ncloc <- '//water.usask.ca/Centre/Reanalysis_data/WATCH/WFD'
#' WATCHgroupWFDobs(siteFile = points, ncLocation = ncloc, outputLocation=outputLocation, 
#' startyear = 1962, endyear = 1962))
#' }
WATCHgroupWFDobs <- function(siteFile, ncLocation, outputLocation, startyear=1901, 
                             endyear=2001, solarMethod='PotSolarInst', 
                             interpolationMethod='linear', quiet=TRUE, logfile=''){
  eol.val <- CRHMr::win.eol()
  land <- land
  # check parameters
  if (is.null(ncLocation ) | (ncLocation == '')){
    cat('Error: missing location of netCDF files\n')
    return(FALSE)
  }
  
  if (is.null(siteFile ) | (siteFile == '')){
    cat('Error: missing site file\n')
    return(FALSE)
  }
  
  sites <- read.csv(siteFile, header=TRUE, stringsAsFactors = FALSE)
  site_count <- nrow(sites)
  
  file_vars <- c('Rainf_WFD_CRU','Snowf_WFD_CRU','Rainf_WFD_GPCC','Snowf_WFD_GPCC','Tair_WFD','Qair_WFD',
            'Wind_WFD','PSurf_WFD', 'SWdown_WFD', 'LWdown_WFD')
  vars <- c('Rainf','Snowf','Rainf','Snowf','Tair','Qair','Wind','PSurf', 'SWdown', 'LWdown')
  num_vars <- length(vars)
  
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
      for (var_num in 1:num_vars){
        variable <- vars[var_num]
        file_var <- file_vars[var_num]
        
        # create file name & open .nc file
        nc_file <- paste(ncLocation, variable,'_WFD/',file_var,'_', year , month.str, '.nc', sep='')
        nc <- RNetCDF::open.nc(nc_file, write=FALSE)
        nctimes <- RNetCDF::var.get.nc(nc, variable='time', unpack = TRUE)
        all.values <- RNetCDF::var.get.nc(nc, variable=variable, unpack = TRUE)
        RNetCDF::close.nc(nc)
        
        count.time <- length(nctimes)
        # depending on number of values, figure out time interval
        if (count.time > 125)
          interval <- '3 hours'
        else
          interval <- '6 hours'
        
        
        # create origin, beginning at midnight on first day of the month
        origin.string <- paste(year,'-',month,'-01 00:00:00',sep='')
        origin.datetime <- timeDate::timeDate(origin.string, format = "%Y-%m-%d %H:%M:%S")
        
        # create time sequence
        times.gmt <- timeDate::timeSequence(from = origin.datetime,
                                            by=interval, length.out=count.time,
                                            format="%Y %m %d %H %M %S", zone='UTC')
        
        for (site_num in 1:site_count){
          # get values for each site
          site_loc <- which(land$Land == sites$Land[site_num])
          lon <- land$Longitude[site_loc]
          lat <- land$Latitude[site_loc]
          timezone <- sites$timezone[site_num]
          
          # calculate hour offset between timezone and GMT
          datetime <- '2000-01-01 01:00'
          s1 <- as.POSIXct(datetime, format='%Y-%m-%d %H:%M', tz='UTC')
          s2 <- as.POSIXct(datetime, format='%Y-%m-%d %H:%M', tz=timezone)
          houroffset <- -1*as.numeric(difftime(s2, s1, units='hours'))
          
          # convert to local time zone
          times = times.gmt + (houroffset * 3600)
          
          selected.values <- all.values[site_loc,]
          out <- data.frame(times, selected.values)
          names(out) <- c('datetime', variable)
          # write to file
          outfileName <- paste(outputLocation, sites$Land[site_num], '_',
                               lon,'_',lat,'_', file_var, '.csv', sep='')
          
          if((year == startyear) & (month == 1))
            write.table(out, outfileName, sep=',', eol=eol.val, row.names=FALSE, col.names = TRUE)
          else
            write.table(out, outfileName, sep=',', eol=eol.val, row.names=FALSE, 
                        col.names = FALSE, append=TRUE)
          
        }
      }
    } 
  }
  
  # now read in all variables for all sites and convert their values
  
  for (site_num in 1:site_count){
    # get values for each site
    site_loc <- which(land$Land == sites$Land[site_num])
    lon <- land$Longitude[site_loc]
    lat <- land$Latitude[site_loc]
    timezone <- sites$timezone[site_num]
    sunTimeOffset <- sites$SolarOffset[site_num]
  
    for (var_num in 1:num_vars){
      # get all variables and assemble dataframe
      variable <- vars[var_num]
      file_var <- file_vars[var_num]
      infileName <- paste(outputLocation, sites$Land[site_num], '_',
                           lon,'_',lat,'_', file_var, '.csv', sep='')
      var_values <- read.csv(file=infileName, header=TRUE, stringsAsFactors = FALSE)
      names(var_values)[2] <- file_var
      var_values$datetime <- as.POSIXct(var_values$datetime, format='%Y-%m-%d %H:%M', tz=timezone)
      
      if(var_num == 1)
        all_values <- var_values
      else
        all_values <- merge(all_values, var_values, by='datetime', all=TRUE)
    }
    
    # convert variables
    all_values$t <- all_values$Tair_WFD - 273.15                    # K -> °C
    all_values$PSurf_WFD <- all_values$PSurf_WFD * 0.01             # Pa -> mb
    all_values$e <- all_values$Qair_WFD * all_values$PSurf_WFD / 
      (0.378 * all_values$Qair_WFD + 0.622)
    
    all_values$ea <- all_values$e * 0.1                             # convert back to kPa
    all_values$p_GPCC <- (all_values$Rainf_WFD_GPCC + 
                            all_values$Snowf_WFD_GPCC) * 10800      # mm/s -> mm
    all_values$p_CRU <- (all_values$Rainf_WFD_CRU + 
                            all_values$Snowf_WFD_CRU) * 10800       # mm/s -> mm
    
    all_values$u <- all_values$Wind
    all_values$Qsi <- all_values$SWdown                             # added by Dhiraj
    all_values$Qli <- all_values$LWdown                             # added by Dhiraj
    
    obsData <- all_values[,c('datetime', 't', 'ea', 'u',
                             'p_GPCC', 'p_CRU', 'Qsi', 'Qli')] #dp, added Qsi and Qli
    
    # write obs file
    obs_filename <-  paste(outputLocation, sites$Land[site_num], '_',
                                          lon,'_',lat, '_', startyear,'-', endyear, 
                           '_WFD_3hour.obs', sep='')
    CRHMr::writeObsFile(obsData, obs_filename, comment='WATCH WFD')
    
    # now erase all other files
    for (var_num in 1:num_vars){
      # get all variables and assemble dataframe
      variable <- vars[var_num]
      file_var <- file_vars[var_num]
      eraseName <- paste(outputLocation, sites$Land[site_num], '_',
                          lon,'_',lat,'_', file_var, '.csv', sep='')
      file.remove(eraseName)
    }
    
    # create hourly obs
    hourly_p_GPCC <- CRHMr::distributeP(obs=obsData, p.cols=4, timestep=1)
    hourly_p_CRU <- CRHMr::distributeP(obs=obsData, p.cols=5, timestep=1)
    
    # shift their time by 2 hours
    hourly_p_GPCC$datetime <- hourly_p_GPCC$datetime + (2 * 3600)
    hourly_p_CRU$datetime <- hourly_p_CRU$datetime + (2 * 3600)
    

    # shift shortwave radiation by 2 hours
    Qsi <- obsData[,c('datetime', 'Qsi')]
    Qsi$datetime <- Qsi$datetime + (2*3600)
    hourly_Qsi <- CRHMr::distributeQsi(Qsi, QsiColnum = 1, latitude=lat, 
                                       sunTimeOffset = 0, timeStep = 1, 
                                       solarMethod = solarMethod, details=FALSE)
    # add sun time offset
    hourly_Qsi$datetime <- hourly_Qsi$datetime + (sunTimeOffset * 3600)
    
    hourly_t_ea_u <- CRHMr::distributeInst(obsData, obsCols=c(1,2,3), timeStep=1, 
                                           interpolationMethod=interpolationMethod,
                                           maxLength=7, 
                                           quiet=quiet, logfile = logfile)
    cleanData <- na.omit(obsData)
    Qli <- cleanData[,c('datetime','Qli')]
    # shift longwave radiation by 2 hours
    Qli$datetime <- Qli$datetime + (2*3600)
    
    hourly_Qli <- CRHMr::distributeQli(Qli, QliColnum = 1, tObs = hourly_t_ea_u, tColnum = 1,
                                       timeStep = 1, interpolationMethod = interpolationMethod)
    
    # merge all together
    all <- merge(hourly_t_ea_u, hourly_p_GPCC , by='datetime')
    all <- merge(all, hourly_p_CRU, by='datetime')
    all <- merge(all, hourly_Qsi, by='datetime')
    all <- merge(all, hourly_Qli, by='datetime')
    names(all) <- c('datetime', 't', 'ea', 'u', 'p_GPCC', 'p_CRU', 'Qsi', 'Qli')
    
    obs_filename <-  paste(outputLocation, sites$Land[site_num], '_',
                           lon,'_',lat,'_', startyear,'-', endyear,'_WFD_hourly.obs', sep='')
    
    # record to logfile
    comment <- paste('WATCHcreateHourlyWFDEIobs startyear:',startyear,
                     ' endyear: ', endyear, ' lon: ',' lat: ',
                     lat, ' houroffset:', houroffset, ' obsFileName:', obs_filename, sep='')
    result <- CRHMr::logAction(comment, logfile)
    
    # write to file
    obsHeader <- paste('Created from WATCH WFD reanalysis data. Lat = ',lat, 
                       ' Lon = ',lon, sep='')

    
    result <- CRHMr::writeObsFile(all, obs_filename, comment=obsHeader, quiet=quiet, logfile = logfile)
    return(result)
    
  }
}