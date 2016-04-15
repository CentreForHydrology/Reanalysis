#' Extracts the WATCH WFD precipitation for a specified domain
#'
#' @description Extacts the WATCH WFD precipitation (rainfall or snowfall) from a monthly netCDF file for a specified domain.
#' @param ncdfFile Required. NetCDF file containing monthly WATCH WFD precipitation (rainfall or snowfall). The file name must begin with the precipitation type (i.e. \option{Rainf} or \option{Snowf} and must end with the year and month, as in \option{190101}.
#' @param minLon Required. Minimum longitude for area to be extracted.
#' @param maxLon Required. Maximum longitude for area to be extracted.
#' @param minLat Required. Minimum latitude for area to be extracted.
#' @param maxLat Required. Maximum latitude for area to be extracted.
#' @param logfile Optional. Name of the file to be used for logging the action. Normally not used 
#' @return If successful, returns a data frame containing all of the precipitation values. The rows contain the values for each location for all time intervals. The first 2 columns contain the Longitude and Latitude, respectively for each location within the domain. The following columns contain the precipitation rate for each time interval. The names of these columns are the time intervals, e.g. \option{1901-01-01_00:00} (note the underscore) in GMT.
#' @author Kevin Shook
#' @export
#'
#' @examples \dontrun{
#' monthly_rain <- WATCHgetWFDarealPrecip('Rainf_WFD_CRU_190101.nc',
#'  minLon=-105, maxLon=-95, minLat = 49, maxLat = 56)}
WATCHgetWFDarealPrecip <- function(ncdfFile='', minLon=0, maxLon=0, 
                                   minLat=0, maxLat=0, logfile=''){
  
  ## check parameters
  if (ncdfFile == ''){
    cat('Error: must specify netCDF file\n')
    return(FALSE)
  }
  
  if (minLon == ''){
    cat('Error: missing output minimum longitude\n')
    return(FALSE)
  }
  if (maxLon == ''){
    cat('Error: missing output maximum longitude\n')
    return(FALSE)
  }
  if (minLat == ''){
    cat('Error: missing output minimum latitude\n')
    return(FALSE)
  }
  if (maxLat == ''){
    cat('Error: missing output maximum latitude\n')
    return(FALSE)
  }
  
  # split up file name to extract year, month and precipitation variable
  
  ncdfFile_len <- nchar(ncdfFile)
  month_end <- ncdfFile_len - 3
  month_start <- month_end - 1
  year_end <- month_start - 1
  year_start <- year_end - 3
  year <- substr(ncdfFile, year_start, year_end)
  month <- substr(ncdfFile, month_start, month_end)
  
  base <- basename(ncdfFile)
  precip_var <- substr(base, 1, 5)
  
  # read in entire netCDF file 
  nc <- RNetCDF::open.nc(ncdfFile, write=FALSE)
  # don't use time variable as it is often wrong
  # nctimes <- RNetCDF::var.get.nc(nc, variable='time', unpack = TRUE)
  all.values <- RNetCDF::var.get.nc(nc, variable=precip_var, unpack = TRUE)
  RNetCDF::close.nc(nc)
  
  # figure out time interval
  time_count <- ncol(all.values)
  if (time_count > 125)
    interval <- 3
  else
    interval <- 6

  # get specified range
  land <- land
  land_selected <- land[(land$Longitude >= minLon) &
                        (land$Longitude <= maxLon) &
                        (land$Latitude >= minLat) &
                        (land$Latitude <= maxLat) ,]
 
  selected_values <- all.values[land_selected$Land,] * 
    interval * 3600 # convert from mm/s to mm
  

  # generate datetime sequence
  origin.string <- paste(year,'-',month,'-01 00:00:00',sep='')
  origin.datetime <- timeDate::timeDate(origin.string, format = "%Y-%m-%d %H:%M:%S")
  
  # create time sequence
  times.gmt <- timeDate::timeSequence(from = origin.datetime, by=interval*3600,
                                      length.out=time_count,
                                      format="%Y-%m-%d %H:%M:%S", zone='UTC')

  # add lat and lon columns to output
  precip <- data.frame(land_selected$Longitude, land_selected$Latitude, selected_values)
  datetimes <- format(times.gmt, format="%Y-%m-%d_%H:%M")
  names(precip)[3:(ncol(precip))] <- datetimes
  names(precip)[1:2] <- c('Longitude', 'Latitude')
  
  comment <- paste('WATCgetWFDarealPrecip netCDF:', ncdfFile,
                   ' minLon:', minLon,
                   ' maxLon:', maxLon,
                   ' minLat:', minLat,
                   ' maxLat:', maxLat,
                   sep='')  
  
  result <- CRHMr::logAction(comment, logfile)
  
  if (result)
    return (precip)
  else
    return(result)
}