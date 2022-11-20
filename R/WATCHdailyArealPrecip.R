#' Calculates daily WATCH areal precipitation
#' @details Calculates the daily total WATCH precipitation for all locations in a region.
#' @param monthlyPrecip All WATCH precipitation for a month, as returned by \code{WATCHgetWFDarealPrecip}.
#' @param timezone Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.
#'
#' @return If successful, returns a data frame consisting if the date,longitude, latitude, and total precipitation (in mm) for each point. If unsuccessful, returns \code{FALSE}.
#' @author Kevin Shook
#' @export
#'
#' @examples \dontrun{
#' daily_rain <- WATCHdailyArealPrecip(monthly_rain, 'CST')}
WATCHdailyArealPrecip <- function(monthlyPrecip, timezone) {
  # check parameters
  if (nrow(monthlyPrecip) == 0) {
    cat('Error: missing precip values\n')
    return(FALSE)
  }
  
  # calculate hour offset between timezone and GMT
  datetime <- '2000-01-01 01:00'
  s1 <- as.POSIXct(datetime, format = '%Y-%m-%d %H:%M', tz = 'UTC')
  s2 <- as.POSIXct(datetime, format = '%Y-%m-%d %H:%M', tz = timezone)
  houroffset <- as.numeric(difftime(s2, s1, units = 'hours'))
  lonlat <- monthlyPrecip[, c(1, 2)]
  monthlyPrecip <- monthlyPrecip[, 3:ncol(monthlyPrecip)]
  monthlyPrecip$ID <- seq(1:nrow(monthlyPrecip))
  lonlat$ID <- monthlyPrecip$ID
  
  # melt dataframe
  melted <- reshape2::melt(monthlyPrecip, id.vars = 'ID')
  names(melted)[2:3] <- c('datetime', 'precip')
  melted$datetime <- as.character(melted$datetime)
  
  datetimes <-
    as.POSIXct(melted$datetime , format = '%Y-%m-%d_%H:%M', tz = 'UTC')
  
  # convert to local time zone
  datetimes <- datetimes + (houroffset * 3600)
  melted$datetime <- lubridate::force_tz(datetimes, tzone = timezone)
  
  # convert datetimes to dates
  melted$date <- as.Date(melted$datetime)
  
  # get daily totals
  daily <-
    aggregate(melted$precip,
              by = list(melted$ID, melted$date),
              FUN = 'sum')
  names(daily) <- c('ID', 'date', 'precip')
  
  # add longitude and latitude
  daily <- merge(daily, lonlat, by = 'ID')
  
  # select variables to be returned
  daily <- daily[, c('date', 'Longitude', 'Latitude', 'precip')]
  names(daily)[2:3] <- c('longitude', 'latitude')
  
  # sort by date
  daily <- daily[order(daily$date), ]
  
  return(daily)
}
