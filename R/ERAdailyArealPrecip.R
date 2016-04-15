#' Calculates daily ERA areal precipitation
#' @description Calculates the daily total ERA precipitation, for all locations.
#' @param ERAarealPrecip Required. The three-hour ERA areal precipitation as returned by \code{ERAgetArealPrecip}.
#'
#' @return If unsuccessful, returns \code{FALSE}. If successful, returns a list containing the following: 
#'\tabular{ll}{
#'\bold{Name}\tab \bold{meaning}\cr
#'\code{daily_precip}\tab 3d array (lon x lat x date) of precip (mm)\cr
#'\code{lonres}\tab longitude resolution (\eqn{^\circ}{degrees})\cr
#'\code{latres}\tab latitude resolution (\eqn{^\circ}{degrees})\cr
#'\code{minLon}\tab mininum longitude of data(\eqn{^\circ}{degrees}W)\cr
#'\code{maxLon}\tab maximum longitude of data(\eqn{^\circ}{degrees}W)\cr
#'\code{minLat}\tab minimum latitude of data(\eqn{^\circ}{degrees})\cr
#'\code{minLat}\tab maximum latitude of data(\eqn{^\circ}{degrees})\cr
#'\code{date}\tab the date of each layer\cr
#' }
#' @export
#'
#' @examples \dontrun{daily_precip <- ERAdailyArealPrecip(threehour_precip)
#' }
ERAdailyArealPrecip <- function(ERAarealPrecip){
  # check parameter
  datetime <- ERAarealPrecip$datetime
  if (length(datetime == 0)){
    cat('Error: missing data\n')
    return(FALSE)
  }
  
  precip <- ERAarealPrecip$precip
  dates_all <- as.Date(datetime)
  dates_unique <- unique(dates_all)
  num_dates <- length(dates_unique)
  
  lonRes <- ERAarealPrecip$lonRes
  latRes <- ERAarealPrecip$latRes
  minLon <- ERAarealPrecip$minLon
  maxLon <- ERAarealPrecip$maxLon  
  minLat <- ERAarealPrecip$minLat
  maxLat <- ERAarealPrecip$maxLat
  
  # get dimensions of precip array
  num_lats <- nrow(precip)
  num_lons <- ncol(precip)
  
  # create data frame to hold totals
  daily <- array(data=0, dim=c(num_lats, num_lons, num_dates))
  for (day_num in 1:num_dates){
    # select all layers in each day
    selected_date <- dates_unique[day_num]
    layers <- which(dates_all == selected_date)
    selected <- precip[,,layers]
    
    # sum all layers
    layer_sum <- apply(selected,1:2,sum)
    daily[,,day_num] <- layer_sum
  }
  
  # create output list
  ret <- list(daily_precip=daily, lonRes=lonRes, latRes=latRes, minLon=minLon, maxLon=maxLon, minLat=minLat, maxLat=maxLat, date=dates_unique)
  return(ret)
  
}