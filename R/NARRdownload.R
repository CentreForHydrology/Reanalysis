#' Downloads NARR NetCDF files
#' @name NARRdownload
#' @description Downloads NARR NetCDF files for a specified variable for a specified time step and range of years.
#' @param interval Optional. Time interval of NARR data. Can be '3h' or 'daily', which is the default.
#' @param startYear Optional. The first year to download. Default is 1979.
#' @param endYear Optional. The last year to download. Default is 1979.
#' @param destination Optional. The destination directory for the downloaded files. The default is the current directory.
#' @param variable Optional. The variable to be downloaded. Acceptable values are 'precip' (the default), 'temp', 'rh', 'wind', 'qli' and 'qli'.
#' @param quiet Optional. Suppresses display of messages, except for errors. Because this function can be very slow to execute, the default value is \code{FALSE}, to provide information on the downloading.
#' @return Writes the specified files to the destination directory. If successful, returns \code{TRUE}. If unsuccessful, returns {FALSE}.
#' @export
#' @examples \dontrun{
#' NARRdownload()}
NARRdownload <- function(interval='daily',
                    startYear = 1979,
                    endYear = 1979,
                    destination='.',
                    variable = 'precip',
                    quiet=FALSE){
  # check parameters

  if(!file.exists(destination)){
    cat('Error: destination directory does not exist\n')
    return(FALSE)
  }

  interval <- stringr::str_to_lower(interval)

  if (stringr::str_detect(interval, stringr::fixed('da')))
    url <- 'ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/'
  else if (stringr::str_detect(interval, stringr::fixed('h')))
    url <- 'ftp://ftp.cdc.noaa.gov/Datasets/NARR/monolevel/'
  else{
    cat('Error: interval not recognized\n')
    return(FALSE)
  }

  variable <- stringr::str_to_lower(variable)
  if (stringr::str_detect(variable, stringr::fixed('pre')))
    variable <- 'apcp'
  else if (stringr::str_detect(variable, stringr::fixed('te')))
    variable <- 'air.2m'
  else if (stringr::str_detect(variable, stringr::fixed('h')))
    variable <- 'rhum'
  else if (stringr::str_detect(variable, stringr::fixed('wi')))
    variable <- 'wind'
  else if (stringr::str_detect(variable, stringr::fixed('qsi')))
    variable <- 'dswrf'
  else if (stringr::str_detect(variable, stringr::fixed('qli')))
    variable <- 'dlwrf'
  else{
    cat('Error: variable not recognized\n')
    return(FALSE)
  }

  for (yearnum in startYear:endYear){
    if (!quiet)
      cat('year:', yearnum, '\n', sep='')


    download.file(url, destination, method='wget', quiet = quiet, mode = "w",
                 cacheOK = TRUE,
                 extra = getOption("download.file.extra"))
  }

  return(TRUE)

}
