#' Downloads NARR NetCDF files
#' @description Downloads NARR NetCDF files for a specified variable for a specified time step and range of years.
#' @param interval Optional. Time interval of NARR data. Can be '3h' or 'daily', which is the default.
#' @param startYear Optional. The first year to download. Default is 1979.
#' @param endYear Optional. The last year to download. Default is 1979.
#' @param destination Optional. The destination directory for the downloaded files. The default is the current directory.
#' @param variable Optional. The variable to be downloaded. Acceptable values are one of "acpcp", "air", "albedo", "apcp",
#'  "bgrun", "bmixl", "cape", "ccond", "cdcon", "cdlyr", "cfrzr", "cicep", "cin", "cnwat", "crain", "csnow", "dlwrf", 
#'  "dpt", "dswrf", "evap", "gflux", "hcdc", "hgt", "hlcy", "hpbl", "lcdc", "lftx4", "lhtfl", "mcdc", "mconv", 
#'  "mslet", "mstav", "pevap", "pottmp", "pr_wtr", "prate", "pres", "prmsl", "rcq", "rcs", "rcsol", "rct", "rhum",
#'  "shtfl", "shum", "snod", "snohf", "snom", "snowc", "soilm", "ssrun", "tcdc", "tke", "ulwrf", "ustm", "uswrf",
#'  "uwnd", "veg", "vis", "vstm", "vvel", "vwnd", "vwsh", "wcconv", "wcinc", "wcuflx", "wcvflx", "weasd", "wvconv",
#'  "wvinc", "wvuflx", "wvvflx".  Default is "acpcp" (precipitation).  See https://www.esrl.noaa.gov/psd/data/gridded/data.narr.monolevel.html
#'  for variables description.
#' @param quiet Optional. Suppresses display of messages, except for errors. Because this function can be very slow to execute, the default value is \code{FALSE}, to provide information on the downloading.
#' @return Writes the specified files to the destination directory. If successful, returns \code{TRUE}. If unsuccessful, returns {FALSE}.
#' @export
#' @examples \dontrun{
#' NARRdownloadNetCDF('3h', 1979, 2015)}
NARRdownloadNetCDF <- function(interval='daily',
                    startYear = 1979,
                    endYear = 1979,
                    destination='.',
                    variable = 'acpcp',
                    quiet=FALSE){
  # check parameters

  if (!file.exists(destination)) {
    cat('Error: destination directory does not exist\n')
    return(FALSE)
  }
  
  # check destination directory for terminal slash
  n <- nchar(destination) 
  last_char <- substr(destination, n, n) 
  if (last_char != '/')
    destination <- paste(destination, '/', sep = '')
  

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

  for (yearnum in startYear:endYear) {
    if (!quiet)
      cat('year:', yearnum, '\n', sep = '')

    # create url for downloading
    file_url <- paste(url, variable,'.', yearnum, '.nc', sep = '')
    destination_file <- paste(destination, variable,'.',yearnum,'.nc', sep = '')
    utils::download.file(file_url, destination_file, method = 'auto', 
                         quiet = quiet, mode = "wb", cacheOK = TRUE, 
                         extra = getOption("download.file.extra"))
  }

  return(TRUE)

}
