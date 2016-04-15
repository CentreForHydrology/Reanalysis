WATCHgetWFDvalues <-
function(location, fname, variable, month, year, houroffset){
  # read in values from monthly WFD netCDF file
  # because the time stamps are often wrong (morons!) it is necessary
  # to ignore them and to create time stamps from scratch

  nc <- RNetCDF::open.nc(fname, write=FALSE)
  nctimes <- RNetCDF::var.get.nc(nc, variable='time', unpack = TRUE)
  all.values <- RNetCDF::var.get.nc(nc, variable==variable, unpack = TRUE)
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

  # convert to MST
  times = times.gmt + (houroffset * 3600)

  selected.values <- all.values[location,]

  # create output dataframe
  out <- data.frame(times, selected.values)
  names(out) <- c('DateTime', variable)
  return(out)
}
