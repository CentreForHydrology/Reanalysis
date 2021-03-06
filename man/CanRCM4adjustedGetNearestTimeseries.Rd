% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CanRCM4adjustedGetNearestTimeseries.R
\name{CanRCM4adjustedGetNearestTimeseries}
\alias{CanRCM4adjustedGetNearestTimeseries}
\title{Reads a time series for the nearest point from adjusted CanRCM4}
\usage{
CanRCM4adjustedGetNearestTimeseries(
  netCDFfile = "",
  longitude = 0,
  latitude = 0,
  startDate = "1979-01-01",
  endDate = "2100-12-01",
  timezone = "",
  logfile = ""
)
}
\arguments{
\item{netCDFfile}{Required. The name of a NetCDF file containing a single variable.}

\item{longitude}{Required. The longitude of the point being sought. Valid values appear to be between -90 and -142, but the input value is \emph{not} checked for validity, in case the model extent changes.}

\item{latitude}{Required. The latitude of the point being sought. Valid values appear to be between 45 and 75, but the input value is \emph{not} checked for validity, in case the model extent changes.}

\item{startDate}{Optional. Beginning date of data to be extracted. A string formatted as "yyyy-mm-dd". The default value of \option{1979-01-01} is the beginning of the data.}

\item{endDate}{Optional. Ending date of data to be extracted. A string formatted as "yyyy-mm-dd". The default value of \option{2100-12-01} is the end of the data.}

\item{timezone}{Required. The name of the timezone of the data as a character string. This should be the timezone of your data, but omitting daylight savings time. Note that the timezone code is specific to your OS. To avoid problems, you should use a timezone without daylight savings time. Under Linux, you can use \option{CST} and \option{MST} for Central Standard or Mountain Standard time, respectively. Under Windows or OSX, you can use \option{etc/GMT+6} or \option{etc/GMT+7} for Central Standard and Mountain Standard time. DO NOT use \option{America/Regina} as the time zone, as it includes historical changes between standard and daylight savings time.}

\item{logfile}{Optional. Name of the file to be used for logging the action. Normally not used.}
}
\value{
If successful, returns a standard \pkg{CRHMr} data frame containing the datetime and the variable. If unsuccessful, returns the value \code{FALSE}.
}
\description{
This function extracts a time series of 3-hourly values from a NetCDF file of CanRCM4 data, which has been bias-corrected using the WFDEI 3-hour reanalysis values. The values are stored at a spatial resolution of 0.125 degees, so the closest point to the specified location will be used. Note that because the reanalysis data omit February 29, the values returned by this function will have \code{NA} values for leap days. You will have to fill these values yourself.\cr Each NetCDF file contains a single variable. Typically the first letters of the file name designate the variable. The variables are converted to values, and variable names, appropriate for CRHM when they are extracted.
\tabular{lllr}{
\bold{NetCDF parameter} \tab \bold{netCDF units} \tab \bold{CRHM Variable} \tab \bold{CRHM units}\cr
pr - 3hr precip\tab mm/s \tab p \tab mm\cr
tas - surface air temp. \tab K \tab t \tab  \eqn{^\circ}{}C\cr
huss - specific humidity \tab dimensionless \tab qair \tab dimensionless \cr
sfcWind - surface (10m) wind speed \tab m/s \tab u10 \tab m/s \cr
ps - surface pressure \tab Pa \tab ps \tab Pa \cr
rsds - incoming SW radiation \tab W/\eqn{^2}{^2}  \tab Qsi \tab W/\eqn{^2}{^2} \cr
rlds - incoming LW radiation \tab W/\eqn{^2}{^2}  \tab Qli \tab W/\eqn{^2}{^2} \cr
}
}
\examples{
\dontrun{ f <- "ps_CanRCM4_hist+fut_1979_2100.nc4"
r <- RCM4adjustedGetNearestTimeseries(f, longitude=-115.27, 
latitude=52.03,  startDate = "1980-01-01", endDate="1980-12-31")
}
}
