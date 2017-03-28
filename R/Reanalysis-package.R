#' @title Contains functions to download and process reanalysis data
#' @docType package
#' @name Reanalysis-package
#'
#' @description
#' The package contains functions to download and process reanalysis data. Functions are provided for ETA, NARR and WATCH datasets. Because there are so many functions and all of the datasets are slightly different, each function has the prefix of its dataset.
#' @references
#' To cite \pkg{Reanalysis} in publications, use the command \code{citation('Reanalysis')} to get the current version of the citation.\cr
#' The CRHM program is described in:\cr
#'\cite{Pomeroy, John W, D M Gray, T Brown, N Hedstrom, W L Quinton, R J Granger, and S K Carey. 2007. \dQuote{The Cold Regions Hydrological Model : A Platform for Basing Process Representation and Model Structure on Physical Evidence}. Hydrological Processes 21 (19): 2650-2567.}\cr
#'The CRHM model may be downloaded from \url{http://www.usask.ca/hydrology/CRHM.php}.\cr
#' @importFrom stats aggregate na.omit
#' @importFrom utils data download.file read.csv read.table write.table
#' @importFrom CRHMr interpolate saturatedVP CRHM_summary logAction win.eol distributeQsi distributeQli vp2rh
#' @importFrom stringr str_detect str_to_lower fixed
#' @importFrom timeDate timeDate timeSequence
#' @importFrom zoo na.approx
#' @importFrom lubridate force_tz
#' @importFrom proj4 project
#' @importFrom RNetCDF open.nc var.get.nc close.nc
#' @importFrom reshape2 melt
#'
NULL
