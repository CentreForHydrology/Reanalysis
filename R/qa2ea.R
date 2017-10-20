#' Converts specific humidity to vapour pressure.
#'
#' @param qair Required. Specific humidity (dimensionless).
#' @param psurf Required. Surface air pressure in Pa.
#'
#' @return If successful, returns the vapour pressure in kPa. If unsuccessful, returns the value \code{FALSE}.
#' @author Kevin Shook
#' @references R code for conversion of air pressure and absolute humidity was taken from project PEcAn The Predictive Ecosystem Analyzer \url{http://pecanproject.github.io}. The source code is available at \url{https://github.com/PecanProject/pecan/blob/master/modules/data.atmosphere/R/metutils.R.}
#' @export
#'
#' @examples ea <- qa2ea(0.001, 101325)

qa2ea <- function(qair, psurf){

  # check parameter values
  if (is.null(qair) | (qair == "")) {
    cat('Error: missing humidities')
    return(FALSE)
  }
  
  if (is.null(psurf) | (psurf == "")) {
    cat('Error: missing surface pressures')
    return(FALSE)
  }
  
  psurf <- psurf * 0.01                           # Pa -> mb
  e <- qair * psurf / (0.378 * qair + 0.622)
  ea <- e * 0.1                                   # convert mb back to kPa
  return(ea)
}