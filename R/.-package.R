

#' Little Rock Lake zooplankton dataset
#' 
#' Contains time series for 10 dominant crustaceous species of zooplanking
#' sampled from Little Rock Lake, Wisconsin. Samples come from two basins: one
#' treated to lower pH and the other an untreated reference.
#' 
#' 
#' @name lrlake
#' @docType data
#' @format A data frame with 592 observations on the following 18 variables.
#' \describe{ \item{list("Month")}{month of year} \item{list("Day")}{day of
#' month} \item{list("Year")}{sample year} \item{list("LRL.Day")}{days from
#' start of study} \item{list("Year.Day")}{day of year}
#' \item{list("Surface.pH")}{pH measured at lake surface}
#' \item{list("Station")}{sample station} \item{list("Basin")}{lake basin,
#' either \code{Reference} or \code{Treatment}}
#' \item{list("Diacyclops")}{species data} \item{list("Mesocyclops")}{species
#' data} \item{list("Tropocyclops")}{species data}
#' \item{list("Leptodiaptomus")}{species data} \item{list("Bosminid")}{species
#' data} \item{list("D..dubia")}{species data}
#' \item{list("D..parvula")}{species data} \item{list("D..catawba")}{species
#' data} \item{list("Diaphanosoma.birgei")}{species data}
#' \item{list("Holopedium.gibberum")}{species data} }
#' @references Keitt, T. H. 2008. Coherent ecological dynamics induced by
#' large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.
#' @source http://lter.limnology.wisc.edu/
#' @keywords datasets
#' @examples
#' 
#' data(lrlake)
#' x = subset(lrlake, Basin == "Reference", LRL.Day)
#' y = subset(lrlake, Basin == "Reference", -(1:8))
#' matplot(x, y, type = "l", lty = 1)
#' 
NULL





#' Wavelet transform of multivariate time series
#' 
#' Computes continuous wavelet transform of multiple irregularly sampled time
#' series.
#' 
#' \tabular{ll}{ Package: \tab mvcwt\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2013-10-27\cr License: \tab GPL\cr } The main functions
#' are \code{\link{mvcwt}}, which computes the wavelet transform of multiple
#' time series, and \code{\link{wmr}}, which computes the wavelet modulus
#' ratio, a measure of time series coherence.
#' 
#' Note that this is a complete rewrite of the code used in the reference
#' below, and as such it is not well tested. It may give different or
#' inaccurate results. I recommend you run tests on known data.
#' 
#' The most recent development version of this package can be found at
#' \url{https://bitbucket.org/tkeitt/mvcwt/overview}.
#' 
#' @name mvcwt-package
#' @docType package
#' @author Timothy H. Keitt (\url{http://www.keittlab.org})
#' 
#' Tim Keitt <tkeitt@@gmail.com>
#' @references Keitt, T. H. 2008. Coherent ecological dynamics induced by
#' large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.
#' @keywords package
#' @examples
#' 
#' run.it = getOption("run.long.examples")     
#' if ( !is.null(run.it) && run.it )
#' {
#'   x = seq(-pi, pi, len = 200)
#'   y1 = sin(8 * x) + sin(32 * x)
#'   y2 = sin(8 * (x + pi/8)) + sin(32 * x)
#'   matplot(x, cbind(y1, y2), type = "l", lty = 1)
#'   w = mvcwt(x, cbind(y1, y2))
#'   plot(w, var = 1:2, scale = 2^seq(log2(min(w$y)), log2(max(w$y)), len = 5))
#'   mr = wmr(w, smoothing = 2)
#'   image(mr, reset.par = FALSE)
#'   contour(mr, levels = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), add = TRUE)
#' }
NULL



