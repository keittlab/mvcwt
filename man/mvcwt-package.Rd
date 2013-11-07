\name{mvcwt-package}
\alias{mvcwt-package}
\docType{package}
\title{Wavelet transform of multivariate time series}
\description{
Computes continuous wavelet transform of multiple irregularly sampled time series.
}
\details{
\tabular{ll}{
Package: \tab mvcwt\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-27\cr
License: \tab GPL\cr
}
The main functions are \code{\link{mvcwt}}, which computes the wavelet transform of multiple time series, and \code{\link{wmr}}, which computes the wavelet modulus ratio, a measure of time series coherence.

Note that this is a complete rewrite of the code used in the reference below, and as such it is not well tested. It may give different or inaccurate results. I recommend you run tests on known data.}
\author{Timothy H. Keitt

Tim Keitt <tkeitt@gmail.com>
}
\references{Keitt, T. H. 2008. Coherent ecological dynamics induced by large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.}
\keyword{ package }
\examples{
run.it = getOption("run.long.examples")     
if ( !is.null(run.it) && run.it )
{
  x = seq(-pi, pi, len = 200)
  y1 = sin(8 * x) + sin(32 * x)
  y2 = sin(8 * (x + pi/8)) + sin(32 * x)
  matplot(x, cbind(y1, y2), type = "l", lty = 1)
  w = mvcwt(x, cbind(y1, y2))
  plot(w, var = 1:2, scale = 2^seq(log2(min(w$y)), log2(max(w$y)), len = 5))
  mr = wmr(w, smoothing = 2)
  image(mr, reset.par = FALSE)
  contour(mr, levels = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), add = TRUE)
}}
