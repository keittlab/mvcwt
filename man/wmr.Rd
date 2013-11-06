\name{wmr}
\alias{wmr}
\title{Compute the wavelet modulus ratio of multivariate data}
\description{Computes the wavelet modulus ratio described in Keitt (2008). A value of one indicated perfect synchrony among time series and a value of zero, perfect compensation.}
\usage{wmr(w, smoothing = 1)}
\arguments{
  \item{w}{an object such as returned by \code{\link{mvcwt}}}
  \item{smoothing}{width of smoothing kernel; larger values give more smoothing}
}
\value{an object of class "mvcwt"}
\references{Keitt, T. H. 2008. Coherent ecological dynamics induced by large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.}
\author{Timothy H. Keitt}
\seealso{\code{\link{mvcwt}}, \code{\link{image.mvcwt}}}
\examples{
run.it = getOption("run.long.examples")     
if ( !is.null(run.it) && run.it )
{
  data(lrlake)
  x = subset(lrlake, Basin == "Treatment", LRL.Day) / 365.25
  y = subset(lrlake, Basin == "Treatment", -(1:8))
  w = mvcwt(x, y, min.scale = 0.25, max.scale = 4)
  mr = wmr(w)
  image(mr, reset.par = FALSE)
  contour(mr, bound = NA, add = TRUE)
}}
\keyword{wavelets}
