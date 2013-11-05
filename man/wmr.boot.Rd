\name{wmr.boot}
\alias{wmr.boot}
\title{Boot strap p-values for wavelet modulus ratio}
\description{Performs a phase-randomization bootstrap estimate of the null hypothesis of independent time series}
\usage{
wmr.boot(w, smoothing = 1, reps = 1000, mr.func = "wmr")
}
\arguments{
  \item{w}{an object such as returned by \code{\link{mvcwt}}}
  \item{smoothing}{degree of smoothing; larger values give greater smoothing}
  \item{reps}{number of repetitions}
  \item{mr.func}{a function taking a "mvcwt" object to be applied to each trial}
}
\details{The phases are randomized \code{reps} times for each combination of input variable and scale. This package depends heavily on the \code{dopar} function in the \code{foreach} package. If you do not have a lot of cores available to you, you may need to let this run overnight.}
\value{an object of class "mvcwt" suitable for use with \code{\link{contour.mvcwt}}.}
\references{Keitt, T. H. 2008. Coherent ecological dynamics induced by large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.}
\author{Timothy H. Keitt}
\seealso{\code{\link{mvcwt}}, \code{\link{wmr}}}
\keyword{statistics}
