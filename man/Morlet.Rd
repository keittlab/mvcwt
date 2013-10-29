\name{Morlet}
\alias{Morlet}
\title{The Morlet function}
\description{Given a sequence of lag distances, this function returns the Morlet wavelet}
\usage{Morlet(lag)}
\arguments{
  \item{lag}{A sequence of lag distances, typically a matrix}
  }
\details{This version of the Morlet is scaled so that the central frequency is exactly 2Pi radians. This is the simple version of the Morlet, sometimes referred to as a psuedo-wavelet as it it not precisely normalized, leading to some leakage into the DC component. It is therefore unsuited to reconstruction using the inverse transform.}
\value{A set of Morlet filters, typically as a matrix}
\author{Timothy H. Keitt}
\seealso{\code{\link{mvcwt}}}
\examples{
  x = seq(-pi, pi, len = 256)
  plot(x, Re(Morlet(x)), col = "darkblue", type = "l")
  lines(x, Im(Morlet(x)), col = "darkred")
  lines(range(x), rep(0, 2), lty = 2)
}
\keyword{wavelets}
