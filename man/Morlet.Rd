\name{Morlet}
\alias{Morlet}
\title{The Morlet function}
\description{Given a sequence of lag distances, this function returns the Morlet wavelet}
\usage{Morlet(lag)}
\arguments{
  \item{lag}{A sequence of lag distances, typically a matrix}
\details{This version of the Morlet is scaled so that the central frequency is exactly 2Pi radians. This is the simple version of the Morlet, sometimes referred to as a psuedo-wavelet as it it not precisely normalized, leading to some leakage into the DC component. It is therefore unsuited to reconstruction using the inverse transform.}

\value{A set of Morlet filters, typically as a matrix}

\author{Timothy H. Keitt}

\seealso{\code{\link{mvcwt}}}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (lag) 
{
    exp(-lag^2/2 + (0+2i) * pi * lag)/pi^4
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{wavelets}
