\name{lrlake}
\alias{lrlake}
\docType{data}
\title{Little Rock Lake zooplankton dataset}
\description{Contains time series for 10 dominant crustaceous species of zooplanking sampled from Little Rock Lake, Wisconsin. Samples come from two basins: one treated to lower pH and the other an untreated reference.}
\usage{data(lrlake)}
\format{
  A data frame with 592 observations on the following 18 variables.
  \describe{
    \item{\code{Month}}{a numeric vector}
    \item{\code{Day}}{a numeric vector}
    \item{\code{Year}}{a numeric vector}
    \item{\code{LRL.Day}}{a numeric vector}
    \item{\code{Year.Day}}{a numeric vector}
    \item{\code{Surface.pH}}{a numeric vector}
    \item{\code{Station}}{a numeric vector}
    \item{\code{Basin}}{a factor with levels \code{Reference} \code{Treatment}}
    \item{\code{Diacyclops}}{a numeric vector}
    \item{\code{Mesocyclops}}{a numeric vector}
    \item{\code{Tropocyclops}}{a numeric vector}
    \item{\code{Leptodiaptomus}}{a numeric vector}
    \item{\code{Bosminid}}{a numeric vector}
    \item{\code{D..dubia}}{a numeric vector}
    \item{\code{D..parvula}}{a numeric vector}
    \item{\code{D..catawba}}{a numeric vector}
    \item{\code{Diaphanosoma.birgei}}{a numeric vector}
    \item{\code{Holopedium.gibberum}}{a numeric vector}
  }
}
\source{http://lter.limnology.wisc.edu/}
\references{Keitt, T. H. 2008. Coherent ecological dynamics induced by large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.}
\examples{
data(lrlake)
x = subset(lrlake, Basin == "Reference", LRL.Day)
y = subset(lrlake, Basin == "Reference", -(1:8))
matplot(x, y, type = "l", lty = 1)
}
\keyword{datasets}
