\name{lrlake}
\alias{lrlake}
\docType{data}
\title{Little Rock Lake zooplankton dataset}
\description{Contains time series for 10 dominant crustaceous species of zooplanking sampled from Little Rock Lake, Wisconsin. Samples come from two basins: one treated to lower pH and the other an untreated reference.}
\usage{data(lrlake)}
\format{
  A data frame with 592 observations on the following 18 variables.
  \describe{
    \item{\code{Month}}{month of year}
    \item{\code{Day}}{day of month}
    \item{\code{Year}}{sample year}
    \item{\code{LRL.Day}}{days from start of study}
    \item{\code{Year.Day}}{day of year}
    \item{\code{Surface.pH}}{pH measured at lake surface}
    \item{\code{Station}}{sample station}
    \item{\code{Basin}}{lake basin, either \code{Reference} or \code{Treatment}}
    \item{\code{Diacyclops}}{species data}
    \item{\code{Mesocyclops}}{species data}
    \item{\code{Tropocyclops}}{species data}
    \item{\code{Leptodiaptomus}}{species data}
    \item{\code{Bosminid}}{species data}
    \item{\code{D..dubia}}{species data}
    \item{\code{D..parvula}}{species data}
    \item{\code{D..catawba}}{species data}
    \item{\code{Diaphanosoma.birgei}}{species data}
    \item{\code{Holopedium.gibberum}}{species data}
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
