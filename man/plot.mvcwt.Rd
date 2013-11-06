\name{plot.mvcwt}
\alias{plot.mvcwt}
\title{Plot wavelet output}
\description{Plot multivariate wavelet output across variables, scales or both.}
\usage{
\method{plot}{mvcwt}(x, var = 1, scale = 1, titles = TRUE, z.fun = "Re", ...)
}
\arguments{
  \item{x}{an object such as produced by \code{\link{mvcwt}}}
  \item{var}{which variables to plot; can be a vector}
  \item{scale}{which scales to plot; can be a vector; closest scale is picked}
  \item{titles}{plot titles on each sub-plot?}
  \item{z.fun}{apply function to data prior to plotting}
  \item{\dots}{additional graphical parameters passed to \code{\link{plot}}}
}
\details{Makes one or more plots on the graphics device. Total number of plots is limited to 10.}
\value{\code{x} is returned invisibly}
\author{Timothy H. Keitt}
\seealso{\code{\link{mvcwt}}}
\examples{
run.it = getOption("run.long.examples")     
if ( !is.null(run.it) && run.it )
{
  data(lrlake)
  x = subset(lrlake, Basin == "Reference", LRL.Day) / 365.25
  y = subset(lrlake, Basin == "Reference", -(1:8))
  w = mvcwt(x, y, min.scale = 0.25, max.scale = 4)
  plot(w, var = 1:10)
}}
\keyword{graphics}
