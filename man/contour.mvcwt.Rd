\name{contour.mvcwt}
\alias{contour.mvcwt}
\title{Make contour plot of a \code{\link{mvcwt}} object}
\description{Draws a contour plot}
\usage{
\method{contour}{mvcwt}(x, z.fun = "Re", bound = 1, reset.par = TRUE, ...)
}
\arguments{
  \item{x}{an object produced by \code{\link{mvcwt}} or \code{\link{wmr}}}
  \item{z.fun}{a function to apply to the data prior to plotting}
  \item{bound}{if finite, draw boundary lines \code{bound * scale} from plot boundaries}
  \item{reset.par}{if true, reset the graphical parameters on exit}
  \item{\dots}{passed to the \code{\link{contour}} function}
}
\details{Draws a contour plot. If you want to add more plotting elements, set \code{reset.pars} to false.}
\value{The object \code{x} is returned invisibly.}
\author{Timothy H. Keitt}
\seealso{\code{\link{mvcwt}}, \code{\link{wmr}}}
\keyword{graphics}
\examples{
run.it = getOption("run.long.examples")     
if ( !is.null(run.it) && run.it )
{
  data(lrlake)
  x = subset(lrlake, Basin == "Treatment", LRL.Day) / 365.25
  y = subset(lrlake, Basin == "Treatment", -(1:8))
  w = mvcwt(x, y, min.scale = 0.25, max.scale = 4)
  mr = wmr(w)
  contour(mr)
}}
