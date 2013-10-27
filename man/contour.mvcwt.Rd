\name{contour.mvcwt}
\alias{contour.mvcwt}
\title{Make contour plot of a \code{\link{mvcwt}} object}
\description{Draws a contour plot}
\usage{
contour.mvcwt(w, z.fun = "Re", bound = 1, reset.par = TRUE, ...)
}
\arguments{
  \item{w}{an object produced by \code{\link{mvcwt}} or \code{\link{wmr}}}
  \item{z.fun}{a function to apply to the data prior to plotting}
  \item{bound}{if finite, draw boundary lines \code{bound * scale} from plot boundaries}
  \item{reset.par}{if true, reset the graphical parameters on exit}
  \item{\dots}{passed to the \code{\link{contour}} function}
}
\details{Draws a contour plot. If you want to add more plotting elements, set \code{reset.pars} to false.}
\value{The object \code{w} is returned invisibly.}
\author{Timothy H. Keitt}
\seealso{\code{\link{mvcwt}}, \code{\link{wmr}}}
\keyword{graphics}
