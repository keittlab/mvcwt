\name{image.mvcwt}
\alias{image.mvcwt}
\title{Draw a heatmap of a \code{\link{mvcwt}} object}
\description{Draws one or more heatmaps}
\usage{
image.mvcwt(w, z.fun = "Re", bound = 1, reset.par = TRUE, ...)
}
\arguments{
  \item{w}{an object as returned by \code{\link{mvcwt}}}
  \item{z.fun}{a function applied to the data before plotting}
  \item{bound}{if finite, draw lines \code{bound * scale} units inside the plot boundaries}
  \item{reset.par}{if true, reset graphical parameters on exit}
  \item{\dots}{additional arguments passed to \code{\link{image}}}
}
\details{This function will draw a series of heatmaps on the graphical device. If you want to add additional graphical elements, set \code{reset.par} to false.}
\value{\code{w} is returned invisibly}
\author{Timothy H. Keitt}
\seealso{\code{\link{mvcwt}}, \code{\link{wmr}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (w, z.fun = "Re", bound = 1, reset.par = TRUE, ...) 
{
    require(grDevices)
    require(RColorBrewer)
    z.fun = match.fun(z.fun)
    opar = par(no.readonly = TRUE)
    if (reset.par) 
        on.exit(par(opar))
    pal = colorRampPalette(rev(brewer.pal(11, "Spectral")))(1024)
    with(w, {
        nvar = ifelse(length(dim(z)) == 3, dim(z)[3], 1)
        par(mfrow = c(nvar, 1), mar = rep(0.2, 4), oma = rep(5, 
            4))
        for (i in 1:nvar) {
            image(x, y, z.fun(z[, , i]), log = "y", col = pal, 
                axes = FALSE, ...)
            if (i\%\%2) 
                axis(2)
            else axis(4)
            if (is.finite(bound)) {
                lines(min(x) + bound * y, y, lty = 2, lwd = 2, 
                  col = "darkgrey")
                lines(max(x) - bound * y, y, lty = 2, lwd = 2, 
                  col = "darkgrey")
            }
            box()
        }
        axis(1)
        mtext("Location", 1, 3, outer = TRUE)
        mtext("Scale", 2, 3, outer = TRUE)
    })
    return(invisible(w))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphics}
