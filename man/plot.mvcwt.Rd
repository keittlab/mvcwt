\name{plot.mvcwt}
\alias{plot.mvcwt}
\title{Plot wavelet output}
\description{Plot multivariate wavelet output across variables, scales or both.}
\usage{
plot.mvcwt(w, var = 1, scale = 1, titles = TRUE, z.fun = "Re", ...)
}
\arguments{
  \item{w}{an object such as produced by \code{\link{mvcwt}}}
  \item{var}{which variables to plot; can be a vector}
  \item{scale}{which scales to plot; can be a vector; closest scale is picked}
  \item{titles}{plot titles on each sub-plot?}
  \item{z.fun}{apply function to data prior to plotting}
  \item{\dots}{additional graphical parameters passed to \code{\link{plot}}}
}
\details{Makes one or more plots on the graphics device. Total number of plots is limited to 10.}
\value{\code{w} is returned invisibly}
\author{Timothy H. Keitt}
\seealso{\code{\link{mvcwt}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (w, var = 1, scale = 1, titles = TRUE, z.fun = "Re", 
    ...) 
{
    z.fun = match.fun(z.fun)
    opar = par(no.readonly = TRUE)
    on.exit(par(opar))
    with(w, {
        nvar = length(var)
        nscal = length(scale)
        par(mfrow = c(nvar * nscal, 1), mar = rep(0.6, 4), oma = rep(5, 
            4), xpd = NA)
        for (s in scale) {
            scale.i = which.min(abs(y - s))
            for (i in 1:nvar) {
                z.out = z.fun(z[, scale.i, var[i]])
                plot(x, z.out, xlab = NA, ylab = NA, axes = FALSE, 
                  type = "l", lwd = 2, ...)
                lines(range(x), rep(median(z.out), 2), lty = 2)
                if (titles) {
                  tstr = paste("Scale =", signif(y[scale.i], 
                    2), "Var =", var[i], sep = " ")
                  title(main = tstr, line = 0.25)
                }
                if (i\%\%2) 
                  axis(2)
                else axis(4)
            }
        }
        axis(1)
        mtext("Location", 1, 3, outer = TRUE)
        mtext("Value", 2, 3, outer = TRUE)
    })
    invisible(w)
  }
}
\keyword{graphics}
