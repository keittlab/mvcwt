\name{plot.mvcwt}
\alias{plot.mvcwt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.mvcwt(w, var = 1, scale = 1, titles = TRUE, z.fun = "Re", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{var}{
%%     ~~Describe \code{var} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{titles}{
%%     ~~Describe \code{titles} here~~
}
  \item{z.fun}{
%%     ~~Describe \code{z.fun} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
