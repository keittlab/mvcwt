\name{wmr.boot}
\alias{wmr.boot}
\title{Boot strap p-values for wavelet modulus ratio}
\description{Performs a phase-randomization bootstrap estimate of the null hypothesis of independent time series}
\usage{
wmr.boot(w, smoothing = 1, reps = 1000)
}
\arguments{
  \item{w}{an object such as returned by \code{\link{mvcwt}}}
  \item{smoothing}{degree of smoothing; larger values give greater smoothing}
  \item{reps}{number of repetitions}
}
\details{The phases are randomized \code{reps} times for each combination of input variable and scale. This package depends heavily on the \code{dopar} function in the \code{foreach} package. If you do not have a lot of cores available to you, you may need to let this run overnight.}
\value{an object of class "mvcwt" suitable for use with \code{\link{contour.mvcwt}}.}
\references{Keitt, T. H. 2008. Coherent ecological dynamics induced by large-scale disturbance. Nature 454:331–4. doi: 10.1038/nature06935.}
\author{Timothy H. Keitt}
\seealso{\code{\link{mvcwt}}, \code{\link{wmr}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (w, smoothing = 1, reps = 1000) 
{
    require(foreach)
    mr.obs = wmr(w, smfac = smfac)
    with(w, {
        nloc = length(x)
        nvars = dim(z)[3]
        nscales = length(y)
        exports = c("reps", "wmr", "lagMat", "regularize", "seqn", 
            "Gauss", "mr.obs", "nscales", "smfac")
        flibs = c("foreach")
        mr.obs$z = foreach(i = 1:nscales, .combine = c, .export = exports, 
            .packages = flibs) \%dopar\% {
            mr.boot = foreach(j = 1:reps, .combine = cbind, .inorder = FALSE) \%do\% 
                {
                  rphase = t(array(runif(nvars, -pi, pi), dim = c(nvars, 
                    nloc)))
                  zp = z[, i, , drop = FALSE] * complex(argument = rphase)
                  as.vector(wmr(list(x = x, y = y[i], z = zp), 
                    smoothing = smoothing)$z)
                }
            foreach(j = 1:nloc, .combine = c) \%do\% {
                ecdf(mr.boot[j, ])(mr.obs$z[j, i, ])
            }
        }
        dim(mr.obs$z) = c(length(x), length(y), 1)
        return(mr.obs)
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{statistics}
