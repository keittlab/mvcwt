\name{wmr}
\alias{wmr}
\title{Compute the wavelet modulus ratio of multivariate data}
\description{Computes the wavelet modulus ratio described in Keitt (1998). A value of one indicated perfect synchrony among time series and a value of zero, perfect compensation.}
\usage{wmr(w, smoothing = 1)}
\arguments{
  \item{w}{an object such as returned by \code{\link{mvcwt}}}
  \item{smoothing}{width of smoothing kernel; larger values give more smoothing}
}
\value{an object of class "mvcwt"}
\references{Keitt, T. H. 2008. Coherent ecological dynamics induced by large-scale disturbance. Nature 454:331–4. doi: 10.1038/nature06935.}
\author{Timothy H. Keitt}
\seealso{\code{\link{mvcwt}}, \code{\link{image.mvcwt}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (w, smoothing = 1) 
{
    with(w, {
        mods = Mod(rowSums(z, dims = 2))
        smod = rowSums(Mod(z), dims = 2)
        lmat = lagMat(x)
        exports = c("Gauss", "smfac")
        foreach(i = 1:length(y), .export = exports) \%dopar\% {
            kern = Gauss(lmat/y[i]/smoothing)
            mods[, i] = kern \%*\% mods[, i]
            smod[, i] = kern \%*\% smod[, i]
        }
        modrat = mods/smod
        i = which(abs(mods - smod) < sqrt(.Machine$double.eps))
        modrat[i] = 1
        dim(modrat) = c(length(x), length(y), 1)
        structure(list(x = x, y = y, z = modrat), class = "mvcwt")
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{wavelets}
