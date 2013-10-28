\name{mvcwt}
\alias{mvcwt}
\title{Computes the wavelet transform of a multivariate time series}
\description{This function takes set a set of seqences as columns in a matrix and computes the continuous wavelet transform on each.}
\usage{
mvcwt(x, y,
      scale.trim = 1,
      scale.exp = 0.5,
      nscales = length(unlist(x)),
      min.scale = 2 * scale.trim * median(diff(unlist(x))),
      max.scale = diff(range(unlist(x)))/scale.trim/2,
      scales = log2Bins(min.scale, max.scale, nscales),
      loc = regularize(as.vector(unlist(x))),
      wave.fun = "Morlet")
}
\arguments{
  \item{x}{sample locations}
  \item{y}{one or more columns of samples corresponding to \code{x}}
  \item{scale.trim}{how many times above and below the minimum and maximum scales to analyze?}
  \item{scale.exp}{scale output by \eqn{1 / scale ^ scale.exp}}
  \item{nscales}{number of scales to analyze}
  \item{min.scale}{minimum scale in units of \code{x}}
  \item{max.scale}{maximum scale in units of \code{x}}
  \item{scales}{a set of scales to analyze; overrides all other scale arguments}
  \item{loc}{the loci at which to evalues the wavelet function}
  \item{wave.fun}{a wavelet function}
}
\details{By default, a regular sequence of \code{x} values are generated and the wavelet function is evaluated at each of these. You can override this by specifying \code{loc}. Also a log-base-2 sequence of scales is generated spanning from the approximate (for irregular data) Nyquist frequency to unit frequency. The manner in which the sequence of scales is determined can be overriden in various ways as shown above. The user can subsitute the wavelet function, but it is assumed that the wavelet transfer function peaks at 2Pi radians.}
\value{an object of class "mvcwt" with fields \code{x}, \code{y} and \code{z}, which give the loci, scales and wavelet coefficient respectively.}
\references{Keitt, T. H. 2008. Coherent ecological dynamics induced by large-scale disturbance. Nature 454:331–4. doi: 10.1038/nature06935.}
\author{Timothy H. Keitt}
\seealso{\code{\link{Morelet}}, \code{\link{wmr}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, scale.trim = 1, scale.exp = 0.5, nscales = length(unlist(x)), 
    min.scale = 2 * scale.trim * median(diff(unlist(x))), max.scale = diff(range(unlist(x)))/scale.trim/2, 
    scales = log2Bins(min.scale, max.scale, nscales), loc = regularize(as.vector(unlist(x))), 
    wave.fun = "Morlet") 
{
    require(foreach)
    wave.fun = match.fun(wave.fun)
    x = as.vector(unlist(x))
    lmat = lagMat(x, loc)
    y = matrix(unlist(y), nrow = length(x))
    w = foreach(s = scales) \%dopar\% {
        Conj(wave.fun(lmat/s)/s^scale.exp) \%*\% y
    }
    w = array(unlist(w), dim = c(length(loc), ncol(y), length(scales)))
    w = aperm(w, perm = c(1, 3, 2))
    structure(list(x = loc, y = scales, z = w), class = "mvcwt")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{wavelets}