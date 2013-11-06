\name{mvcwt}
\alias{mvcwt}
\alias{print.mvcwt}
\title{Computes the wavelet transform of a multivariate time series}
\description{This function takes set a set of seqences as columns in a matrix and computes the continuous wavelet transform on each.}
\usage{
mvcwt(x, y,
      scale.trim = 1,
      scale.exp = 0.5,
      nscales = length(as.vector(unlist(x))),
      min.scale = 2 * scale.trim * median(diff(as.vector(unlist(x)))),
      max.scale = diff(range(as.vector(unlist(x)))) / scale.trim / 2,
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
\details{By default, a regular sequence of \code{x} values are generated and the wavelet function is evaluated at each of these. You can override this by specifying \code{loc}. A log-base-2 sequence of scales is generated spanning from the approximate (for irregular data) Nyquist frequency to unit frequency. The manner in which the sequence of scales is determined can be overriden in various ways as shown above. The user can subsitute the wavelet function, but it is assumed that the wavelet transfer function peaks at 2Pi radians, i.e., \code{wave.fun(lag / s)} is centered at scale \code{s}. The output energy or wavelet variance is scaled according to \code{1 / scale ^ scale.exp}. The default \code{scale.exp = 0.5} gives a flat white-noise scalogram, but the choice is arbitrary and you may wish to emphsize or demphasize coarse scale variation.}
\value{an object of class "mvcwt" with fields \code{x}, \code{y} and \code{z}, which give the loci, scales and wavelet coefficient respectively.}
\note{This function uses \code{dopar} from the \code{foreach} package. Registering a parallel backend is highly recommneded for large datasets. Note that \code{x} and \code{y} are coersed to vector and matrix respectively. If you feed this function something that cannot be coersed correctly, you will get strange results.}
\references{Keitt, T. H. 2008. Coherent ecological dynamics induced by large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.}
\author{Timothy H. Keitt}
\seealso{\code{\link{Morlet}}, \code{\link{wmr}}}
\examples{
run.it = getOption("run.long.examples")     
if ( !is.null(run.it) && run.it )
{
  data(sunspot.year)
  w = mvcwt(time(sunspot.year), sunspot.year, min.scale = 5, max.scale = 40)
  image(w, "Mod", reset.par = FALSE)
  contour(w, "Mod", bound = NA, drawlabels = FALSE, add = TRUE)
}}
\keyword{wavelets}