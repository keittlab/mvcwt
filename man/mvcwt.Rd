\name{mvcwt}
\alias{mvcwt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvcwt(x, y, scale.trim = 1, scale.exp = 0.5, nscales = length(unlist(x)), min.scale = 2 * scale.trim * median(diff(unlist(x))), max.scale = diff(range(unlist(x)))/scale.trim/2, scales = log2Bins(min.scale, max.scale, nscales), loc = regularize(as.vector(unlist(x))), wave.fun = "Morlet")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{scale.trim}{
%%     ~~Describe \code{scale.trim} here~~
}
  \item{scale.exp}{
%%     ~~Describe \code{scale.exp} here~~
}
  \item{nscales}{
%%     ~~Describe \code{nscales} here~~
}
  \item{min.scale}{
%%     ~~Describe \code{min.scale} here~~
}
  \item{max.scale}{
%%     ~~Describe \code{max.scale} here~~
}
  \item{scales}{
%%     ~~Describe \code{scales} here~~
}
  \item{loc}{
%%     ~~Describe \code{loc} here~~
}
  \item{wave.fun}{
%%     ~~Describe \code{wave.fun} here~~
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
