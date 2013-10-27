\name{wmr.boot}
\alias{wmr.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wmr.boot(w, smoothing = 1, reps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{smoothing}{
%%     ~~Describe \code{smoothing} here~~
}
  \item{reps}{
%%     ~~Describe \code{reps} here~~
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
