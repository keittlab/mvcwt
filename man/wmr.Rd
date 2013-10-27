\name{wmr}
\alias{wmr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wmr(w, smoothing = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{smoothing}{
%%     ~~Describe \code{smoothing} here~~
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
