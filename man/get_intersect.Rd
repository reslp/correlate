\name{get_intersect}
\alias{get_intersect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get_intersect
%%  ~~function to do ... ~~
}
\description{
This is the core function of the correlate package. It will create the intersect probability  for pairwise character combinations from stochastic maps.
}
\usage{
get_intersect(ntrees = 1, nmaps = 1, smap1 = simmap, smap2 = simmap_multi, chars2 = c("A","B"), chars1 = c(0, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ntrees}{
  Specifies how many trees are analyzed.
}
  \item{nmaps}{
Specifies how many maps per tree should be analyzed.
}
  \item{smap1}{
stochastic maps object for first set of characters.
}
  \item{smap2}{
stochastic maps object for second set of characters.
}
  \item{chars1}{
Specifies which characters from first character set should be analyzed.
}
  \item{chars2}{
Specifies which characters from second character set should be analyzed.
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
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
