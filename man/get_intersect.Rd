\name{get_intersect}
\alias{get_intersect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating character overlap
}
\description{
This function calculates the overlapping parts of two (or more) sets of stochastically mapped categorical characters. As input it takes two sets of stochastic maps created with the R package phytools. It will return the proportion of the tree where two characters overlap.
}
\usage{
get_intersect(ntrees = 1, nmaps = 1, smap1 = simmap, smap2 = simmap_multi, chars2 = c(0, 1), chars1 = c(0, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ntrees}{
 The number of trees which should be analysed. default = 1
}
  \item{nmaps}{
The number of stochastic character maps which should be analysed per tree. default = 1
}
  \item{smap1}{
List of multiPhylo Objects containing stochastic maps for first character. Created with phytools.
}
  \item{smap2}{
List of multiPhylo Objects containing stochastic maps for second character. Created with phytools.
}
  \item{chars2}{
Vector containing character states which should be analysed for the first set of characters.
}
  \item{chars1}{
Vector containing character states which should be analysed for the second set of characters.
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

}
\author{
Philipp Resl
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
add example code

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{correlate}% use one of  RShowDoc("KEYWORDS")
\keyword{intersect}% __ONLY ONE__ keyword per line
