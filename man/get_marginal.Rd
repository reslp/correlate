\name{get_marginal}
\alias{get_marginal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the marginal probability of character occurrence.
}
\description{
This function calculates the marginal probability of observing a particular character on the provided phylogenetic tree(s). Will return a dataframe containing the marginal probabilities of character occurrence for each character.
}
\usage{
get_marginal(smap)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smap}{
List of multiPhylo objects containing stochastic character maps created with phytools.
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
Philipp Resl
}
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
\keyword{marginal}% use one of  RShowDoc("KEYWORDS")
\keyword{probability}% __ONLY ONE__ keyword per line
