\name{get_conditional}
\alias{get_conditional}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate conditional probability of character co-occurrence
}
\description{
This function calculates the conditional probability distribution of observing two character states along the provided phylogeneti tree(s).
}
\usage{
get_conditional(matrix, probs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{
matrix with intersect probabilities obtained from the function get_intersect.
}
  \item{probs}{
dataframe with marginal probabilities obtained from the function get_marginal.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{conditional}% use one of  RShowDoc("KEYWORDS")
\keyword{probability}% __ONLY ONE__ keyword per line
\keyword{distribution}% __ONLY ONE__ keyword per line
\keyword{Violin}% __ONLY ONE__ keyword per line
