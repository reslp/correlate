% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.pagel.R
\name{plot.pagel}
\alias{plot.pagel}
\title{Plots the results of a multi-tree Pagel correlation analysis.}
\usage{
\method{plot}{pagel}(pagels, what = "", ...)
}
\arguments{
\item{what}{Select which parameter to plot. Options = "lik", "pvalue". Default = "pvalue"}

\item{pagel}{The results of pagel_multi.}
}
\description{
This function plots the results of a multi-tree Pagel correlation analysis created with pagel_multi. Depends on reshape2 and ggplot2.
}
\examples{
plot(pagels, what="lik")
}
