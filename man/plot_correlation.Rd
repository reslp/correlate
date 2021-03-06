% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_correlation.R
\name{plot_correlation}
\alias{plot_correlation}
\title{Plots a conditional probability from correlate distribution as a violin plot.}
\usage{
plot_correlation(cond = cond, title = "", xlab = "", ylab = "")
}
\arguments{
\item{cond}{dataframe with conditional probabilities created with get_conditional.}

\item{title}{Title of plot. Default = no title}

\item{xlab}{Label for x-axis. Default = no label}

\item{ylab}{Label for y-axis. Default = no label}
}
\description{
This function plots a conditional probability distribution created with get_conditional as a violin plot. This functions depends on reshape2 and ggplot2.
}
\examples{
plot_correlation(cond, title, xlab, ylab)
}
