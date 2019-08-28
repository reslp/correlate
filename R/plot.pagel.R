#' Plots the results of a multi-tree Pagel correlation analysis.
#' 
#' This function plots the results of a multi-tree Pagel correlation analysis created with pagel_multi. Depends on reshape2 and ggplot2.
#' @param pagel The results of pagel_multi.
#' @param what Select which parameter to plot. Options = "lik", "pvalue". Default = "pvalue"
#' @examples
#' plot(pagels, what="lik")
plot.pagel <- function(pagels, what="", ...) {
  library(reshape2)
  library(ggplot2)
  if (class(pagels) != "pagel") {
    cat("Object has to be of class \"pagel\" to plot.\n")
    return
  }
  if (what == "lik") {
    df <- data.frame(pagels$loglik_inds,pagels$loglik_deps)
    df <- melt(df)
    df$tree <- seq(1:length(pagels$loglik_inds))
    names(df) <- c("model","loglik","tree")
    p <- ggplot(df, aes(tree,loglik))
    p + geom_point(aes(colour=factor(model)), size=3)
  }
  else {
    l <- melt(pagels$p_values)
    ggplot(l) + geom_histogram(aes(x=value),breaks=seq(0, 1, by=0.02), ...) + scale_y_continuous(name="count")+ scale_x_continuous(name="p-value")+ geom_vline(xintercept = 0.05, color="red")
  }
}
