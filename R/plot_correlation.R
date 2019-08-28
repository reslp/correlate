#' Plots a conditional probability from correlate distribution as a violin plot.
#' 
#' This function plots a conditional probability distribution created with get_conditional as a violin plot. This functions depends on reshape2 and ggplot2.
#' @param cond dataframe with conditional probabilities created with get_conditional.
#' @param title Title of plot. Default = no title
#' @param xlab Label for x-axis. Default = no label
#' @param ylab Label for y-axis. Default = no label
#' @examples
#' plot_correlation(cond, title, xlab, ylab)
plot_correlation <- function(cond=cond, title="", xlab="", ylab="") {
  require(reshape2)
  require(ggplot2)
  df <- melt(cond)
  p <- ggplot(data = df, aes(x=factor(variable), y=value,fill=variable))  + geom_violin() + stat_summary(fun.y="mean",geom="point", color="black")
  if (xlab != "") {p <- p + xlab(xlab)}
  if (ylab != "") {p <- p + ylab(ylab)}
  p <- p + labs(title = title) +theme(plot.title = element_text(size = rel(0.9)))
  return(p)
}
