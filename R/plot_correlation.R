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
