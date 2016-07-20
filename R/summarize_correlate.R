#' @export
summarize_correlate <- function (cond = cond) {
  require(reshape2)
  cat ("Summary statistics of character correlation based on ", nrow(cond), " trees.\n")
  cat("condition\tmean\tvariance\tstddev\n")
  list <- vector()
  for (row in 1:ncol(cond)) {
    names <- unlist(strsplit(colnames(cond)[row], "&"))
    names <- paste("P(",names[1], "|", names[2], ")", sep="")
    list <- c(list, names)
    cat(names, ":\t", mean(cond[,row]), "\t",var(cond[,row]),"\t",sd(cond[,row]), "\n")
  }
  colnames(cond) <- list
  d <- melt(cond)
  pairwise.wilcox.test(d$value, d$variable, p.adjust = "none")
}
