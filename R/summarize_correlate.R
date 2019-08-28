#' Creates a summary overview of a correlate analysis.
#' 
#' This functions provides summatized output of a correlate output in tabular format. It will apply pairwise wilcox tests to compare the means of conditional probability distributions for character co-occurence.
#' @param cond Results created by get_conditional.
#' @examples
#' summarize_correlate(cond)
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
