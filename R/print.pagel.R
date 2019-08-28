#' Prints the results of a multi-tree Pagel analysis.
#' 
#' This functions provides print functionality for a multi-tree Pagel analyses results created with pagel_multi.
#' @param pagels Results created by pagel_multi.
#' @examples
#' print(pagels)
print.pagel <- function(pagels) {
  cat("Summary of Pagel's correlation test\n")
  if (pagels$multi == T) {cat("Test performed on multiple trees.\n") }
  else {cat("Test performed on single tree.\n")}
  cat("\nHypothesis test results:\n")
  cat("\nLogLik(s) independent:")
  print(pagels$loglik_inds)
  cat("\nLogLik(s) dependent:")
  print(pagels$loglik_deps)
  cat("\nLogLik ratio(s):")
  print(pagels$ratios)
  cat("\nP-value(s):")
  print(pagels$p_values)
}
