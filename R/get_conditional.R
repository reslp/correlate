#' Calculate conditional probability of character co-occurrence
#' 
#' This function calculates the conditional probability distribution of observing two character states along the provided phylogenetic tree(s).
#'
#' @param matrix matrix with intersect probabilities obtained from the function get_intersect.
#' @param probs dataframe with marginal probabilities obtained from the function get_marginal.
#' @examples
#' get_conditional(matrix, probs)
get_conditional <- function(matrix, probs) {
  increment <- ncol(matrix) / ncol(probs)
  which <- 1
  for (i in 1:ncol(matrix)){
    print(i)
    if (iseven(i) == TRUE) {
      matrix[i] <- matrix[i]/probs[which]
      which <- which + 1
    }
    else {
      matrix[i] <- matrix[i]/probs[which]
    }
  }
  return(matrix)
}
