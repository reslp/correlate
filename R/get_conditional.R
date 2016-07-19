## calculates conditional from dfs created with correlate() and get_prob()
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
