#' Calculates the overlap of two stochastic character state distributions. Used internally by correlate.
#' 
#' @param c1_min Starting point (along branch) of character state 1.
#' @param c1_max End point (along branch) of character state 1.
#' @param c2_min Starting point (along branch) of character state 2.
#' @param c2_max End point (along branch) of character state 2.
#' @examples
#' calculate_overlap(0.1, 0.5, 0.005, 0.4)
calculate_overlap <- function(c1_min=0, c1_max=0, c2_min=0, c2_max=0){
  require(phytools)
  if (c1_max <= c2_min || c1_min >= c2_max) { return(0)} #case 1
  if ((c1_max - c2_max) <= 0 && (c1_min - c2_min) <= 0) { #case 3
    return (c1_max-c2_min)}
  if ((c1_max - c2_max) > 0 && (c1_min - c2_min) > 0) { #case 4
    return (c2_max-c1_min)} 
  if ((c1_min - c2_min) <= 0 && (c1_max - c2_max) >= 0) { #case 5
    return (c2_max-c2_min)}
  if ((c2_min - c1_min) < 0 && (c2_max - c1_max) > 0) { #case2
    return (c1_max-c1_min)}
  return(0) #case 6
}  