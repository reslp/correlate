#' @export
#' 
require(phytools)

calculate_overlap <- function(c1_min=0, c1_max=0, c2_min=0, c2_max=0){
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