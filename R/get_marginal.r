#' @export
#function get_marginal, will extract P(character) from the provided stochastic maps per tree
#needed for later calculation of conditional
get_marginal <- function(sm){
  if (class(sm) != "list") { #single trees, behaviour unclear for otherwise formatted phylo objects
    sim <- vector("list", 1)
    sim[[1]] <- sm
    sm <- sim
  }
  ntrees <- length(sm)
  nmaps <- length(sm[[1]])
  time_matrix <- vector()
  edge_length <- vector()
  time_summary <- vector()
  characters <- colnames(sm[[1]][[1]]$mapped.edge)
  df <- data.frame(matrix(0,nrow=ntrees,ncol=length(characters)))
  whichchar <- 0
  for (c in characters) {
    whichchar <- whichchar + 1
    colnames(df)[whichchar] <- c
    for (i in 1:ntrees) {
      edge_length <- sum(sm[[i]][[1]]$edge.length)
      time_matrix <- vector()
      for (j in 1:nmaps) {
        time_matrix  <- c(time_matrix,sum(as.numeric(sm[[i]][[j]]$mapped.edge[,whichchar]))) #1=character 0, 2=character 1,...
      }
      time_summary <- c(time_summary,(sum(time_matrix,na.rm=T)/nmaps)) ## maybe na.rm=T NA needed in sum
      df[i,whichchar] <- sum(time_summary/edge_length)
      time_matrix <- vector()
      time_summary <- vector()
    }

  }
  return(df)
}

iseven <- function(x) x %% 2 == 0
