require(phytools)


#calculate_overlap <- function(c1_min=0, c1_max=0, c2_min=0, c2_max=0){
#  if (c1_max <= c2_min || c1_min >= c2_max) { return(0)}
#  if ((c1_max - c2_max) <= 0 && (c1_min - c2_min) <= 0) {
#    return (c1_max-c2_min)}
#  if ((c1_max - c2_max) > 0 && (c1_min - c2_min) > 0) {
#    return (c2_max-c1_min)}
#  if ((c1_min - c2_min) <=0 && (c1_max - c2_max) >= 0) {
#    return (c2_max-c2_min)}
#  if ((c2_min - c1_min) < 0 && (c2_max - c1_max) > 0) {
#    return (c1_max-c1_min)}
#  return(0)
#}  

calculate_overlap <- function(c1_min=0, c1_max=0, c2_min=0, c2_max=0){
  if (c1_max <= c2_min || c1_min >= c2_max) { return(0)} #case 1
  if ((c1_max - c2_max) <= 0 && (c1_min - c2_min) <= 0) { #case 3
    return (c1_max-c2_min)}
  if ((c1_max - c2_max) >= 0 && (c1_min - c2_min) >= 0) { #case 4
    return (c2_max-c1_min)} 
  if ((c1_min - c2_min) <= 0 && (c1_max - c2_max) >= 0) { #case 5
    return (c2_max-c2_min)}
  if ((c2_min - c1_min) <= 0 && (c2_max - c1_max) >= 0) { #
    return (c1_max-c1_min)}
  return(0) #case 6
}  

correlate <- function(ntrees=1,nmaps=1,smap1=simmap, smap2=simmap_multi, chars2= c("wood", "bark", "rock", "soil", "lichen"), chars1 = c(0,1)) {
### function currently only accepts lists of multiPhylo objects!!
if (class(smap1) != "list") { #single trees, behaviour unclear for otherwise formatted phylo objects
    sim <- vector("list", 1)
    sim[[1]] <- smap1 
} else { sim <- smap1 }

if (class(smap2) != "list") { #single trees, behaviour unclear for otherwise formatted phylo objects
    simm <- vector("list", 1)
    simm[[1]] <- smap2 
} else { simm <- smap2 }
  number_rows <- length(sim[[1]][[1]]$edge)
ncols <- length(chars1)*length(chars2)
df <- data.frame(matrix(0, ncol=ncols, nrow=ntrees))

whichmatrix <- 0
whichchar <- 0
for(chs2 in chars2) {
char2 <- chs2
for(chs1 in chars1) {
char1 <- chs1
whichchar <- whichchar + 1
colnames(df)[whichchar] <- paste(char1,"_",char2,sep="") 
for(tr in 1:ntrees) {
shared_time_tree <- vector()
for(mp in 1:nmaps) {
cat("For ", chs2, " and ", chs1, " Working on Tree ", tr, " Map ", mp, "\n")
for (i in 1:number_rows) { 
  char1_min <- 0
  char1_max <- sim[[tr]][[mp]]$edge.length[i]
  char2_min <- 0
  char2_max <- simm[[tr]][[mp]]$edge.length[i]  
  char_1_states <- unlist(sim[[tr]][[mp]]$maps[i])
  char_2_states <- unlist(simm[[tr]][[mp]]$maps[i])  
  if ((length(char_1_states) == 1) && (length(char_2_states) == 1) && names(char_1_states)==char1 && names(char_2_states)==char2) { 
    ## only one state mapped on edge
    #in this case correlation is 100% of edge length
    #print("Only one state")
    shared_time_tree <- c(shared_time_tree, calculate_overlap(char1_min, char1_max, char2_min, char2_max))
  } else if ((length(char_1_states) != 1) && (length(char_2_states) == 1) && names(char_1_states)==char1 && names(char_2_states)==char2) {
      ## two states of char 1 mapped on edge
      #print("Two states in char 1")
      char1_max <- 0
      for (j in 1:length(char_1_states)) {
        if (names(char_1_states[j]) == char1) {
          char1_max <- char1_max + as.numeric(char_1_states[j])
        } 
      }
      shared_time_tree <- c(shared_time_tree, calculate_overlap(char1_min, char1_max, char2_min, char2_max))
    } else if ((length(char_1_states) == 1) && (length(char_2_states) != 1) && names(char_1_states)==char1 && names(char_2_states)==char2) {
      #print("Two states in char 2")  
      ## two states of char 2 mapped on edge
        char2_max <- 0
        for (j in 1:length(char_2_states)) {
          if (names(char_2_states[j]) == char2) {
            char2_max <- char2_max + as.numeric(char_2_states[j])
            
          } 
        }
        shared_time_tree <- c(shared_time_tree, calculate_overlap(char1_min, char1_max, char2_min, char2_max))
    } else if ((length(char_1_states) > 1) && (length(char_2_states) > 1)) {
      #print("Both states >1")
        for (k in 1:length(char_1_states)) {
          if (names(char_1_states[k])==char1) {
            if (k == 1) {
              char1_min <- 0
              char1_max <- as.numeric(char_1_states[k])
            } else {
                char1_min <- as.numeric(sum(char_1_states[1:k-1]))
                char1_max <- as.numeric(sum(char_1_states[1:k]))
              }  
          for (l in 1:length(char_2_states)) {
            if(names(char_2_states[l])==char2) {
              if (l == 1) {
                char2_min <- 0
                char2_max <- as.numeric(char_2_states[l])
                shared_time_tree <- c(shared_time_tree, calculate_overlap(char1_min, char1_max, char2_min, char2_max))
                #cat("Total overlap:", simm[[tr]][[mp]]$mapped.edge[i], "\n")
                #cat("calculated overlap:", shared_time_tree, "\n")
              } else {
                      char2_min <- as.numeric(sum(char_2_states[1:l-1]))
                      char2_max <- as.numeric(sum(char_2_states[1:l])) 
                      shared_time_tree <- c(shared_time_tree, calculate_overlap(char1_min, char1_max, char2_min, char2_max))
                    } 
              
              #### here happens the important stuff          
            }
            }
          }
        }
      }   else { #shared probability on edge is 0
                shared_time_tree <- c(shared_time_tree, 0)
                }
} #loop i
} #loop mp

df[tr,whichchar] <- (sum(shared_time_tree)/nmaps)/sum(sim[[tr]][[1]]$edge.length) #divided by number of maps and tree length
} #loop tr
} #loop chs1
} #loop chs2
return(df)
} #function calc_overlap
# calculate relative times

#function get probability, will extract P(character) from the provided stochastic maps per tree
#needed for later calculation of conditional
get_prob <- function(sm){
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
    #print("char")
    #print(c)
    for (i in 1:ntrees) {
      #print("i")
      #print(i)
      edge_length <- sum(sm[[i]][[1]]$edge.length)
      #print(sum(edge_length))
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

## calculates conditional from dfs created with correlate() and get_prob()
get_conditional <- function(matrix, probs) {
  #number <- ncol(matrix) / ncol(probs)
  #probs_expanded <- do.call("cbind", replicate(number, probs, simplify=F))
  #conditional <- matrix / probs_expanded
  increment <- ncol(matrix) / ncol(probs)
  which <- 1
  for (i in 1:ncol(matrix)){
    print(i)
    if (iseven(i) == TRUE) {
      #print("even")
      #print(i)
      matrix[i] <- matrix[i]/probs[which]
      which <- which + 1
    }
    else { 
      #print("else")
      matrix[i] <- matrix[i]/probs[which]
      }
  }
  return(matrix)
}

plot_correlation <- function(cond=cond, title="") {
  require(reshape2)
  require(ggplot2)
  df <- melt(cond)
  p <- ggplot(data = df, aes(x=variable, y=value))  + geom_violin(aes(fill=variable)) + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black") + labs(title = "title") +theme(plot.title = element_text(size = rel(0.9)))
  p
  return(p)
}

#colnames(conditional)[1] <- names(matrix)[1]
#conditional[,1] <- matrix[1] / probs[1]

#colnames(conditional)[2] <- names(matrix)[2]
#conditional[,2] <- matrix[2] / probs[1]

#colnames(conditional)[3] <- names(matrix)[3]
#conditional[,3] <- matrix[3] / probs[2]

#colnames(conditional)[4] <- names(matrix)[4]
#conditional[,4] <- matrix[4] / probs[2]
