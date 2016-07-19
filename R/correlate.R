## correlate
## written by Philipp Resl
## last update May 28,  2016

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

get_intersect <- function(ntrees=1,nmaps=1,smap1=simmap, smap2=simmap_multi, chars2= c(0,1), chars1 = c(0,1)) {
### function currently only accepts lists of multiPhylo objects, stochastic maps of single trees will be converted to lists
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
pb <- txtProgressBar(min=1, max=length(chars2)*length(chars1)*ntrees*nmaps, style=3)
progress <- 0
cat("Calculating intersect Probability: \n")
cat("character set1: ", chars1, "\n")
cat("character set2: ", chars2, "\n")
cat("Number of trees: ", ntrees, "\n")
cat("Number of maps per tree: ", nmaps, "\n")
for(chs2 in chars2) 
  {
    char2 <- chs2
    for(chs1 in chars1) 
      {
        char1 <- chs1
        whichchar <- whichchar + 1
        colnames(df)[whichchar] <- paste(char1,"&",char2,sep="") 
        for(tr in 1:ntrees) 
          {
            
            shared_time_tree <- vector()
            for(mp in 1:nmaps) 
              {
                progress <- progress+1
                setTxtProgressBar(pb, progress)
                for (i in 1:number_rows) 
                  { 
                    char1_min <- 0
                    char1_max <- sim[[tr]][[mp]]$edge.length[i]
                    char2_min <- 0
                    char2_max <- simm[[tr]][[mp]]$edge.length[i]  
                    char_1_states <- unlist(sim[[tr]][[mp]]$maps[i])
                    char_2_states <- unlist(simm[[tr]][[mp]]$maps[i])  
                    if ((length(char_1_states) == 1) && (length(char_2_states) == 1) && names(char_1_states)==char1 && names(char_2_states)==char2) 
                      { 
                        # only one state mapped on edge
                        # in this case correlation is 100% of edge length
                        shared_time_tree <- c(shared_time_tree, calculate_overlap(char1_min, char1_max, char2_min, char2_max))
                      } 
                    else if ((length(char_1_states) != 1) && (length(char_2_states) == 1) && names(char_1_states)==char1 && names(char_2_states)==char2) 
                      {
                        # two states of char 1 mapped on edge
                        char1_max <- 0
                        for (j in 1:length(char_1_states)) 
                          {
                            if (names(char_1_states[j]) == char1) 
                              {
                                char1_max <- char1_max + as.numeric(char_1_states[j])
                              } 
                          }
                        shared_time_tree <- c(shared_time_tree, calculate_overlap(char1_min, char1_max, char2_min, char2_max))
                      }
                    else if ((length(char_1_states) == 1) && (length(char_2_states) != 1) && names(char_1_states)==char1 && names(char_2_states)==char2) 
                      {
                        #print("Two states in char 2")  
                        ## two states of char 2 mapped on edge
                        char2_max <- 0
                        for (j in 1:length(char_2_states)) 
                          {
                            if (names(char_2_states[j]) == char2) 
                              {
                                char2_max <- char2_max + as.numeric(char_2_states[j])
                              } 
                          }
                        shared_time_tree <- c(shared_time_tree, calculate_overlap(char1_min, char1_max, char2_min, char2_max))
                      } 
                    else if ((length(char_1_states) > 1) && (length(char_2_states) > 1)) 
                      {
                        #print("Both states >1")
                        for (k in 1:length(char_1_states)) 
                          {
                            if (names(char_1_states[k])==char1) 
                              {
                                if (k == 1) 
                                  {
                                    char1_min <- 0
                                    char1_max <- as.numeric(char_1_states[k])
                                  } 
                                else 
                                  {
                                    char1_min <- as.numeric(sum(char_1_states[1:k-1]))
                                    char1_max <- as.numeric(sum(char_1_states[1:k]))
                                  }  
                                for (l in 1:length(char_2_states)) 
                                  {
                                    if(names(char_2_states[l])==char2) 
                                      {
                                        if (l == 1) 
                                          {
                                            char2_min <- 0
                                            char2_max <- as.numeric(char_2_states[l])
                                            shared_time_tree <- c(shared_time_tree, calculate_overlap(char1_min, char1_max, char2_min, char2_max))
                                          } 
                                        else 
                                          {
                                            char2_min <- as.numeric(sum(char_2_states[1:l-1]))
                                            char2_max <- as.numeric(sum(char_2_states[1:l])) 
                                            shared_time_tree <- c(shared_time_tree, calculate_overlap(char1_min, char1_max, char2_min, char2_max))
                                          } 
                                      }
                                  }
                              }
                          }
                      }   
                      else 
                        { #shared probability on edge is 0
                          shared_time_tree <- c(shared_time_tree, 0)
                        }
                  } #loop i
              } #loop mp
              df[tr,whichchar] <- (sum(shared_time_tree)/nmaps)/sum(sim[[tr]][[1]]$edge.length) #divided by number of maps and tree length
          } #loop tr
      } #loop chs1
  } #loop chs2
  #my_intersect <- list(df)
  #class(my_intersect) <- "correlate_intersect"
  return(df)
} #function calc_overlap
# calculate relative times

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

## fit Pagels test for binary character correlation for multiple trees
# borrows a few ideas and lines of code from the phytools function fitPagel

pagel_multi <- function(trees, a, b) {
  if (!is.factor(a)) a <- as.factor(a)
  names_a <- levels(a)
  if (!is.factor(b)) b <- as.factor(b)
  names_b <- levels(b)
  comb_names <- setNames(factor(paste(a,b,sep="_"), levels=sapply(names_a,paste,names_b,sep="_")),names(a))
  
  #prepare transition rate matrices
  Q_ind <- matrix(c(0, 1, 2, 0, 3, 0, 0, 2, 4, 0, 0, 1, 0, 4, 3, 0), 4, 4, byrow = TRUE)
  rownames(Q_ind) <- levels(comb_names)
  colnames(Q_ind) <- levels(comb_names)
  Q_dep <- matrix(c(0, 1, 2, 0, 3, 0, 0, 4, 5, 0, 0, 6, 0, 7, 8, 0), 4, 4, byrow = TRUE)
  rownames(Q_dep) <- levels(comb_names)
  colnames(Q_dep) <- levels(comb_names)
  ##fit models to trees:
  lik.ratios <- vector() 
  p.values <- vector()
  lik.inds <- vector()
  lik.deps <- vector()
  cat("Pagel's correlation test:\n")
  cat("Analyzing ", length(trees), "trees\n")
  pb <- txtProgressBar(min=1, max=length(trees), style=3)
  for (i in 1:length(trees)) {
    setTxtProgressBar(pb, i)
    ###set up independent model
    ace_ind <- ape::ace(comb_names,trees[[i]],type="discrete",model=Q_ind)
    print(ace_ind$rates)
    ###dependent model
    ace_dep <- ape::ace(comb_names,trees[[i]],type="discrete",model=Q_dep)
    print(ace_dep$rates)
    lik.ind <- logLik(ace_ind)
    lik.dep <- logLik(ace_dep)
    
    lik.inds <- c(lik.inds, lik.ind)
    lik.deps <- c(lik.deps, lik.dep)
    
    lik.ratios <- c(lik.ratios, 2*(lik.dep-lik.ind))
    p.values <- c(p.values, pchisq(2*(lik.dep-lik.ind), df=length(levels(a))+length(levels(b)),lower.tail=F))
  }
  my_pagel <- list(loglik_inds = lik.inds, loglik_deps = lik.deps, ratios = lik.ratios, p_values=p.values, multi=TRUE)
  class(my_pagel) <- "pagel"
  return(my_pagel)
}

## fit Pagels test for binary character correlation
# borrows a few ideas and lines of code from the phytools function fitPagel

pagel <- function(tree, a, b) {
  if (!is.factor(a)) a <- as.factor(a)
  names_a <- levels(a)
  if (!is.factor(b)) b <- as.factor(b)
  names_b <- levels(b)
  comb_names <- setNames(factor(paste(a,b,sep="_"), levels=sapply(names_a,paste,names_b,sep="_")),names(a))
  
  ###independent model
  Q_ind <- matrix(c(0, 1, 2, 0, 3, 0, 0, 2, 4, 0, 0, 1, 0, 4, 3, 0), 4, 4, byrow = TRUE)
  rownames(Q_ind) <- levels(comb_names)
  colnames(Q_ind) <- levels(comb_names)
  #fit ace
  ace_ind <- ape::ace(comb_names,tree,type="discrete",model=Q_ind)
  #get transition rates
  rates_ind <- ace_ind$index.matrix
  rates_ind[rates_ind==0] <-NA
  for (i in 1:length(ace_ind$rates)) Q_ind <- ifelse(rates_ind==i,ace_ind$rates[i],Q_ind)
  Q_ind[is.na(Q_ind)] <- 0
  #calculate likelihood
  lik.ind <- logLik(ace_ind)
  
  ###dependent model
  Q_dep <- matrix(c(0, 1, 2, 0, 3, 0, 0, 4, 5, 0, 0, 6, 0, 7, 8, 0), 4, 4, byrow = TRUE)
  rownames(Q_dep) <- levels(comb_names)
  colnames(Q_dep) <- levels(comb_names)
  #fit
  ace_dep <- ape::ace(comb_names,tree,type="discrete",model=Q_dep)
  #get transition rates
  rates_dep <- ace_dep$index.matrix
  rates_dep[rates_dep==0] <-NA
  for (i in 1:length(ace_dep$rates)) Q_dep <- ifelse(rates_dep==i,ace_dep$rates[i],Q_dep)
  Q_dep[is.na(Q_dep)] <- 0
  #calculate likelihood
  lik.dep <- logLik(ace_dep)
  
  ##likelihood ratio test
  lik.ratio = 2*(lik.dep-lik.ind)
  
  ##p-value
  p = pchisq(lik.ratio, df=length(levels(a))+length(levels(b)),lower.tail=F)
  
  my_pagel <- list(Q_ind = Q_ind, Q_dep = Q_dep, loglik_ind = lik.ind, loglik_dep = lik.dep, ratio = lik.ratio, p_value=p, multi=FALSE)
  class(my_pagel) <- "pagel"
  return(my_pagel)
}

plot.pagel <- function(pagels, what="",...) {
  library(reshape2)
  library(ggplot2)
  if (class(pagels) != "pagel") {
    cat("Object has to be of class \"pagel\" to plot.\n")
    return
  }
  if (what == "lik") {
    df <- data.frame(pagels$loglik_inds,pagels$loglik_deps)
    df <- melt(df)
    df$tree <- seq(1:length(pagels$loglik_inds))
    names(df) <- c("model","loglik","tree")
    p <- ggplot(df, aes(tree,loglik))
    p + geom_point(aes(colour=factor(model)), size=3)
  }
  else {
    l <- melt(pagels$p_values)
    ggplot(l) + geom_histogram(aes(x=value),breaks=seq(0, 1, by=0.02), ...) + scale_y_continuous(name="count")+ scale_x_continuous(name="p-value")+ geom_vline(xintercept = 0.05, color="red")
  }
}

print.pagel <- function(pagels) {
  cat("Summary of Pagels correlation test\n")
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


