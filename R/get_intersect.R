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
