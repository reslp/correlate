setwd('/Users/sinnafoch/Dropbox/Philipp/correlate_package')
setwd('/Volumes/sinnafoch/sinnafoch/Dropbox/Philipp/correlate_package') #laptop
source("./functions_correlation.r")

####simulate character history
library(diversitree)
set.seed(2)
phy <- tree.bd(c(.1, .03), max.taxa=20) #lambda=0.1, mu=0.03

mbt <- max(branching.times(phy))
mbt

# rescale tree to edge of 1.0
phy$edge.length <- phy$edge.length / mbt
phy
states <- sim.character(phy, c(12, 12), x0=0, model="mk2") #0->1 = 0.1, 1->0 = 0.2
states
states2 <- sim.character(phy, c(.5, .5), x0=0, model="mk2") #0->1 = 0.01, 1->0 = 0.02
states2
states2[states2==1] <- "A"
states2[states2==0] <- "B"
states2


foo <- vector("list", 3)
class(foo) <- "multiPhylo"

foo[[1]] <- phy
foo[[2]] <- phy
foo[[3]] <- phy


phy
?make.simmap()
#####perform stochastic mapping
library(phytools)
sm1 <- list()
for (i in 1:length(foo)) {
  sm1[[i]] <- make.simmap(foo[[i]], states, model="SYM", nsim=100)
}
sm2 <- list()
for (i in 1:length(foo)) {
  sm2[[i]] <- make.simmap(foo[[i]], states2, model="SYM", nsim=100)
}
sm2


#### get shared state information from two stochastic maps
### char codings are: wood, soil, bark, rock, lichenicolous
source("./functions_correlation.r")
matrix <- correlate(ntrees=100, nmaps=100, chars2=c("bark","lichenicolous", "rock", "soil", "wood"), chars1=c(0,1),smap1=simmap_bi,smap2=simmap_multi)
out <- get_prob(simmap_multi) #get probabilities of each state = P(B)
out <- out[1:10,] #subsample for first 10 trees
cond <- get_conditional(matrix=matrix, probs=out)
cond

plot_correlation(cond)






