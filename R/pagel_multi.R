#' Runs Pagel's correlation method for multiple trees.
#' 
#' Runs Pagel's test for binary character correlation for multiple trees. Borrows a few ideas and lines of code from the phytools function fitPagel.
#' @param tree Multiphylo object with phylogenetic trees.
#' @param a Character distribution for first character.
#' @param b Character distribution for second character.
#' @examples
#' pagel_multi(trees, a, b)
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
