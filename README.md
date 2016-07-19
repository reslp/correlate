correlate
=========

correlate performes character correlation analyses of categorical (discrete) characters based on stochastic character mappings.

Change History
=================
July 19, 2016

	 - new version 0.2
	 - added function to perform Pagel'S 1994 correlation test on multiple trees
	 - functions for correlation test have been renamed for consistency
	 - small bug fixes
	


April 17, 2016:

	- added documentation for some functions


April 15, 2016:

	- added correct behavior for single trees
	- added function summarize_correlate() which produces basic summary statistics.
	- added possibility to label axis in plot_correlation()
	- some small code clean-ups


DESCRIPTION
===========

This small package performs character correlation analyses based on stochastic character maps from two sets of categorical characters. This also works for sets of trees (e.g. from a Bayesian posterior tree sample), to account for phylogenetic uncertainty. The functions in correlate will calculate the conditional probability for character pairs occurring on a particular point on the phylogenetic tree and output them as summarized probability distributions in the analyzed tree space. This method refers to concepts described in [this paper](http://sysbio.oxfordjournals.org/content/52/2/131.short) by Huelsenbeck et al. (2003) 


REQUIREMENTS
============

- R 3.2.2
- R dependecies: devtools, phytools,reshape, ggplot2 (if you want to use the included plotting function)
- one or more phylogenetic trees
- characters to map


INSTALL
========

First you need to install and load the devtools R package:

`install.packages("devtools") `

`library(devtools) `


Now you can install correlate:

`install_github("reslp/correlate")`





USAGE
=====
I am working on a detailed documentation of the functions in correlate. Meanwhile here is a very brief describtions of what steps are needed for correlate to run:

1. Lets assume you have a set of 10 stochastic maps for each of 10 trees created with phytools for two character sets named `map1`and `map2` 
2. First load correlate: `library(correlate)`
3. Now you can call the function correlate:
`corr_matrix <- get_intersect(ntrees=10, nmaps=10, chars2=c("bark","lichenicolous", "rock", "soil", "wood"), chars1=c(0,1),smap1=map1,smap2=map2)` 
ntrees and nmaps specify the number of trees and stochastic maps, chars1 and chars2 which characters should be analyzed in a pairwise matter. smap1 and smap2 specify the stoachstic maps.
4. To calculate the conditional you will also need the probability for each character to appear alone on the tree(s): `prob <- get_marginal(simmap_multi)`. This is only needed for one of the mappings, depending on the question you want to ask. See the [wiki](https://en.wikipedia.org/wiki/Conditional_probability) on conditional probaility for further information.
5. Now you can calculate the conditional probability: `cond <- get_conditional(matrix=corr_matrix, probs=prob)`
6. correlate also provides the output of the obtained conditional probability distributions as Violin plots: `plot_correlation(cond, title="title")`
(https://github.com/reslp/correlate/blob/master/correlate_example.png). This shows the conditional for binary characters 0 and 1 with the multistate character. For example the very left violin plot 0_bark could be read as: The probability of bark given that character 0 is observed.
7. By running the built-in function `summarize_correlate(cond)` you will get basic information about the distributions including, means, variance and standard deviation. It will also compare distributions with pairwise wilcox tests.


Pagel'S (1994) test:

This test can be performed with the function `pagel_multi`. It takes a list of trees and two character distributions for the tree tips as named vectors as input: `pagel_multi(trees, character_a, character_b)`. It will fit two models to each tree in the list and perform chi-squared tests on log-likelihood differences of the models. The function returns an object of class "pagel", which includes log-likelihoods for both models, likelihood ratios and p-values of the chi-squared test. This object may be plotted with plot(pagel). Plotting requires ggplot2.



LIMITATIONS
===========
Although functional, correlate is still in a very early development stage, so be careful when using this experimental software. I give no warranty! Also, there are several features I would like to include. Among those are:

- Correlation of multiple (more then two) sets of characters
- Improved plotting
- higher speed and memory efficiency



COPYRIGTH AND LICENSE
=====================

Copyright (C) 2015 Philipp Resl

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program in the file LICENSE. If not, see http://www.gnu.org/licenses/.
