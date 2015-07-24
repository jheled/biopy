Version 0.1.6 contains the following utilities which may be of interest to non programmers.

## BEAST utilities ##

**starbeast\_demog\_log.py**

> A specialized utility for `*`BEAST users. The population size numbers as reported in the BEAST log can be misleading, since they may refer to diffrent populations at diffrent stages of the MCMC chain. This utility generates a new trace file, viewable in tracer, tracing population size(s) for each population encountered during the run.

**summary\_tree.py**

> Create a summary tree from posterior trees. BEAST tree annotator utility first picks up a topology, then independently assigns heights to each internal node. Sometimes this can be suboptimal - in the worst case the tree contains negative branch lengths.  summary\_tree takes a more global approach and looks for a tree which minimizes the overall distance to the whole set of posterior trees, using a tree distance measure.

**starbeast\_posterior\_popsizes.py**

> Annotate a `*`BEAST summary tree (such as the one from summary\_tree above) with posterior estimates of population sizes.

**sptree\_plot.py**

> Generate a figure showing all posterior species trees from a `*`BEAST analysis. The figure can be used to visually examine the uncertainty of both divergence times and population sizes.

**multispecies\_coalescent\_estimate.py**

> Heuristically estimate the birth rate and effective population size from a multispecies data in a `*`BEAST XML file.

## Sampling trees ##

**sample\_ranked\_tree.py**

> Uniformly sample a ranked labelled tree from the space of all ranked labelled trees with a fixed **unranked** topology.

**sample\_coalescent\_tree.py**

> Draw a gene tree using the coalescent and a given an effective population size function. Constant, stepwise constant and linear piecewise functions are supported.

**sample\_species\_tree.py**

> Draw a species tree and associated population sizes using a birth/death process and simple strategies for drawing population sizes.

**genetree\_in\_sptree\_sim.py**

> Generate gene trees compatible with species tree using the multispecies coalescent.

## Likelihood calculations ##

**genetree\_in\_sptree\_probs.py**

> Compute probabilities of all possible gene trees evolving according to the multispecies coalescent inside species tree. Due to combinatorial explosion, this is good only for **very** small cases.


## Sequence generation ##

**simulate\_sequences.py**

> Simulation sequences for tips of a gene tree using a standard GTR substitution model. Juker-Cantor and HKY are provided for convenience.