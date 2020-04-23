# clust-perturb

`clust.perturb` is an R package designed to test the robustness of graph-based clusters by randomly rewiring a portion of the underlying network and reclustering. Robust clusters “stick together” even after network rewiring, while non-robust clusters do not. Robustness is quantified using maximum Jaccard index.

## System requirements

`clust.perturb` relies on functions from the following R packages:

```
	igraph (>= 1.2.4.2)
	readr (>= 1.3.1)
	Matrix (>= 1.2-17)
	cluster (>= 2.1.0)
```

## Installation

To install `clust.perturb`, first install the devtools package, if it is not already installed: 

```r
install.packages("devtools") 
```

Then, install `clust.perturb` from GitHub: 

```r
devtools::install_github("GregStacey/clust-perturb")
```

This should take no more than a few minutes.


## Usage

The main function `clust.perturb` takes a network as input. The network is formatted as a two-column  dataframe, which specifices the edge list of the network. Therefore, you must have a network in the form of an edge list in order to use `clust.perturb`. Two input parameters control the functioning of `clust.perturb`: `noise` is a scalar between 0 and 1 that specifies the fraction of the network that will be rewired on each iteration, and `iters` is a positive integer that specifies the number of iterations. 

`noise`: In order to resolve robust and non-robust clusters, `noise` should be set to a "medium amount", such that some clusters remain intact while others are disrupted. This could be accomplished by choosing `noise` such that the variance of the `repJ` value returned by `clust.perturb` is close to maximum, although in practice the range `0.1<noise<0.2` is often sufficient.

`iters: `Since `clust.perturb` clusters the network `iters+1` times, `clust.perturb` will be time-intensive if the underlying clustering algorithm is time-intensive and/or the network is large. Therefore, users can benefit from choosing a value for `iters` that is minimally sufficient to estimate robustness. In practice, `iters=5` is often sufficient to estimate `repJ` within +/-0.05 (95% CI).

`edge.list.format` and `cluster.format`: `clust.perturb` perturbs the network by shuffling a portion of the inputted edge list. Therefore, clustering algorithms that take input other than an edge list will require a conversion function, specified by the argument `edge.list.format`, to convert the edge list into the input format required by the clustering function. `edge.list.format` must take exactly one argument: the 2-column edge list datafram. Similarly, `cluster.format` formats the output of the clustering function into a common format, namely a character vector of semicolon-separated nodes. `cluster.format` must take exactly two arguments: the first is the output object of the clustering algorithm, and the second is the sorted character vector of unique nodes from the original network, i.e. `unique(c(network[,1], network[,2]))`. For clustering algorithms such as `mcl` that return node indices but not names, this second argument is used to include node names in clusters; otherwise this second argument can be ignored in `cluster.format`, although `cluster.format` must take two arguments even if the second is not used.

However, four clustering algorithms are included in `clust.perturb`, namely `k-med`, `MCL`, `walktrap`, and `hierarchical`. When using these clustering algorithms, it is not necessary to pass conversion functions. Simply run `clust.perturb` on the network. We provide a test network `corum_5000.csv`, which is 5000 edges selected from the binarized [CORUM 3.0 network](https://mips.helmholtz-muenchen.de/corum/#download).

```r
network = as.data.frame(read_csv("corum_5000.csv"))
clusts = clust.perturb(network, clustering.algorithm = "hierarchical")
```

To confirm that the default `noise=0.1` is appropriate, visualize the `repJ` values and confirm they roughly span the range 0 to 1. If `repJ` values are closer to 1, a higher `noise` might be required. Conversely if `repJ` is low, try a smaller `noise` parameter.

```r
clusts1 = clust.perturb(network, clustering.algorithm = "hierarchical", noise = 0.01) # lower noise
clusts2 = clust.perturb(network, clustering.algorithm = "hierarchical", noise = 0.2) # higher noise
```

`cluster.perturb` is a general purpose wrapper for many clustering algorithms. To use an arbitrary algorithm so, the arguments `clustering.algorithm`, `edge.list.format`, and `cluster.format` must be functions. Here is an example that explicitly sets these functions for hierarchical clustering

```r
library(igraph)
library(RFLPtools)

# hierarchical clustering algorithm
alg = function(x) cutree(hclust(x, method="average"), k = 50)

# function for converting edge-list to dist object
# hclust requires a dist object
ef = function(x) {
  # make dense, triu adjacency matrix
  am = as_adjacency_matrix(graph_from_edgelist(as.matrix(x), directed = F), sparse = F)
  am[am>1] = 1
  # convert to dist object
  d = sim2dist(am, maxSim = 1)
  return(d)
}

# function for converting hierarchical output to semicolon-separated character vector
cf = function(x, y) {
  tmp = character()
  unqclusts = unique(x)
  for (ii in 1:length(unqclusts)) {
    tmp[ii] = paste(y[x == unqclusts[ii]], collapse = ";")
  }
  return(tmp)
}

# cluster and test robustness
clusts3 = clust.perturb(network, clustering.algorithm = alg, edge.list.format = ef, cluster.format = cf)
```

Note that this is identical to the built-in `hierarchical` method, i.e. `clusts3 = clust.perturb(network, clustering.algorithm = "hierarchical")`.

Exploring some robust and non-robust clusters, we can see that the `repJ` value for the second cluster, `clusts3$repJ[2]`, is relatively high, typically >0.7. This cluster is is "O15143;O15144;O15145;O15511;P59998;P61158;P61160", which completely matches to the Arp2/3 protein complex in CORUM. That is, it completely matches a ground-truth cluster. Therefore, this is a correctly assigned cluster that completely matches a community in the underlying network, and consequentially has high reproducibility.

The `repJ` value for the first cluster `clusts3$repJ[1]` is relatively low, typically <0.3. This cluster is a large agglomeration of proteins from multiple complexes from CORUM, and loosely corresponds to clusters such as the Polycomb repressive complex 1 and SNF2h-cohesin-NuRD complex. Therefore, this cluster represents an incorrect assignment of proteins from multiple ground-truth clusters. Consequentially, small variations in the network can disrupt this assignment, which is represented as low reproducibility as quantified by `repJ`.
