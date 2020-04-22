# clust-perturb

clust.perturb is an R package designed to test the robustness of graph-based clusters by randomly rewiring a portion of the underlying network and reclustering. Robust clusters “stick together” even after network rewiring, while non-robust clusters do not. Robustness is quantified using maximum Jaccard index.

## System requirements

Augur relies on functions from the following R packages:

```
	igraph (>= 1.2.4.2)
	readr (>= 1.3.1)
	Matrix (>= 1.2-17)
	cluster (>= 2.1.0)
```

## Installation

To install clust.perturb, first install the devtools package, if it is not already installed: 

```r
> install.packages("devtools") 
```

Then, install clust.perturb from GitHub: 

```r
> devtools::install_github("GregStacey/clust-perturb")
```

This should take no more than a few minutes.
