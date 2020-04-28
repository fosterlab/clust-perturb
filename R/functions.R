#' @export

# dependencies
require(igraph)
require(readr)
require(Matrix)
require(cluster)

# functions

calcJ = function(this.cluster, clusters) {
  # J, Ji
  # assignment reproducibility
  # essentially the maximum Jaccard index
  #
  # this.cluster = single string, with semicolon-delimited IDs
  # clusters = vector of strings, all semicolon-delimited IDs
  
  this.cluster = unlist(strsplit(this.cluster, ";"))
  JJ = numeric(length(clusters))
  for (ii in 1:length(clusters)) {
    that.cluster = unlist(strsplit(clusters[ii], ";"))
    JJ[ii] = length(intersect(this.cluster, that.cluster)) / 
      length(unique(c(this.cluster, that.cluster)))
  }
  Ji = max(JJ, na.rm=T)
  tmp = list(Ji = Ji, best.cluster = clusters[which.max(JJ)])
  
  return(tmp)
}

shufflecorum = function(ints.corum, ff){
  unqprots = sort(unique(c(ints.corum[,1], ints.corum[,2])))
  unqprots = unqprots[!unqprots==""]
  
  #
  ints.shuffle = ints.corum
  N.replace = round(nrow(ints.corum) * ff)
  I.replace = sample(nrow(ints.corum), N.replace)
  
  # indices of unshuffled
  ia0 = match(ints.corum[,1], unqprots)
  ib0 = match(ints.corum[,2], unqprots)
  
  # indices of shuffled
  ia = ia0
  ib = ib0
  ia[I.replace] = sample(length(unqprots), N.replace, replace = T)
  ib[I.replace] = sample(length(unqprots), N.replace, replace = T)
  
  # ensure no self-interactions
  while (sum(ia==ib)>0) {
    I = ia==ib
    ib[I] = sample(length(unqprots), sum(I), replace = T)
  }
  
  #
  ints.shuffle[,1] = unqprots[ia]
  ints.shuffle[,2] = unqprots[ib]
  
  # quality control: make sure you shuffled N.replace interactions
  N.diff = sum(!ia0==ia | !ib0==ib)
  #if (abs(N.replace - N.diff)>5) print(paste("shuffling missed", N.replace-N.diff, "interactions"))
  
  return(ints.shuffle) 
}


pam.edge.list.format = function(ints) {
  unqprots = unique(c(ints[,1], ints[,2]))
  nn = length(unqprots)
  
  # create adjacency matrix in the style of `dist` object
  # MAKE SURE YOU'RE GETTING THIS RIGHT!!!
  I.row = match(ints[,1], unqprots) # row
  I.col = match(ints[,2], unqprots) # column
  
  unqprots = unique(c(ints[,1], ints[,2]))
  nn = length(unqprots)
  
  # create adjacency matrix in the style of `dist` object
  # MAKE SURE YOU'RE GETTING THIS RIGHT!!!
  I.row = match(ints[,1], unqprots) # row
  I.col = match(ints[,2], unqprots) # column
  
  I.fill = numeric(length(I.row))
  for (ii in 1:length(I.row)) {
    a = I.row[ii]
    b = I.col[ii]
    # ensure I.col < I.row, i.e. upper triangular
    if (I.row[ii] < I.col[ii]) {
      a = I.col[ii]
      b = I.row[ii]
    }
    
    I.fill[ii] = a - 1
    if (b>1) {
      colsum = 0
      for (jj in 1:(b-1)) {
        colsum = colsum + (nn-jj) - 1
      }
      I.fill[ii] = colsum + a - 1
    }
  }
  
  # dummy dist object
  x = matrix(stats::runif(nn * 10), nrow = nn, ncol=10)
  d = stats::dist(x)
  attr(d, 'Upper') = T
  d[1:length(d)] = 1
  d[I.fill] = 0
  
  return(d)
}

pam.cluster.format = function(clusts, unqprots) {
  # compile `clusts` into lists of proteins
  unqclusts = unique(clusts$clustering)
  Nmembers = numeric(length(unqclusts))
  clusts.prots = character(length(unqclusts))
  for (ii in 1:length(unqclusts)) {
    I = clusts$cluster==unqclusts[ii]
    clusts.prots[ii] = paste(unqprots[I], collapse=";")
    Nmembers[ii] = sum(I)
  }
  #clusts.prots = clusts.prots[Nmembers>=3]
  #clusts.prots = as.list(clusts.prots)
  
  return(clusts.prots)
}


mcl.edge.list.format = function(ints.corum) {
  G = igraph::graph.data.frame(ints.corum,directed=FALSE)
  A = igraph::as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)
}

mcl.cluster.format = function(tmp, unqnodes) {
  tmp = tmp$Cluster
  clusts = character()
  unqclusts = unique(tmp)
  for (ii in 1:length(unqclusts)) {
    I = tmp == unqclusts[ii]
    #if (sum(I)<3) next
    clusts[ii] = paste(unqnodes[I], collapse = ";")
  }
  clusts = clusts[!clusts==""]
  clusts = clusts[!is.na(clusts)]
  return(clusts)
}

# hierachical
hierarch.edge.list.format = function(ints) {
  return(pam.edge.list.format(ints))
}

hierarch.cluster.format = function(tmp, unqprots) {
  clusts = character()
  unqclusts = unique(tmp)
  for (ii in 1:length(unqclusts)) {
    I = tmp == unqclusts[ii]
    clusts[ii] = paste(unqprots[I], collapse = ";")
  }
  clusts = clusts[!clusts==""]
  clusts = clusts[!is.na(clusts)]
  return(clusts)
}
