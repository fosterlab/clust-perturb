
#' Perturb clusters
#' 
#' Test cluster robustness through random network rewiring
#' 
#' clust.perturb is a general-purprose wrapper for any clustering algorithm. Four default
#' clustering functions are included (MCL, walktrap, hierarchical, and k-medoids) with the
#' option of passing any clustering function. clust.perturb takes input networks as 
#' an unweighted edge list formatted as a 2 column dataframe. Because clustering functions can have
#' different input and output formats, in order to handle arbitrary clustering functions,
#' clust.perturb also takes two conversion functions. 
#' The first, edge.list.format converts the network edge list into the format required by the 
#' clustering algorithm, for example a dist object as required by MCL. The second, cluster.format, 
#' converts the output of the clustering algorithm into a common format, namely a character vector, each
#' element of which is a cluster with semicolon-separated nodes (e.g. c("A;B", "C;D;E"))
#' 
#' clust.perturb returns two metrics for each cluster. repJ measures a cluster's
#' reproducibility, and calculated as the average maximum Jaccard index over noise iterations.
#' fnode, which is calculated for each node in a cluster, counts the frequency with which that 
#' node is reclustered in the closest-match cluster in each noise iteration, divided by the
#' number of iterations.
#' 
#' @param network data frame with two columns. Each row is an edge between two nodes.
#' @param clustering.algorithm a character string specifying one of four clustering
#' algorithms ("mcl", "walktrap", "hierarchical", "k-med"), or a function responsible for
#' clustering
#' @param noise scalar with value between 0 and 1. Specifyies the amount of noise to 
#' add to the network. 0 specifies no noise, and 1 specifies total rewiring. Typical values
#' are between 0.1 and 0.5.
#' @param iters positive integer specifying number of iterations. Typical values are between
#' 3 and 100, with 5-10 iterations often sufficient for estimation.
#' @param edge.list.format NULL or a function that transforms network into format required
#' by clustering.algorithm. If a function, must take exactly one argument.
#' @param cluster.format NULL or a function that transforms output returned by 
#' clustering.algorithm into a character vector, where each element is a cluster whose 
#' format is semicolon-separated nodes. If a function, must take exactly two arguments. The
#' second argument must be a sorted character vector of unique nodes in the original network.
#' @return data frame containing clusters and their repJ scores, fnode scores for each node
#' in each cluster, and the best-matching clusters in each noise iteration.
#' @examples
#' library(igraph)
#' 
#' # walktrap clustering algorithm with random network
#' # make random network
#' network = data.frame(x = sample(1:100, 1000, replace=TRUE), 
#'   y = sample(1:100, 1000, replace=TRUE))
#' # cluster and measure robustness
#' clusts = clust.perturb(network, clustering.algorithm="walktrap")
#' 
#' 
#' # test robustness at low, medium, and high noise levels
#' # demonstrates that an appropriate noise level is one that gives the best resolution of repJ
#' # read network
#' # cluster and measure robustness
#' clusts1 = clust.perturb(network, clustering.algorithm="hierarchical", 
#'   noise=0.001) # low noise
#' clusts2 = clust.perturb(network, clustering.algorithm="hierarchical", 
#'   noise=0.15) # medium noise
#' clusts3 = clust.perturb(network, clustering.algorithm="hierarchical", 
#'   noise=0.75) # high noise
#' # plot
#' plot(sort(clusts1$repJ)) 
#' lines(sort(clusts2$repJ))
#' lines(sort(clusts3$repJ))
#' 
#' 
#' # clustering algorithm with custom conversion functions
#'
#' # use clustering algorithm MCL, explicitly show conversion functions
#' library(MCL)
#' clustalg = function(x) mcl(x, addLoops = FALSE)
#' 
#' # edge.list.format converts dataframe edge.list to adjacency matrix, as required by MCL
#' edgelist.func = function(ints.corum) {
#'   G = graph.data.frame(ints.corum,directed=FALSE)
#'   A = as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)
#' }
#' 
#' # cluster.format converts converts MCL output to character vector of semicolon-separated nodes
#' # cluster.format requires a second argument, unqnodes, which is the sorted vector of unique 
#' # nodes in the network, i.e. unqnodes = unique(c(network[,1], network[,2]))
#' clust.func = function(tmp, unqnodes) {
#'   tmp = tmp$Cluster
#'   clusts = character()
#'   unqclusts = unique(tmp)
#'   for (ii in 1:length(unqclusts)) {
#'     I = tmp == unqclusts[ii]
#'     if (sum(I)<3) next
#'     clusts[ii] = paste(unqnodes[I], collapse = ";")
#'   }
#'   clusts = clusts[!clusts==""]
#'   clusts = clusts[!is.na(clusts)]
#'   return(clusts)
#' }
#' # cluster and measure robustness
#' clusts = clust.perturb(network, clustering.algorithm=clustalg, 
#'   edge.list.format=edgelist.func, cluster.format=clust.func)
#' @export


clust.perturb = function(network, 
                         clustering.algorithm, 
                         noise = 0.1, 
                         iters = 3, 
                         edge.list.format = NULL,
                         cluster.format = NULL) {
  
  if (is.character(clustering.algorithm)) {
    clustering.algorithm = tolower(clustering.algorithm)
    
    if (!clustering.algorithm %in% c("mcl", "hierarchical", "walktrap", "k-med")) {
      stop("clustering.algorithm must either be a function or one of four character strings ('mcl', 'hierarchical', 'walktrap', 'k-med')")
    }
    
    if (!is.null(edge.list.format) | !is.null(cluster.format)) {
      warning(paste("default clustering algorithm specified (", clustering.algorithm, ")"))
      warning(paste("using default edge.list.format and cluster.format functions for ", clustering.algorithm))
    }
    
    if ((clustering.algorithm) == "mcl") {
      tmp = function(x) mymcl(x, infl = 2)
      cluster.format = mcl.cluster.format
      edge.list.format = mcl.edge.list.format
    } else if ((clustering.algorithm) == "hierarchical") {
      tmp = function(x) stats::cutree(stats::hclust(x, method="average"), k=50)
      cluster.format = hierarch.cluster.format
      edge.list.format = hierarch.edge.list.format
    } else if ((clustering.algorithm) == "walktrap") {
      tmp = igraph::walktrap.community
      cluster.format = NULL
      edge.list.format = function(x) igraph::graph_from_edgelist(as.matrix(x), directed = F)
    } else if ((clustering.algorithm) == "k-med") {
      tmp = function(x) cluster::pam(x, 50)
      cluster.format = pam.cluster.format
      edge.list.format = pam.edge.list.format
    }
    clustering.algorithm = tmp
  } 
  
  if (!is.function(clustering.algorithm)) {
    stop("clustering.algorithm must be a character (one of either 'mcl', 'hierarchical', 'walktrap', or 'k-med') or a function")
  }
  
  # cluster without noise
  cc = 0
  unqprots = unique(c(network[,1], network[,2]))
  network.input = network
  if (!is.null(edge.list.format)) network.input = edge.list.format(network)
  tmp = clustering.algorithm(network.input)
  if (!is.null(cluster.format)) tmp = cluster.format(tmp, unqprots)
  # store clusters
  clusters0 = data.frame(cluster = character(length(tmp)),
                         repJ = rep(NA, length(tmp)),
                         fnode = rep(NA, length(tmp)),
                         best.match = rep(NA, length(tmp)),
                         stringsAsFactors = F)
  for (jj in 1:length(tmp)) clusters0$cluster[jj] = paste(tmp[[jj]], collapse=";")
  rm(tmp)

  # cluster with noise
  clusters.noise = data.frame(iter = numeric(10^6),
                              noise_mag = numeric(10^6),
                              cluster = character(10^6), stringsAsFactors = F)
  cc = 0
  for (iter in 1:iters) {
    for (ii in 1:length(noise)) {
      #print(paste("clustering iter",iter,"at noise=", noise[ii]))
      
      # add noise to network
      ints.shuffle = shufflecorum(network, noise[ii])
      unqprots = unique(c(ints.shuffle[,1], ints.shuffle[,2]))
      
      # transform network to required format (if needed)
      if (!is.null(edge.list.format)) ints.shuffle = edge.list.format(ints.shuffle)

      # cluster
      these.clusters = clustering.algorithm(ints.shuffle)
      
      # transform clusters to list (if needed)
      if (!is.null(cluster.format)) these.clusters = cluster.format(these.clusters, unqprots)
      
      # store these.clusters
      for (jj in 1:length(these.clusters)) {
        cc = cc+1
        clusters.noise$iter[cc] = iter
        clusters.noise$noise_mag[cc] = noise[ii]
        clusters.noise$cluster[cc] = paste(these.clusters[[jj]], collapse=";")
      }
    }
  }
  clusters.noise = clusters.noise[1:cc,]
  
  # calculate repJ for every cluster in clusters0
  for (ii in 1:nrow(clusters0)) {
    #print(paste("cluster", ii))
    tmp.j = numeric(iters * length(noise))
    best.match = character(iters * length(noise))
    cc = 0
    for (jj in 1:iters) {
      for (kk in 1:length(noise)) {
        noise.clusters = clusters.noise$cluster[clusters.noise$iter==jj & 
                                                  clusters.noise$noise_mag==noise[kk]]
        cc = cc+1
        tmp = calcJ(clusters0$cluster[ii], noise.clusters)
        tmp.j[cc] = tmp[[1]]
        best.match[cc] = tmp[[2]]
      }
    }
    clusters0$repJ[ii] = mean(tmp.j, na.rm=T)
    clusters0$best.match[ii] = paste(best.match, collapse = " | ")
  }
  
  # calculate fnode for every cluster in clusters0
  clusters0$fnode = character(nrow(clusters0))
  for (jj in 1:nrow(clusters0)) {
    noiseclusts = unlist(strsplit(clusters0$best.match[jj], " | ", fixed=T))
    prots = unlist(strsplit(clusters0$cluster[jj], ";"))
    fnode = numeric(length(prots))
    for (kk in 1:length(fnode)) {
      fnode[kk] = sum(grepl(prots[kk], noiseclusts))
    }
    clusters0$fnode[jj] = paste(fnode/length(noiseclusts), collapse = ";")
  }
  
  return(clusters0)
}

