
#' Perturb clusters
#' 
#' Network to be clustered is formated as a data frame edge list (edge.list). Users can 
#' pass a custom function edge.list.format to handle different formats required by
#' clustering.algorithm. edge.list.format must transform a data frame edge list (two 
#' columns unweighted, or three columsn weighted) into the required format.
#' 
#' Assumes clustering.algorithm returns a list of clustered nodes. If clustering.algorithm
#' returns a different format, users can pass a custom function cluster.format that
#' transforms the output of clustering.algorithm into a list of clustered nodes.
#' 
#' @param edge.list data frame with two or three columns, with first two columns providing 
#' the list of edges and an optional third column giving edge weight.
#' @param clustering.algorithm either a character string specifying one of four clustering
#' algorithms ("mcl", "walktrap", "hierarchical", "k-med") or a function responsible for
#' clustering
#' @param noise vector with values between 0 and 1 specifying the amount(s) of noise to 
#' add to the network.
#' @param iters positive integer specifying number of iterations.
#' @param edge.list.format optional function that transforms edge list into format required
#' by clustering.algorithm.
#' @param cluster.format optional function that transforms output returned by 
#' clustering.algorithm into a data frame with three or four columns, where columns one 
#' and two are an edge list, column three is an id specifying the cluster assignment, and 
#' @return data frame containing all clusters and their density score.
#' @examples
#' first example
#' 
#' second example


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
      tmp = function(x) cutree(hclust(x, method="average"), k=50)
      cluster.format = hierarch.cluster.format
      edge.list.format = hierarch.edge.list.format
    } else if ((clustering.algorithm) == "walktrap") {
      tmp = walktrap.community
      cluster.format = NULL
      edge.list.format = function(x) graph_from_edgelist(as.matrix(x), directed = F)
    } else if ((clustering.algorithm) == "k-med") {
      tmp = function(x) pam(x, 50)
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
  
  # calculate J for every cluster in clusters0
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

