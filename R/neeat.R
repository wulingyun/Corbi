#' Network enhanced enrichment analysis tool
#' 
#' Functional enrichment anslysis of gene sets by integrating network information
#' 
#' This is a novel functional enrichment analysis tool based on the gene set and network
#' information.
#' 
#' @param core.sets Logical matrix indicated the core genes associated with specific functions or pathways
#' @param gene.set Logical vector indicated the gene set for evaluating
#' @param net The adjacent matrix of network
#' @param subnet The adjacent matrix of sub-network for evaluating in the NEEAT-subnet model
#' @param method A string indicated the NEEAT model, including "gene", "net" and "subnet"
#' @param rho The weight parameter for depths
#' @param n.perm The number of permutations for calculating p-values
#' @param max.depth Integer for the maximum depth considered in the NEEAT models
#' @param n.cpu The number of CPUs/cores used in the parallel computation
#' @return This function will return a matrix of same columns as \code{core.sets}, and each column
#' containing the following components for the correponding core gene set \code{core.sets[,i]}:
#' \itemize{
#'   \item \code{z.score} The Z-score for \code{gene.set}
#'   \item \code{p.value} The statistic significance for \code{gene.set} under specified NEEAT model
#'   \item \code{raw.score} The raw score for \code{gene.set} under specified NEEAT model
#'   \item \code{avg.score} The average score for \code{n.perm} random permutations
#'   \item \code{var.score} The variance of scores for \code{n.perm} random permutations
#' }
#' 
#' @import Matrix parallel
#'
#' @export
neeat <- function(core.sets, gene.set = NULL, net = NULL, subnet = NULL, method = "gene", rho = 0.5, n.perm = 10000, max.depth = 10, n.cpu = 1)
{
  if (n.cpu > 1) {
    cl <- makeCluster(n.cpu)
    jobs <- clusterSplit(cl, 1:dim(core.sets)[2])
    fun <- function (x, ...) neeat_internal(core.sets[, x], ...)
    result <- clusterApply(cl, jobs, fun, gene.set, net, subnet, method, rho, n.perm, max.depth)
    stopCluster(cl)
    matrix(unlist(result), nrow = 5, dimnames = list(c("z.score", "p.value", "raw.score", "avg.score", "var.score")))
  }
  else
    neeat_internal(core.sets, gene.set, net, subnet, method, rho, n.perm, max.depth)
}

neeat_internal <- function(core.sets, gene.set, net, subnet, method, rho, n.perm, max.depth)
{
  if (method == "gene" && !is.null(gene.set)) {
    gene.set <- as.logical(gene.set)
    if (is.null(net)) {
      n.gene <- length(gene.set)
      net <- sparseMatrix(n.gene, n.gene, x = 0)
    }
    net.edges <- net_edges(net)
    sapply(1:dim(core.sets)[2], function(i) neeat_gene(core.sets[,i], gene.set, net.edges, rho, n.perm, max.depth))
  }
  else if (method == "net" && !is.null(net)) {
    if (!is.null(gene.set)) {
      core.sets <- core.sets[gene.set, ]
      net <- net[gene.set, gene.set]
    }
    sapply(1:dim(core.sets)[2], function(i) neeat_net(core.sets[,i], net, rho, n.perm, max.depth))
  }
  else if (method == "subnet" && !is.null(gene.set) && !is.null(net)) {
    gene.set <- as.logical(gene.set)
    net.edges <- net_edges(net)
    if (is.null(subnet)) {
      subnet.edges <- net_edges(net, gene.set)
    }
    else {
      subnet.edges <- net_edges(subnet)
    }
    sapply(1:dim(core.sets)[2], function(i) neeat_subnet(core.sets[,i], gene.set, net.edges, subnet.edges, rho, n.perm, max.depth))
  }
  else {
    stop("Incorrect parameters!")
  }
}

net_edges <- function(net, gene.set = NULL)
{
  edges <- which(net != 0, arr.ind = T)
  if (!is.null(gene.set))
    edges <- matrix(edges[gene.set[edges[,1]] & gene.set[edges[,2]], ], ncol=2)
  edges <- edges[order(edges[,1]), ]
  index <- findInterval(0:dim(net)[1], edges[,1])
  list(edges=edges, index=index)
}

neeat_score <- function(w.depth, n.depth, raw.depth, n.perm)
{
  raw.score <- sum(w.depth * raw.depth)
  
  perm.depth <- rmultihyper(n.perm, n.depth, sum(raw.depth))
  perm.score <- colSums(w.depth * perm.depth)
  
  avg.score <- mean(perm.score)
  var.score <- var(perm.score)
  z.score <- (raw.score - avg.score) / sqrt(var.score)
  p.value <- sum(perm.score >= raw.score) / n.perm

  c(z.score=z.score, p.value=p.value, raw.score=raw.score, 
    avg.score=avg.score, var.score=var.score)
}

neeat_gene <- function(core.set, gene.set, net.edges, rho, n.perm, max.depth)
{
  core.set <- as.logical(core.set)

  depth <- .Call(NE_GetDepths, net.edges$edges, net.edges$index, core.set, max.depth)
  max.depth <- max(0, depth)

  w.depth <- c(0, rho^(0:max.depth))
  n.depth <- .Call(NE_CountDepths, depth, max.depth)
  raw.depth <- .Call(NE_CountDepths, depth[gene.set], max.depth)
  
  neeat_score(w.depth, n.depth, raw.depth, n.perm)
}

neeat_net <- function(core.set, net, rho, n.perm, max.depth)
{
  core.set <- as.logical(core.set)

  w.depth <- rho^(0:2)

  nc <- sum(core.set)
  nn <- length(core.set) - nc
  n.depth <- c(nc*(nc-1)/2, nc*nn, nn*(nn-1)/2)

  nc <- nnzero(net[core.set, core.set])
  nn <- nnzero(net[!core.set, !core.set])
  raw.depth <- c(nc/2, (nnzero(net)-nc-nn)/2, nn/2)

  neeat_score(w.depth, n.depth, raw.depth, n.perm)
}

edge_depth <- function(node.depth, net.edges)
{
  edges <- net.edges$edges
  edges <- matrix(edges[edges[,1] < edges[,2], ], ncol=2)
  depth <- node.depth[edges[,1]] + node.depth[edges[,2]]
  depth[depth < 0] <- -1
  depth
}

neeat_subnet <- function(core.set, gene.set, net.edges, subnet.edges, rho, n.perm, max.depth)
{
  core.set <- as.logical(core.set)
  
  node.depth <- .Call(NE_GetDepths, net.edges$edges, net.edges$index, core.set, max.depth)
  node.depth[node.depth < 0] <- -length(gene.set)
  
  depth <- edge_depth(node.depth, net.edges)
  max.depth <- max(0, depth)
  
  w.depth <- c(0, rho^(0:max.depth))
  n.depth <- .Call(NE_CountDepths, depth, max.depth)
  raw.depth <- .Call(NE_CountDepths, edge_depth(node.depth, subnet.edges), max.depth)
  
  neeat_score(w.depth, n.depth, raw.depth, n.perm)
}
