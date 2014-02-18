#' Network enhanced enrichment analysis tool
#' 
#' Functional enrichment anslysis of gene sets by integrating network information
#' 
#' This is a novel functional enrichment analysis tool based on the gene set and network
#' information.
#' 
#' 
#' @import Matrix
#'
#' @export
neeat <- function(gene.set, core.sets, net, bgnet, method = "gene", rho = 0.5, n.perm = 10000, max.depth = 10)
{
  if (method == "gene") {
    net.edges <- net_edges(net)
    result <- sapply(1:dim(core.sets)[2], function(i) neeat_gene(gene.set, core.sets[,i], net.edges, rho, n.perm, max.depth))
  }
  else if (method == "net") {
    core.sets <- core.sets[gene.set,]
    net <- net[gene.set, gene.set]
    result <- sapply(1:dim(core.sets)[2], function(i) neeat_net(core.sets[,i], net, rho, n.perm, max.depth))
  }
  else if (method == "subnet") {
    net[!gene.set, ] <- 0
    net[, !gene.set] <- 0
    net.edges <- net_edges(net)
    bgnet.edges <- net_edges(bgnet)
    result <- sapply(1:dim(core.sets)[2], function(i) neeat_subnet(gene.set, core.sets[,i], net.edges, bgnet.edges, rho, n.perm, max.depth))
  }
  result
}

net_edges <- function(net)
{
  edges <- which(net != 0, arr.ind = T)
  edges <- edges[order(edges[,1]),]
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

neeat_gene <- function(gene.set, core.set, net.edges, rho, n.perm, max.depth)
{
  gene.set <- as.logical(gene.set)
  core.set <- as.logical(core.set)

  depth <- .Call(NE_GetDepths, net.edges$edges, net.edges$index, core.set, max.depth)
  max.depth <- max(depth)

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

neeat_subnet <- function(gene.set, core.set, net.edges, bgnet.edges, rho, n.perm, max.depth)
{
  
}
