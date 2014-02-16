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
neeat <- function(gene.set, core.sets, net, rho = 0.5, n.perm = 10000)
{
  net.edges <- which(net != 0, arr.ind = T)
  net.edges <- net.edges[order(net.edges[,1]),]
  net.index <- findInterval(0:dim(core.sets)[1], net.edges[,1])
  sapply(1:dim(core.sets)[2], function(i) neeat_g(gene.set, core.sets[,i], net.edges, net.index, rho, n.perm))
}
 
neeat_g <- function(gene.set, core.set, net.edges, net.index, rho = 0.5, n.perm = 10000)
{
  gene.set <- as.logical(gene.set)
  core.set <- as.logical(core.set)

  depth <- .Call(NE_Depths, net.edges, net.index, core.set)
  max.depth <- max(depth)
  n.depth <- sapply(-1:max.depth, function(d) sum(depth == d))
  w.depth <- c(0, rho^(0:max.depth))
  
  raw.depth <- sapply(-1:max.depth, function(d) sum(depth[gene.set] == d))
  raw.score <- sum(w.depth * raw.depth)

  perm.depth <- rmultihyper(n.perm, n.depth, sum(gene.set))
  perm.score <- colSums(w.depth * perm.depth)

  avg.score <- mean(perm.score)
  var.score <- var(perm.score)
  z.score <- (raw.score - avg.score) / sqrt(var.score)
  p.value <- sum(perm.score >= raw.score) / n.perm
  c(z.score=z.score, p.value=p.value, raw.score=raw.score, 
    avg.score=avg.score, var.score=var.score)
}
