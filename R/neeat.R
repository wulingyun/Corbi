#' Network enhanced enrichment analysis tool
#' 
#' Functional enrichment anslysis of gene sets by integrating network information
#' 
#' This is a novel functional enrichment analysis tool based on the gene set and network
#' information.
#' 
#' 
#'
#' @export
neeat_g <- function(gene.set, core.set, net, rho = 0.5, n.perm = 10000)
{
  gene.set <- as.logical(gene.set)
  core.set <- as.logical(core.set)

  depth <- .Call(NE_Depths, net, core.set)
  max.depth <- max(depth)
  n.depth <- sapply(-1:max.depth, function(d) sum(depth == d))
  w.depth <- c(0, rho^(0:max.depth))
  
  raw.depth <- sapply(-1:max.depth, function(d) sum(depth[gene.set] == d))
  raw.score <- sum(w.depth * raw.depth)

  perm.depth <- rmultihyper(n.perm, n.depth, sum(gene.set))
  perm.score <- colSums(w.depth * perm.depth)

  exp.score <- mean(perm.score)
  g.score <- raw.score / exp.score
  p.value <- sum(exp.score >= raw.score)
  list(g.score=g.score, p.value=p.value, raw.score=raw.score, exp.score=exp.score)
}
