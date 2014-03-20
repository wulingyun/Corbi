#' Network enhanced enrichment analysis tool
#' 
#' Functional enrichment anslysis of gene sets by integrating network information
#' 
#' This is a novel functional enrichment analysis tool based on the gene set and network
#' information.
#' 
#' @param core.sets Logical matrix indicated the core genes associated with specific functions or pathways.
#' @param gene.set Logical vector indicated the gene set for evaluating.
#' @param net The adjacent matrix of network.
#' @param subnet The adjacent matrix of sub-network for evaluating in the NEEAT-subnet model.
#' @param method A string indicated the NEEAT model, including "gene", "net" and "subnet".
#' Use "hyper" for traditional hypergeometric test, in which the network information is ignored.
#' @param rho The weight parameter for depths.
#' @param max.depth Integer for the maximum depth considered in the NEEAT models.
#' @param n.perm The number of permutations for calculating p-values.
#' @param use.multinom Logical variable indicated whether use \code{\link{rmultinom}} to 
#' approximate \code{\link{rmultihyper}}.
#' @param z.threshold The threshold for filtering out small Z-scores. The p-values are calculated only
#' for the Z-scores \code{>= z.threshold}, otherwise p-values are set as 1.
#' @param adjust.p The method for adjusting p-values for multiple testing. Use "none" for bypassing.
#' See \code{\link{p.adjust}} for available methods.
#' @param n.cpu The number of CPUs/cores used in the parallel computation.
#' @param batch.size The desired size of batches in the parallel computation.
#' 
#' @return This function will return a matrix of same columns as \code{core.sets}, and each column
#' containing the following components for the correponding core gene set \code{core.sets[,i]}:
#'   \item{\code{z.score}}{The Z-score for \code{gene.set}}
#'   \item{\code{p.value}}{The statistic significance for \code{gene.set} under specified NEEAT model}
#'   \item{\code{raw.score}}{The raw score for \code{gene.set} under specified NEEAT model}
#'   \item{\code{avg.score}}{The average score for random permutations of \code{gene.set}}
#'   \item{\code{var.score}}{The variance of scores for random permutations of \code{gene.set}}
#'
#' @seealso \code{\link{get_core_sets}}
#' 
#' @import Matrix parallel
#'
#' @export
neeat <- function(core.sets, gene.set = NULL, net = NULL, subnet = NULL, method = "gene", 
                  rho = 0.5, max.depth = 10, n.perm = 10000, use.multinom = FALSE,
                  z.threshold = 2.0, adjust.p = "BH", n.cpu = 1, batch.size = 5000)
{
  options <- list(rho = rho,
                  max.depth = max.depth,
                  n.perm = n.perm,
                  use.multinom = use.multinom,
                  z.threshold = z.threshold)
  
  if (n.cpu > 1) {
    if (.Platform$OS.type == "windows")
      cl <- makeCluster(n.cpu)
    else
      cl <- makeForkCluster(n.cpu)
    n.cl <- length(cl)
    n.job <- dim(core.sets)[2]
    n.batch <- max(1, round(n.job / (batch.size * n.cl))) * n.cl
    jobs <- splitIndices(n.job, n.batch)
    fun <- function (x, ...) neeat_internal(core.sets[, x], ...)
    result <- clusterApply(cl, jobs, fun, gene.set, net, subnet, method, options)
    stopCluster(cl)
    result <- matrix(unlist(result), nrow = 5)
  }
  else {
    result <- neeat_internal(core.sets, gene.set, net, subnet, method, options)
  }
  result[2,] <- p.adjust(result[2,], method = adjust.p)
  rownames(result) <- c("z.score", "p.value", "raw.score", "avg.score", "var.score")
  result
}

neeat_internal <- function(core.sets, gene.set, net, subnet, method, options)
{
  if (method == "gene" && !is.null(gene.set)) {
    gene.set <- as.logical(gene.set)
    if (is.null(net)) {
      n.gene <- length(gene.set)
      net <- sparseMatrix(n.gene, n.gene, x = 0)
    }
    net.edges <- net_edges(net)
    sapply(1:dim(core.sets)[2], function(i) neeat_gene(core.sets[,i], gene.set, net.edges, options))
  }
  else if (method == "net" && !is.null(net)) {
    if (!is.null(gene.set)) {
      core.sets <- core.sets[gene.set, ]
      net <- net[gene.set, gene.set]
    }
    sapply(1:dim(core.sets)[2], function(i) neeat_net(core.sets[,i], net, options))
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
    sapply(1:dim(core.sets)[2], function(i) neeat_subnet(core.sets[,i], gene.set, net.edges, subnet.edges, options))
  }
  else if (method == "hyper" && !is.null(gene.set)) {
    gene.set <- as.logical(gene.set)
    N <- length(gene.set)
    n <- sum(gene.set)
    M <- colSums(core.sets)
    m <- colSums(core.sets & gene.set)
    sapply(1:dim(core.sets)[2], function(i) neeat_hyper(N, n, M[i], m[i]))
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

neeat_score <- function(w.depth, n.depth, raw.depth, options)
{
  p <- n.depth / sum(n.depth)
  w <- w.depth * p
  n <- sum(raw.depth)
  N <- sum(n.depth)
  cov <- -tcrossprod(w)
  diag(cov) <- w.depth * w * (1-p)

  raw.score <- sum(w.depth * raw.depth)
  avg.score <- sum(w) * n
  var.score <- sum(cov) * n

  if (options$use.multinom) {
    simulate <- rmultinom
  }
  else {
    var.score <- var.score * (N - n) / (N - 1)
    simulate <- rmultihyper
  }

  if (var.score > 0)
    z.score <- (raw.score - avg.score) / sqrt(var.score)
  else
    z.score <- 0

  if (z.score >= options$z.threshold) {
    perm.depth <- simulate(options$n.perm, sum(raw.depth), n.depth)
    perm.score <- colSums(w.depth * perm.depth)
    p.value <- sum(perm.score >= raw.score) / options$n.perm
  }
  else {
    p.value <- 1.0
  }

  c(z.score=z.score, p.value=p.value, raw.score=raw.score, 
    avg.score=avg.score, var.score=var.score)
}

neeat_gene <- function(core.set, gene.set, net.edges, options)
{
  core.set <- as.logical(core.set)

  depth <- .Call(NE_GetDepths, net.edges$edges, net.edges$index, core.set, options$max.depth)
  max.depth <- max(0, depth)

  w.depth <- c(0, options$rho^(0:max.depth))
  n.depth <- .Call(NE_CountDepths, depth, max.depth)
  raw.depth <- .Call(NE_CountDepths, depth[gene.set], max.depth)
  
  neeat_score(w.depth, n.depth, raw.depth, options)
}

neeat_net <- function(core.set, net, options)
{
  core.set <- as.logical(core.set)
  
  max.depth <- min(2, options$max.depth)
  w.depth <- c(0, options$rho^(0:max.depth))
  
  nc <- sum(core.set)
  nn <- length(core.set) - nc
  n.depth <- c(nc*(nc-1)/2, nc*nn, nn*(nn-1)/2)
  n.depth <- c(sum(n.depth[-(0:max.depth+1)]), n.depth[0:max.depth+1])

  nc <- nnzero(net[core.set, core.set])
  nn <- nnzero(net[!core.set, !core.set])
  raw.depth <- c(nc/2, (nnzero(net)-nc-nn)/2, nn/2)
  raw.depth <- c(sum(raw.depth[-(0:max.depth+1)]), raw.depth[0:max.depth+1])
  
  neeat_score(w.depth, n.depth, raw.depth, options)
}

edge_depth <- function(node.depth, net.edges, max.depth)
{
  edges <- net.edges$edges
  edges <- matrix(edges[edges[,1] < edges[,2], ], ncol=2)
  depth <- node.depth[edges[,1]] + node.depth[edges[,2]]
  depth[depth < 0] <- -1
  depth[depth > max.depth] <- -1
  depth
}

neeat_subnet <- function(core.set, gene.set, net.edges, subnet.edges, options)
{
  core.set <- as.logical(core.set)
  
  node.depth <- .Call(NE_GetDepths, net.edges$edges, net.edges$index, core.set, (options$max.depth+1)/2)
  node.depth[node.depth < 0] <- -Inf
  
  depth <- edge_depth(node.depth, net.edges, options$max.depth)
  max.depth <- max(0, depth)
  
  w.depth <- c(0, options$rho^(0:max.depth))
  n.depth <- .Call(NE_CountDepths, depth, max.depth)
  raw.depth <- .Call(NE_CountDepths, edge_depth(node.depth, subnet.edges, max.depth), max.depth)
  
  neeat_score(w.depth, n.depth, raw.depth, options)
}

neeat_hyper <- function(N, n, M, m)
{
  avg.score <- n * M / N
  var.score <- n * (M / N) * ((N - M) / N) * ((N - n) / (N - 1))
  
  if (var.score > 0)
    z.score <- (m - avg.score) / sqrt(var.score)
  else
    z.score <- 0
  
  p.value <- 1 - phyper(m-1, M, N-M, n)
  
  c(z.score=z.score, p.value=p.value, raw.score=m, 
    avg.score=avg.score, var.score=var.score)
}



#' Extract core sets information from GO annotation
#'
#' Generate the core sets matrix from an object of class "Go3AnnDbBimap" in Bioconductor annotation package.
#'
#' This function generates the \code{core.sets} matrix required by \code{\link{neeat}} method.
#' 
#' @param go.map An object of class "Go3AnnDbBimap".
#' @param evidence Vector of string to filter the GO annotations, could be "ALL" or a set of evidence codes.
#' @param category Vector of string to filter the GO categories, could be "ALL" or a set of GO categories.
#' @param gene.set Vector of string to filter the genes.
#' 
#' @return This function returns a sparse matrix as required by \code{\link{neeat}} method.
#' 
#' @seealso \code{\link{neeat}}
#'
#' @examples
#' 
#' \dontrun{
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("org.Hs.eg.db")
#' library(org.Hs.eg.db)
#' x <- get_core_sets(org.Hs.egGO2ALLEGS)
#' }
#'
#' @import Matrix
#' 
#' @export
get_core_sets <- function(go.map, evidence = "ALL", category = "ALL", gene.set = NULL)
{
  term.table <- AnnotationDbi::toTable(go.map)
  term.table[,1] <- toupper(term.table[,1])
  
  if (evidence != "ALL")
    term.table <- term.table[term.table[,3] %in% evidence, ]
  
  if (category != "ALL")
    term.table <- term.table[term.table[,4] %in% category, ]
  
  if (!is.null(gene.set)) {
    term.table <- term.table[term.table[,1] %in% gene.set, ]
    all.gene <- gene.set    
  }
  else
    all.gene <- unique(term.table[,1])
  
  all.term <- unique(term.table[,2])
  
  gene.id <- seq_along(all.gene)
  names(gene.id) <- all.gene
  term.id <- seq_along(all.term)
  names(term.id) <- all.term
  
  sparseMatrix(gene.id[term.table[,1]], term.id[term.table[,2]], x = T, dims = c(length(all.gene), length(all.term)), dimnames = list(all.gene, all.term))
}
