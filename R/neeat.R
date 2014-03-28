#' Network enhanced enrichment analysis tool
#' 
#' Functional enrichment anslysis of gene sets by integrating network information
#' 
#' This is a novel functional enrichment analysis tool based on the gene set and network
#' information.
#' 
#' @param core.sets Logical matrix indicated the core genes associated with specific functions or pathways.
#' The rows correspond to genes, while the columns represent the functions or pathways.
#' @param gene.sets Logical vector or matrix indicated the gene set(s) for evaluating.
#' The rows correspond to genes, while the columns represent the gene sets.
#' @param net The adjacent matrix of network.
#' @param subnet The adjacent matrix of sub-network for evaluating in the NEEAT-subnet model.
#' @param depths The node depths defined in the NEEAT model, which can be pre-computed by calling
#' \code{\link{neeat_depths}} in advance in order to reduce the computation time. If not provided,
#' it will be computed on-demand.
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
#' @return This function will return a 3-dimensional array of dimensions \code{c(5, dim(gene.sets)[2], dim(core.sets)[2])},
#' and each column \code{[,i,j]} containing the following components for the correponding gene set \code{gene.sets[,i]}
#' and core set \code{core.sets[,j]}:
#'   \item{\code{z.score}}{The Z-score for the gene set}
#'   \item{\code{p.value}}{The statistic significance for the gene set under specified NEEAT model}
#'   \item{\code{raw.score}}{The raw score for the gene set under specified NEEAT model}
#'   \item{\code{avg.score}}{The average score for random permutations of the gene set}
#'   \item{\code{var.score}}{The variance of scores for random permutations of the gene set}
#'
#' @seealso \code{\link{get_core_sets}}, \code{\link{neeat_depths}}
#' 
#' @import Matrix parallel
#'
#' @export
neeat <- function(core.sets, gene.sets = NULL, net = NULL, subnet = NULL, depths = NULL,
                  method = "gene", rho = 0.5, max.depth = 10, n.perm = 10000, use.multinom = FALSE,
                  z.threshold = 2.0, verbose = FALSE, adjust.p = "BH", n.cpu = 1, batch.size = 5000)
{
  options <- list(rho = rho,
                  max.depth = max.depth,
                  n.perm = n.perm,
                  use.multinom = use.multinom,
                  z.threshold = z.threshold,
                  verbose = verbose)

  if (is.null(gene.sets))
    gene.sets <- matrix(T, dim(core.sets)[1], 1, dimnames = list(rownames(core.sets), "net"))
  
  if (is.null(dim(core.sets)))
    core.sets <- Matrix(as.logical(core.sets))
  if (is.null(dim(gene.sets)))
    gene.sets <- Matrix(as.logical(gene.sets))
  if (!is.null(depths) && is.null(dim(depths)))
    depths <- Matrix(depths)

  if (n.cpu > 1) {
    if (.Platform$OS.type == "windows")
      cl <- makeCluster(n.cpu)
    else
      cl <- makeForkCluster(n.cpu)
    n.cl <- length(cl)
    n.job <- dim(core.sets)[2]
    n.batch <- max(1, round(n.job / (batch.size * n.cl))) * n.cl
    jobs <- splitIndices(n.job, n.batch)
    result <- clusterApply(cl, jobs, neeat_internal, core.sets, gene.sets, net, subnet, depths, method, options)
    stopCluster(cl)
    result <- unlist(result)
  }
  else {
    result <- neeat_internal(seq_len(dim(core.sets)[2]), core.sets, gene.sets, net, subnet, depths, method, options)
  }
  if (options$verbose)
    result <- array(result, dim = c(5, dim(gene.sets)[2], dim(core.sets)[2]),
                    dimnames <- list(c("z.score", "p.value", "raw.score", "avg.score", "var.score"),
                                     colnames(gene.sets), colnames(core.sets)))
  else
    result <- array(result, dim = c(2, dim(gene.sets)[2], dim(core.sets)[2]),
                    dimnames <- list(c("z.score", "p.value"),
                                     colnames(gene.sets), colnames(core.sets)))
  result[2,,] <- p.adjust(result[2,,], method = adjust.p)
  result
}

neeat_internal <- function(cs.ids, core.sets, gene.sets, net, subnet, depths, method, options)
{
  gs.ids <- seq_len(dim(gene.sets)[2])
  if (method == "gene") {
    if (is.null(net)) {
      n.gene <- dim(gene.sets)[1]
      net <- sparseMatrix(n.gene, n.gene, x = F)
    }
    net.edges <- net_edges(net)
    fun <- function(i)
    {
      cs <- column(core.sets, i)
      if (is.null(depths)) {
        dp <- get_depths(cs, net.edges, options$max.depth)
      }
      else {
        dp <- column(depths, i)
      }
      sapply(gs.ids, function(j) neeat_gene(cs, column(gene.sets, j), net.edges, dp, options))
    }
  }
  else if (method == "net" && !is.null(net)) {
    n.gene <- colSums(gene.sets)
    fun <- function(i)
    {
      cs <- column(core.sets, i)
      sapply(gs.ids, function(j) neeat_net(cs, column(gene.sets, j), net, n.gene[j], options))
    }
  }
  else if (method == "subnet" && !is.null(net)) {
    net.edges <- net_edges(net)
    if (is.null(subnet)) {
      subnet <- net
    }
    subnet.edges <- lapply(gs.ids, function(i) net_edges(subnet, column(gene.sets, i)))
    fun <- function(i)
    {
      cs <- column(core.sets, i)
      if (is.null(depths)) {
        dp <- get_depths(cs, net.edges, options$max.depth)
      }
      else {
        dp <- column(depths, i)
      }
      sapply(gs.ids, function(j) neeat_subnet(cs, column(gene.sets, j), net.edges, subnet.edges[[j]], dp, options))
    }
  }
  else if (method == "hyper") {
    N <- dim(core.sets)[1]
    M <- colSums(core.sets)
    n <- colSums(gene.sets)
    fun <- function(i)
    {
      cs <- column(core.sets, i)
      sapply(gs.ids, function(j) neeat_hyper(N, n[j], M[i], sum(cs & column(gene.sets, j))))
    }
  }
  else {
    stop("Incorrect parameters!")
  }
  sapply(cs.ids, fun)
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

  if (options$verbose)
    c(z.score, p.value, raw.score, avg.score, var.score)
  else
    c(z.score, p.value)
}

neeat_gene <- function(core.set, gene.set, net.edges, depth, options)
{
  if (is.null(depth))
    depth <- .Call(NE_GetDepths, net.edges$edges, net.edges$index, core.set, options$max.depth)
  max.depth <- max(0, depth)

  w.depth <- c(0, options$rho^(0:max.depth))
  n.depth <- .Call(NE_CountDepths, depth, max.depth)
  raw.depth <- .Call(NE_CountDepths, depth[gene.set], max.depth)
  
  neeat_score(w.depth, n.depth, raw.depth, options)
}

neeat_net <- function(core.set, gene.set, net, n.gene, options)
{
  max.depth <- min(2, options$max.depth)
  w.depth <- c(0, options$rho^(0:max.depth))
  
  cs.0 <- core.set & gene.set
  cs.1 <- !core.set & gene.set
  
  nc <- sum(cs.0)
  nn <- n.gene - nc
  n.depth <- c(nc*(nc-1)/2, nc*nn, nn*(nn-1)/2)
  n.depth <- c(sum(n.depth[-(0:max.depth+1)]), n.depth[0:max.depth+1])

  n0 <- nnzero(net, cs.0, cs.0) / 2
  n1 <- (nnzero(net, cs.0, cs.1) + nnzero(net, cs.1, cs.0)) / 2
  n2 <- nnzero(net, cs.1, cs.1) / 2
  raw.depth <- c(n0, n1, n2)
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

neeat_subnet <- function(core.set, gene.set, net.edges, subnet.edges, node.depth, options)
{
  if (is.null(node.depth))
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


#' Compute the node depths
#' 
#' Compute the node depths for given core sets and network
#' 
#' This function calculates the node depths, which can be used for calling \code{\link{neeat}} in the
#' batch mode to reduce the computation time.
#' 
#' @param core.sets Logical matrix indicated the core genes associated with specific functions or pathways.
#' @param net The adjacent matrix of network.
#' @param max.depth Integer for the maximum depth that will be computed.
#' 
#' @return This function returns a matrix with the same dimensions as \code{core.sets}, where the element
#' \code{[i, j]} represents the depth of node \code{i} under the condition (function or pathway) \code{j}.
#' 
#' @seealso \code{\link{neeat}}
#' 
#' @export
neeat_depths <- function(core.sets, net, max.depth)
{
  net.edges <- net_edges(net)
  fun <- function(i) get_depths(column(core.sets, i), net.edges, max.depth)
  sapply(seq_len(dim(core.sets)[2]), fun)
}

get_depths <- function(core.set, net.edges, max.depth)
  .Call(NE_GetDepths, net.edges$edges, net.edges$index, core.set, max.depth)


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
