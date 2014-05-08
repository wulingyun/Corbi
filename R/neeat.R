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
                  method = "gene", rho = 0.001, max.depth = 2, n.perm = 0, use.multinom = FALSE,
                  z.threshold = 0, verbose = FALSE, adjust.p = "BH", n.cpu = 1, batch.size = 5000)
{
  neeat.par <- new.env()
  neeat.par$rho = rho
  neeat.par$max.depth = max.depth
  neeat.par$n.perm = n.perm
  neeat.par$use.multinom = use.multinom
  neeat.par$z.threshold = z.threshold
  neeat.par$verbose = verbose

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
    result <- clusterApply(cl, jobs, neeat_internal, core.sets, gene.sets, net, subnet, depths, method, neeat.par)
    stopCluster(cl)
    result <- unlist(result)
  }
  else {
    result <- neeat_internal(seq_len(dim(core.sets)[2]), core.sets, gene.sets, net, subnet, depths, method, neeat.par)
  }
  if (neeat.par$verbose)
    output.names <- c("z.score", "p.value", "raw.score", "avg.score", "var.score")
  else
    output.names <- c("z.score", "p.value")
  result <- array(result, dim = c(length(output.names), dim(gene.sets)[2], dim(core.sets)[2]),
                  dimnames <- list(output.names, colnames(gene.sets), colnames(core.sets)))
  result[2,,] <- p.adjust(result[2,,], method = adjust.p)
  result
}

neeat_internal <- function(cs.ids, core.sets, gene.sets, net, subnet, depths, method, neeat.par)
{
  gs.ids <- seq_len(dim(gene.sets)[2])
  if (method == "gene") {
    if (is.null(net)) {
      n.gene <- dim(gene.sets)[1]
      net <- sparseMatrix(n.gene, n.gene, x = F)
    }
    net.edges <- net_edges(net)
    neeat_cs <- function(i)
    {
      cs <- column(core.sets, i)
      if (is.null(depths)) {
        node.depth <- get_depths(cs, net.edges, neeat.par$max.depth)
      }
      else {
        node.depth <- column(depths, i)
      }
      max.depth <- max(0, node.depth)
      w.depth <- c(0, neeat.par$rho^(0:max.depth))
      n.depth <- .Call(NE_CountDepths, node.depth, max.depth)
      neeat_score_init(w.depth, n.depth, neeat.par)
      neeat_gs <- function(j)
      {
        raw.depth <- .Call(NE_CountDepths, node.depth[column(gene.sets, j)], max.depth)
        neeat_score(w.depth, n.depth, raw.depth, neeat.par)
      }
      sapply(gs.ids, neeat_gs)
    }
  }
  else if (method == "net" && !is.null(net)) {
    n.gene <- colSums(gene.sets)
    neeat_cs <- function(i)
    {
      cs <- column(core.sets, i)
      sapply(gs.ids, function(j) neeat_net(cs, column(gene.sets, j), net, n.gene[j], neeat.par))
    }
  }
  else if (method == "subnet" && !is.null(net)) {
    net.edges <- net_edges(net)
    if (is.null(subnet)) {
      subnet <- net
    }
    subnet.edges <- lapply(gs.ids, function(i) net_edges(subnet, column(gene.sets, i)))
    neeat_cs <- function(i)
    {
      cs <- column(core.sets, i)
      if (is.null(depths)) {
        node.depth <- get_depths(cs, net.edges, neeat.par$max.depth)
      }
      else {
        node.depth <- column(depths, i)
      }
      node.depth[node.depth < 0] <- -Inf
      edge.depth <- edge_depth(node.depth, net.edges, neeat.par$max.depth)
      max.depth <- max(0, edge.depth)
      w.depth <- c(0, neeat.par$rho^(0:max.depth))
      n.depth <- .Call(NE_CountDepths, edge.depth, max.depth)
      neeat_score_init(w.depth, n.depth, neeat.par)
      neeat_gs <- function(j)
      {
        raw.depth <- .Call(NE_CountDepths, edge_depth(node.depth, subnet.edges[[j]], max.depth), max.depth)
        neeat_score(w.depth, n.depth, raw.depth, neeat.par)
      }
      sapply(gs.ids, neeat_gs)
    }
  }
  else if (method == "hyper") {
    N <- dim(core.sets)[1]
    M <- colSums(core.sets)
    n <- colSums(gene.sets)
    neeat_cs <- function(i)
    {
      cs <- column(core.sets, i)
      sapply(gs.ids, function(j) neeat_hyper(N, n[j], M[i], sum(cs & column(gene.sets, j)), neeat.par))
    }
  }
  else {
    stop("Incorrect parameters!")
  }
  sapply(cs.ids, neeat_cs)
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

neeat_score_init <- function(w.depth, n.depth, neeat.par)
{
  p <- n.depth / sum(n.depth)
  w <- w.depth * p
  cov <- -tcrossprod(w)
  diag(cov) <- w.depth * w * (1-p)
  
  neeat.par$N <- sum(n.depth)
  neeat.par$avg.score <- sum(w)
  neeat.par$var.score <- sum(cov)
  neeat.par$perm.score <- list()
}

neeat_score <- function(w.depth, n.depth, raw.depth, neeat.par)
{
  n <- sum(raw.depth)

  raw.score <- sum(w.depth * raw.depth)
  avg.score <- neeat.par$avg.score * n
  var.score <- neeat.par$var.score * n

  if (neeat.par$use.multinom) {
    permutation <- rmultinom
  }
  else {
    var.score <- var.score * (neeat.par$N - n) / (neeat.par$N - 1)
    permutation <- rmultihyper
  }

  if (!is.nan(var.score) && var.score > 0)
    z.score <- (raw.score - avg.score) / sqrt(var.score)
  else
    z.score <- 0

  if (z.score >= neeat.par$z.threshold) {
    if (neeat.par$n.perm > 0) {
      if(n > length(neeat.par$perm.score) || is.null(neeat.par$perm.score[[n]]))
        neeat.par$perm.score[[n]] <- colSums(w.depth * permutation(neeat.par$n.perm, n, n.depth))
      p.value <- sum(neeat.par$perm.score[[n]] >= raw.score) / neeat.par$n.perm
    }
    else {
      p.value <- pmultihyper(raw.score, n, n.depth, w.depth)
    }
  }
  else {
    p.value <- 1.0
  }

  if (neeat.par$verbose)
    c(z.score, p.value, raw.score, avg.score, var.score)
  else
    c(z.score, p.value)
}

neeat_net <- function(core.set, gene.set, net, n.gene, neeat.par)
{
  max.depth <- min(2, neeat.par$max.depth)
  w.depth <- c(0, neeat.par$rho^(0:max.depth))
  
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
  
  neeat_score_init(w.depth, n.depth, neeat.par)
  neeat_score(w.depth, n.depth, raw.depth, neeat.par)
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

neeat_hyper <- function(N, n, M, m, neeat.par)
{
  avg.score <- n * M / N
  var.score <- n * (M / N) * ((N - M) / N) * ((N - n) / (N - 1))
  
  if (var.score > 0)
    z.score <- (m - avg.score) / sqrt(var.score)
  else
    z.score <- 0
  
  p.value <- phyper(m-1, M, N-M, n, lower.tail=F)

  if (neeat.par$verbose)
    c(z.score, p.value, m, avg.score, var.score)
  else
    c(z.score, p.value)
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
    all.gene <- gene.set    
  }
  else
    all.gene <- unique(term.table[,1])
  
  toMatrix(term.table, all.gene)
}


#' Turn a two-columns mapping table to matrix
#'
#' Generate the sparse mapping matrix from a two-columns table.
#'
#' This function generates the sparse logical matrix from the two-columns table returned by \code{\link{toTable}} 
#' in Bioconductor package \code{\link{AnnotationDbi}}.
#' 
#' @param table The two-columns table as matrix or data frame. The other columns will not be used, if available.
#' @param rows The row names for the mapping matrix.
#' @param cols The column names for the mapping matrix.
#' 
#' @return This function returns a sparse logical matrix.
#' 
#' @seealso \code{\link{toTable}}
#'
#' @examples
#' 
#' \dontrun{
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("org.Hs.eg.db")
#' library(org.Hs.eg.db)
#' x <- toMatrix(toTable(org.Hs.egGO2ALLEGS))
#' }
#'
#' @import Matrix
#' 
#' @export
toMatrix <- function(table, rows = unique(table[,1]), cols = unique(table[,2]))
{
  rid <- seq_along(rows)
  names(rid) <- rows
  cid <- seq_along(cols)
  names(cid) <- cols
  table <- table[table[,1] %in% rows & table[,2] %in% cols, ]
  sparseMatrix(rid[table[,1]], cid[table[,2]], x = T, dims = c(length(rows), length(cols)), dimnames = list(rows, cols))
}


#' Transform the gene lists into the matrix of gene sets
#'
#' Generate the gene sets matrix from several gene lists.
#'
#' This function generates the \code{gene.sets} matrix required by \code{\link{neeat}} method.
#' The name of each gene set is given by the name of corresponding component of \code{gene.lists},
#' and assigned to the column of \code{gene.sets}.
#' 
#' @param gene.lists A list of vectors of gene ids.
#' @param all.gene A vector of all gene ids.
#' 
#' @return This function returns a sparse matrix as required by \code{\link{neeat}} method.
#' 
#' @seealso \code{\link{neeat}}
#'
#' @import Matrix
#' 
#' @export
get_gene_sets <- function(gene.lists, all.gene)
{
  gene.set <- Matrix(F, nrow=length(all.gene), ncol=length(gene.lists), dimnames=list(all.gene, names(gene.lists)))
  for (i in 1:length(gene.lists)) {
    genes <- as.character(gene.lists[[i]])
    genes <- genes[genes %in% all.gene]
    gene.set[genes, i] <- T
  }
  gene.set
}
