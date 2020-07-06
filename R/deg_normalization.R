#' Calculate expression ratio variances
#' 
#' Calculate the variances of expression ratios for each pair of genes.
#' 
#' @param expr.matrix The expression matrix. Each row represents a gene and each column represents a sample.
#' @param log.expr Logical variable indicating whether the input expression matrix is in logarithmic scale.
#' @param use.parallel Logical variable indicating to use the BiocParallel package to accelerate computation.
#' 
#' @return This function will return a numeric matrix with element [i,j] represents the variance of 
#' expression ratios for gene pairs (i, j).
#'   
#' @export
get_ratio_variance <- function(expr.matrix, log.expr = FALSE, use.parallel = FALSE)
{
  if (use.parallel && requireNamespace("BiocParallel")) {
    lapply <- BiocParallel::bplapply
  }

  if (!log.expr) expr.matrix <- log(expr.matrix)
  if (use.parallel) {
    n.genes <- dim(expr.matrix)[1]
    V <- lapply(1:ceiling((n.genes-1)/2), function (i) .Call(ND_RatioVarianceParI, expr.matrix, i))
    .Call(ND_ParMerge, V, n.genes, Inf, TRUE)
  }
  else {
    .Call(ND_RatioVariance, expr.matrix)
  }
}

#' Identify the stable genes
#' 
#' Identify the stable genes by using URG (uniform ratio graph) method.
#' 
#' @param expr.matrix The expression matrix. Each row represents a gene and each column represents a sample.
#' @param max.stable The maximal number of genes that are selected as the stable genes.
#' @param p.edge The percentage of gene pairs that are selected into the uniform ratio graph.
#' @param log.expr Logical variable indicating whether the input expression matrix is in logarithmic scale.
#' @param use.parallel Logical variable indicating to use the BiocParallel package to accelerate computation.
#' 
#' @return This function will return an integer vector containing the indexes of stable genes.
#' 
#' @export
get_stable_genes <- function(expr.matrix, max.stable, p.edge = 0.25, log.expr = FALSE, use.parallel = FALSE)
{
  ratio_var <- get_ratio_variance(expr.matrix, log.expr, use.parallel)
  cutoff <- stats::quantile(ratio_var[lower.tri(ratio_var)], p.edge)
  adj_matrix <- ifelse(ratio_var <= cutoff, 1, 0)

  g <- igraph::graph_from_adjacency_matrix(adj_matrix)
  components <- igraph::components(g)
  max_component <- which(components$membership == which.max(components$csize))
  node_degree <- igraph::degree(g, max_component)
  
  n_stable <- min(length(max_component), max.stable)
  cutoff <- sort(node_degree, decreasing = TRUE)[n_stable]
  max_component[node_degree >= cutoff]
}

#' Calculate normalization factors for URG method
#' 
#' Calculate the normalization factor for each sample by using URG (uniform ratio graph) method.
#' 
#' @param expr.matrix The expression matrix. Each row represents a gene and each column represents a sample.
#' @param p.edge The percentage of gene pairs that are selected into the uniform ratio graph.
#' @param p.gene The maximal percentage of genes that are selected as the stable genes.
#' @param log.expr Logical variable indicating whether the input expression matrix is in logarithmic scale.
#' @param zero.as.dropout Logical variable indicating whether the zero expressions are regarded as dropouts.
#' @param use.parallel Logical variable indicating to use the BiocParallel package to accelerate computation.
#' 
#' @return This function will return a numeric vector with each element [i] represents the normalization
#' factor of sample (i).
#' 
#' @seealso \code{\link{URG_normalize}}
#' 
#' @references Xinhan Ye, Ling-Yun Wu. URG: a new normalization method for 
#' gene expression data based on graph model. Manuscript.
#'   
#' @export
URG_getFactor <- function(expr.matrix, p.edge = 0.25, p.gene = 0.4, log.expr = FALSE, zero.as.dropout = FALSE, use.parallel = FALSE)
{
  if (!log.expr) expr.matrix <- log(expr.matrix)
  if (!zero.as.dropout) {
    zero.expr = min(expr.matrix[is.finite(expr.matrix)], 0, na.rm = TRUE) - log(2)
    expr.matrix[!is.finite(expr.matrix)] <- zero.expr
  }

  stable_genes <- get_stable_genes(expr.matrix, round(dim(expr.matrix)[1] * p.gene), p.edge, TRUE, use.parallel)
  stable_expr <- expr.matrix[stable_genes, , drop = FALSE]
  stable_expr[!is.finite(stable_expr)] <- NA
  geom <- exp(apply(stable_expr, 2, mean, na.rm = TRUE))
  mean(geom) / geom
}

#' Normalize using given factors
#' 
#' Normalize the expression matrix by using the given factor for each sample.
#' 
#' @param expr.matrix The expression matrix. Each row represents a gene and each column represents a sample.
#' @param factor The numeric vector of normalization factors.
#' @param log.expr Logical variable indicating whether the input expression matrix is in logarithmic scale.
#' 
#' @return This function will return a numeric matrix with the same dimension of \code{expr.matrix}.
#' 
#' @seealso \code{\link{URG_getFactor}}
#' 
#' @export
URG_normalize <- function(expr.matrix, factor, log.expr = FALSE)
{
  n <- nrow(expr.matrix)
  if (log.expr) {
    norm.matrix <- expr.matrix + rep(log(factor), each = n)
  }
  else {
    norm.matrix <- expr.matrix * rep(factor, each = n)
  }
  norm.matrix
}
