#' Calculate expression ratio variances
#' 
#' Calculate the variances of expression ratios for each pair of genes.
#' 
#' @param expr.matrix The expression matrix. Each row represents a gene and each column represents a sample.
#' @param log.expr Logical variable indicating whether the input expression matrix is in logarithmic scale.
#' 
#' @return This function will return a numeric matrix with element [i,j] represents the variance of 
#' expressioin ratios for gene pairs (i, j).
#'   
#' @export
get_ratio_variance <- function(expr.matrix, log.expr = FALSE)
{
  if (!log.expr) expr.matrix <- log(expr.matrix)
  .Call(ND_RatioVariance, expr.matrix)
}

#' Calculate normalization factors for URG method
#' 
#' Calculate the normalization factor for each sample by using URG (uniform ratio graph) method.
#' 
#' @param expr.matrix The expression matrix. Each row represents a gene and each column represents a sample.
#' @param p.edge The percentage of gene pairs that are selected into the uniform ratio graph.
#' @param p.gene The maximal percentage of genes that are selected as the stable genes.
#' @param log.expr Logical variable indicating whether the input expression matrix is in logarithmic scale.
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
URG_getFactor <- function(expr.matrix, p.edge = 0.25, p.gene = 0.4, log.expr = FALSE)
{
  if (!log.expr) {
    expr.matrix[expr.matrix <= 0] <- 0.01
    expr.matrix <- log(expr.matrix)
    log.expr <- TRUE
  }
  ratio_var <- get_ratio_variance(expr.matrix, log.expr = log.expr)
  cutoff <- stats::quantile(ratio_var[lower.tri(ratio_var)], p.edge)
  adj_matrix <- ifelse(ratio_var <= cutoff, 1, 0) 

  g <- igraph::graph_from_adjacency_matrix(adj_matrix)
  max_component <- which(igraph::components(g)$membership == 1)
  node_degree <- igraph::degree(g, max_component)
  
  n_stable <- min(length(max_component), round(dim(adj_matrix)[1] * p.gene))
  cutoff <- sort(node_degree, decreasing = TRUE)[n_stable]
  stable_genes <- max_component[node_degree >= cutoff]

  stable_expr <- expr.matrix[stable_genes, , drop = FALSE]

  geom <- exp(apply(stable_expr, 2, mean))
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
