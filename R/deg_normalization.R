


#' @export
get_ratio_variance <- function(expr.matrix, log.expr = FALSE)
{
  if (!log.expr) expr.matrix <- log(expr.matrix)
  .Call(ND_RatioVariance, expr.matrix)
}

#' @export
URG_getFactor <- function(expr.matrix, p_var = 0.25, p_gene = 0.4, log.expr = FALSE)
{
  if (!log.expr) {
    expr.matrix[expr.matrix <= 0] <- 0.01
    expr.matrix <- log(expr.matrix)
    log.expr <- TRUE
  }
  ratio_var <- get_ratio_variance(expr.matrix, log.expr = log.expr)
  cutoff <- quantile(ratio_var[lower.tri(ratio_var)], p_var)
  adj_matrix <- ifelse(ratio_var <= cutoff, 1, 0) 

  g <- igraph::graph_from_adjacency_matrix(adj_matrix)
  max_component <- which(igraph::components(g)$membership == 1)
  node_degree <- igraph::degree(g, max_component)
  
  n_stable <- min(length(max_component), round(dim(adj_matrix)[1] * p_gene))
  cutoff <- sort(node_degree, decreasing = TRUE)[n_stable]
  stable_genes <- max_component[node_degree >= cutoff]

  stable_expr <- expr.matrix[stable_genes, , drop = FALSE]

  geom <- exp(apply(stable_expr, 2, mean))
  mean(geom) / geom
}

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
