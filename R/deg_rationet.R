
#' @export
netDEG <- function(adj.matrix, p.edge = NULL) {
  n.gene <- nrow(adj.matrix)
  n.edge <- sum(adj.matrix)
  in.degree <- colSums(adj.matrix)
  out.degree <- rowSums(adj.matrix)
  score <- out.degree - in.degree

  if (is.null(p.edge)) p.edge <- n.edge / choose(n.gene, 2)
  pvalue <- .Call(ND_PvalueNetDEG, abs(score), n.gene, p.edge)
  up.pvalue <- ifelse(score > 0, pvalue, 1-pvalue)
  down.pvalue <- ifelse(score < 0, pvalue, 1-pvalue)
  pvalue <- 2.0 * pvalue

  return(list(score=score, pvalue=pvalue, up.pvalue=up.pvalue, down.pvalue=down.pvalue))
}


#' @export
getRatioDistribution <- function(expr, p = 0.1)
  .Call(ND_RatioDistribution, expr, p)


#' @export
getRatioNet <- function(ratio.dist, expr)
  .Call(ND_RatioNet, ratio.dist$LB, ratio.dist$UB, expr)


#' @export
random_net <- function(size, p.edge, sparse = TRUE) {
  max.edge <- choose(size, 2)
  edge.id <- which(runif(max.edge) <= p.edge)
  edge.dir <- runif(length(edge.id)) <= 0.5
  edge.1 <- edge.id[edge.dir]
  edge.2 <- edge.id[!edge.dir]

  rows <- rep.int(2:size, times=1:(size-1))
  cols <- sequence(1:(size-1))
  i <- c(rows[edge.1], cols[edge.2])
  j <- c(cols[edge.1], rows[edge.2])

  if (sparse)
  {
    adj.matrix <- sparseMatrix(i, j, x=TRUE, dims=c(size, size))
  }
  else
  {
    adj.matrix <- matrix(FALSE, size, size)
    adj.matrix[cbind(i, j)] <- TRUE
  }
  adj.matrix
}
