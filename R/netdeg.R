
#' @export
netDEG <- function(adj.matrix) {
  n.gene <- nrow(adj.matrix)
  n.edge <- sum(adj.matrix)
  in.degree <- colSums(adj.matrix)
  out.degree <- rowSums(adj.matrix)
  score <- out.degree - in.degree
  up.score <- ifelse(score > 0, score, 0)
  down.score <- ifelse(score < 0, -score, 0)
  net.degree <- abs(score)
  total.degree <- in.degree + out.degree
  
  p.edge <- n.edge/choose(n.gene, 2)
  pvalue <- .Call(ND_PvalueNetDEG, net.degree, n.gene, p.edge)

  return(list(pvalue=pvalue, score=net.degree, up.score=up.score, down.score=down.score))
}


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
