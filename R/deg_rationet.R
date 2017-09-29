
#' @export
netDEG <- function(ratio.dist, expr)
{
  net <- getDiffRatioNet(ratio.dist, expr)
  PvalueNetDEG(net, ratio.dist$p.edge)
}


#' @export
PvalueNetDEG <- function(adj.matrix, p.edge = NULL)
{
  n.gene <- nrow(adj.matrix)
  in.degree <- colSums(adj.matrix)
  out.degree <- rowSums(adj.matrix)
  score <- out.degree - in.degree

  if (is.null(p.edge)) p.edge <- sum(adj.matrix) / choose(n.gene, 2)

  coefs <- list(alpha=c(0.5049524, 0.5049657), norm2=c(1.0, 200.0, 16.663143, 1.110288))
  z <- poly(p.edge, degree = 2, coefs = coefs)
  sd <- (0.4087 + 2.0832*z[,1] - 0.5833*z[,2]) * n.gene
  up.pvalue <- pnorm(-score, mean = 0, sd = sd, lower.tail = TRUE)
  down.pvalue <- pnorm(score, mean = 0, sd = sd, lower.tail = TRUE)
  pvalue <- pmin(up.pvalue, down.pvalue) * 2

  return(list(score=score, pvalue=pvalue, up.pvalue=up.pvalue, down.pvalue=down.pvalue))
}


#' @export
getRatioDistribution <- function(expr, p.edge = 0.1)
  .Call(ND_RatioDistribution, expr, p.edge)


#' @export
getDiffRatioNet <- function(ratio.dist, expr)
  .Call(ND_DiffRatioNet, ratio.dist$LB, ratio.dist$UB, expr)


#' @export
random_net <- function(size, p.edge, sparse = TRUE)
{
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
