
#' @export
netDEG <- function(ratio.dist, expr)
{
  net <- getDiffRatioNet(ratio.dist, expr)
  PvalueNetDEG(net, ratio.dist)
}


#' @export
PvalueNetDEG <- function(adj.matrix, ratio.dist)
{
  score <- getAdjustedDiff(adj.matrix)
  pvalue <- pnbinom(abs(score), size = ratio.dist$NB['size'], mu = ratio.dist$NB['mu'], lower.tail = FALSE)
  p = pvalue * 0.5
  up = ifelse(score > 0, p, 1-p)
  down = ifelse(score < 0, p, 1-p)

  return(list(up = up, down = down, pvalue = pvalue))
}


#' @import MASS
#'
#' @export
getRatioDistribution <- function(expr, p.edge = 0.1)
{
  dist <- .Call(ND_RatioDistribution, expr, p.edge)
  diff <- as.vector(sapply(1:dim(expr)[2], function(i) getAdjustedDiff(getDiffRatioNet(dist, expr[,i]))))
  dist$NB <- fitdistr(abs(diff), "negative binomial", lower = c(1e-10, 1e-10))$estimate
  dist
}


#' @export
getDiffRatioNet <- function(ratio.dist, expr)
{
  .Call(ND_DiffRatioNet, ratio.dist$LB, ratio.dist$UB, expr)
}


#' @export
getAdjustedDiff <- function(net, p = 0.5)
{
  n.gene <- dim(net)[1]
  d.out <- rowSums(net)
  d.in <- colSums(net)
  d.sum <- d.out + d.in
  d.diff <- d.out - d.in
  d.diff - ceiling(median(d.diff[d.sum <= quantile(d.sum, p)]))
}


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
