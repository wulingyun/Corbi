
#' @export
netDEG <- function(ratio.dist, expr)
{
  net <- getDiffRatioNet(ratio.dist, expr)
  PvalueNetDEG(net, ratio.dist)
}


#' @export
PvalueNetDEG <- function(adj.matrix, ratio.dist)
{
  n.gene <- nrow(adj.matrix)
  in.degree <- colSums(adj.matrix)
  out.degree <- rowSums(adj.matrix)
  score <- out.degree - in.degree
  up.gene <- score >= 0
  down.gene <- score <= 0
  up.pvalue <- down.pvalue <- rep(1, n.gene)
  up.pvalue[up.gene] <- ratio.dist$rate['up'] * pnbinom(score[up.gene], size = ratio.dist$up['size'], mu = ratio.dist$up['mu'], lower.tail = FALSE)
  down.pvalue[down.gene] <- ratio.dist$rate['down'] * pnbinom(-score[down.gene], size = ratio.dist$down['size'], mu = ratio.dist$down['mu'], lower.tail = FALSE)
  up.pvalue[!up.gene] <- 1 - down.pvalue[!up.gene]
  down.pvalue[!down.gene] <- 1 - up.pvalue[!down.gene]
  pvalue <- 2 * pmin(up.pvalue, down.pvalue)
  return(list(score=score, pvalue=pvalue, up.pvalue=up.pvalue, down.pvalue=down.pvalue))
}


#' @import MASS
#'
#' @export
getRatioDistribution <- function(expr, p.edge = 0.1)
{
  dist <- .Call(ND_RatioDistribution, expr, p.edge)
  net <- lapply(1:dim(expr)[2], function(i) getDiffRatioNet(dist, expr[,i]))
  diff <- as.vector(sapply(net, rowSums) - sapply(net, colSums))
  diff.up <- diff[diff >= 0]
  diff.down <- -diff[diff <= 0]
  dist$up <- fitdistr(diff.up, "negative binomial")$estimate
  dist$down <- fitdistr(diff.down, "negative binomial")$estimate
  dist$rate <- c(up = length(diff.up)/length(diff), down = length(diff.down)/length(diff))
  dist
}


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
