
#' @export
netDEG <- function(ref.ratio.dist, expr)
{
  PvalueNetDEG(ref.ratio.dist, getDiffRatioNet(ref.ratio.dist, expr))
}


#' @export
PvalueNetDEG <- function(ref.ratio.dist, net)
{
  score <- .netDEG.getAdjustedDiff(net)
  pvalue <- pnbinom(abs(score), size = ref.ratio.dist$NB['size'], mu = ref.ratio.dist$NB['mu'], lower.tail = FALSE)
  p = pvalue * 0.5
  up = ifelse(score > 0, p, 1-p)
  down = ifelse(score < 0, p, 1-p)

  return(list(up = up, down = down, pvalue = pvalue))
}


#'
#'
#' @export
getRatioDistribution <- function(ref.expr.matrix, p.edge = 0.1)
{
  dist <- .Call(ND_RatioDistribution, ref.expr.matrix, p.edge)
  diff <- as.vector(sapply(1:dim(ref.expr.matrix)[2], function(i) .netDEG.getAdjustedDiff(getDiffRatioNet(dist, ref.expr.matrix[,i]))))
  dist$NB <- MASS::fitdistr(abs(diff), "negative binomial", lower = c(1e-10, 1e-10))$estimate
  dist
}


#' @export
getDiffRatioNet <- function(ref.ratio.dist, expr)
{
  .Call(ND_DiffRatioNet, ref.ratio.dist$LB, ref.ratio.dist$UB, expr)
}


.netDEG.getAdjustedDiff <- function(net, p = 0.5)
{
  n.gene <- dim(net)[1]
  d.out <- rowSums(net)
  d.in <- colSums(net)
  d.sum <- d.out + d.in
  d.diff <- d.out - d.in
  d.diff - ceiling(median(d.diff[d.sum <= quantile(d.sum, p)]))
}
