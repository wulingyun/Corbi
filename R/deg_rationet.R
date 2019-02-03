#' netDEG: Differentially expressed gene identification method
#' 
#' @export
netDEG <- function(ref.expr.matrix, expr.matrix, p.edge = 0.1, log.expr = FALSE, use.parallel = FALSE)
{
  if (use.parallel && requireNamespace("BiocParallel"))
  {
    lapply <- BiocParallel::bplapply
  }

  if (!log.expr) {
    ref.expr.matrix <- log(ref.expr.matrix)
    expr.matrix <- log(expr.matrix)
  }
  n <- dim(expr.matrix)[1]
  m <- dim(expr.matrix)[2]
  ref.ratio.dist <- get_ratio_distribution(ref.expr.matrix, p.edge, log.expr = T)
  p <- lapply(1:m, function(i) netDEG_pvalue(ref.ratio.dist, expr.matrix[,i], log.expr = T))

  up <- sapply(p, function(i) i$up)
  down <- sapply(p, function(i) i$down)
  twoside <- sapply(p, function(i) i$twoside)
  dimnames(up) <- dimnames(down) <- dimnames(twoside) <- dimnames(expr.matrix)

  g.up <- sapply(1:n, function(i) p_combine(up[i,])$p)
  g.down <- sapply(1:n, function(i) p_combine(down[i,])$p)
  g.twoside <- sapply(1:n, function(i) p_combine(twoside[i,])$p)
  names(g.up) <- names(g.down) <- names(g.twoside) <- rownames(expr.matrix)

  s.up <- sapply(1:m, function(i) p_combine(up[,i])$p)
  s.down <- sapply(1:m, function(i) p_combine(down[,i])$p)
  s.twoside <- sapply(1:m, function(i) p_combine(twoside[,i])$p)
  names(s.up) <- names(s.down) <- names(s.twoside) <- colnames(expr.matrix)
  
  zero.genes <- rowSums(is.finite(ref.expr.matrix)) == 0 | rowSums(is.finite(expr.matrix)) == 0
  n0 <- sum(zero.genes)
  m1 <- ref.expr.matrix[zero.genes, , drop = F]
  m1 <- ifelse(is.finite(m1), exp(m1), 0)
  m2 <- expr.matrix[zero.genes, , drop = F]
  m2 <- ifelse(is.finite(m2), exp(m2), 0)
  g.up[zero.genes] <- sapply(1:n0, function(i) p_combine(rep(stats::t.test(m1[i, ], m2[i, ], alternative = "less")$p.value, m))$p)
  g.down[zero.genes] <- sapply(1:n0, function(i) p_combine(rep(stats::t.test(m1[i, ], m2[i, ], alternative = "greater")$p.value, m))$p)
  g.twoside[zero.genes] <- sapply(1:n0, function(i) p_combine(rep(stats::t.test(m1[i, ], m2[i, ], alternative = "two.sided")$p.value, m))$p)

  return(list(up = up, down = down, twoside = twoside,
              gene = list(up = g.up, down = g.down, twoside = g.twoside),
              sample = list(up = s.up, down = s.down, twoside = s.twoside)))
}

#' Calculate netDEG statistics and p-values
#'
#' @export
netDEG_pvalue <- function(ref.ratio.dist, expr.val, log.expr = FALSE)
{
  if (!log.expr) expr.val <- log(expr.val)
  score <- get_diff_ratio_net(ref.ratio.dist, expr.val, log.expr = T)$diff
  pvalue <- stats::pnbinom(abs(score), size = ref.ratio.dist$NB['size'], mu = ref.ratio.dist$NB['mu'], lower.tail = FALSE)
  p = pvalue * 0.5
  up = ifelse(score > 0, p, 1-p)
  down = ifelse(score < 0, p, 1-p)
  return(list(up = up, down = down, twoside = pvalue))
}


#' Calculate expression ratio distribution
#' 
#' @export
get_ratio_distribution <- function(ref.expr.matrix, p.edge = 0.1, log.expr = FALSE, use.parallel = FALSE)
{
  if (use.parallel && requireNamespace("BiocParallel"))
  {
    lapply <- BiocParallel::bplapply
  }

  if (!log.expr) ref.expr.matrix <- log(ref.expr.matrix)
  dist <- .Call(ND_RatioDistribution, ref.expr.matrix, p.edge)
  diff <- unlist(lapply(1:dim(ref.expr.matrix)[2], function(i) get_diff_ratio_net(dist, ref.expr.matrix[,i], log.expr = T)$diff))
  diff <- diff[!is.na(diff)]
  dist$NB <- MASS::fitdistr(abs(diff), "negative binomial", lower = c(1e-10, 1e-10))$estimate
  dist
}

#' Calculate expression ratio distribution
#' 
#' @export
get_ratio_distribution2 <- function(ref.expr.matrix, p.edge = 0.1, p.trim = 0.3, log.expr = FALSE, use.parallel = FALSE)
{
  if (use.parallel && requireNamespace("BiocParallel"))
  {
    lapply <- BiocParallel::bplapply
  }
  
  if (!log.expr) ref.expr.matrix <- log(ref.expr.matrix)
  dist <- .Call(ND_RatioDistribution2, ref.expr.matrix, p.edge, p.trim)
  diff <- unlist(lapply(1:dim(ref.expr.matrix)[2], function(i) get_diff_ratio_net(dist, ref.expr.matrix[,i], log.expr = T)$diff))
  diff <- diff[!is.na(diff)]
  dist$NB <- MASS::fitdistr(abs(diff), "negative binomial", lower = c(1e-10, 1e-10))$estimate
  dist
}

#' Construct differential expression ratio network
#' 
#' @export
get_diff_ratio_net <- function(ref.ratio.dist, expr.val, log.expr = FALSE)
{
  if (!log.expr) expr.val <- log(expr.val)
  edges <- .Call(ND_DiffRatioNet, ref.ratio.dist$LB, expr.val)
  n <- length(expr.val)
  net <- sparseMatrix(i = edges$i, j = edges$j, dims = c(n, n))
  d <- get_adjusted_deg_diff(net, expr.val)
  list(net = net, diff = d$diff, degree = d$degree)
}


#' Calculate adjusted degree differences for given network
#'
get_adjusted_deg_diff <- function(net, log.expr.val, p = 0.5)
{
  g.NA <- !is.finite(log.expr.val)
  d.out <- rowSums(net)
  d.in <- colSums(net)
  d.out[g.NA] <- NA
  d.in[g.NA] <- NA
  d.sum <- d.out + d.in
  d.diff <- d.out - d.in
  adj.diff <- d.diff - ceiling(stats::median(d.diff[d.sum <= stats::quantile(d.sum, p, na.rm = T)], na.rm = T))
  list(diff = adj.diff, degree = list(diff = d.diff, sum = d.sum))
}


#' Calculate combined p-value by Fisher's method
#' 
#' @export
p_combine <- function(p)
{
  p <- p[!is.na(p)]
  p[p > 1] <- 1
  p[p < .Machine$double.xmin] <- .Machine$double.xmin
  chisq <- (-2) * sum(log(p))
  df <- 2 * length(p)
  list(chisq = chisq, df = df, p = stats::pchisq(chisq, df, lower.tail = FALSE))
}
