#' netDEG: Differentially expressed gene identification method
#' 
#' Perform netDEG for two group samples.
#' 
#' @param ref.expr.matrix The reference expression matrix. Each row represents a gene and each column represents a sample.
#' @param expr.matrix The test expression matrix. Each row represents a gene and each column represents a sample.
#' @param p.edge The expected probability of edges in the expression ratio network for a normal sample.
#' @param summarize Character vector indicating how to summarize the results. Available methods are \code{c("gene", "sample")}.
#' @param log.expr Logical variable indicating whether the input expression matrix is in logarithmic scale.
#' @param use.parallel Logical variable indicating to use the BiocParallel package to accelerate computation.
#' 
#' @return This function will return a list with the following components:
#'   \item{up}{A numeric matrix with same dimension as \code{expr.matrix}, containing the p-values of up-regulation test.}
#'   \item{down}{A numeric matrix with same dimension as \code{expr.matrix}, containing the p-values of down-regulation test.}
#'   \item{twoside}{A numeric matrix with same dimension as \code{expr.matrix}, containing the p-values of twoside test.}
#'   \item{gene}{A list containing the gene-wise summaried results, containing three components: \code{up}, \code{down}, 
#'   and \code{twoside}. Available if the corresponding method is specified in \code{summarize} argument.}
#'   \item{sample}{A list containing the sample-wise summaried results, containing three components: \code{up}, \code{down},
#'   and \code{twoside}. Available if the corresponding method is specified in \code{summarize} argument.}
#' 
#' @export
netDEG <- function(ref.expr.matrix, expr.matrix, p.edge = 0.1, summarize = c("gene", "sample"),
                   log.expr = FALSE, use.parallel = FALSE)
{
  if (use.parallel && requireNamespace("BiocParallel"))
  {
    lapply <- BiocParallel::bplapply
  }

  if (!log.expr) {
    ref.expr.matrix <- log(ref.expr.matrix)
    expr.matrix <- log(expr.matrix)
  }
  n.samples <- dim(expr.matrix)[2]
  dist <- get_ratio_distribution(ref.expr.matrix, p.edge, log.expr = T)
  p <- lapply(1:n.samples, function(i) netDEG_pvalue(dist, expr.matrix[,i], log.expr = T))
  rm(dist)
  
  up <- sapply(p, function(i) i$up)
  down <- sapply(p, function(i) i$down)
  twoside <- sapply(p, function(i) i$twoside)
  dimnames(up) <- dimnames(down) <- dimnames(twoside) <- dimnames(expr.matrix)

  results <- list(up = up, down = down, twoside = twoside)
  
  if ("gene" %in% summarize)
  {
    n.genes <- dim(expr.matrix)[1]
    n.refs <- dim(ref.expr.matrix)[2]
    dist <- get_ratio_distribution(expr.matrix, p.edge, log.expr = T)
    p <- lapply(1:n.refs, function(i) netDEG_pvalue(dist, ref.expr.matrix[,i], log.expr = T))
    rm(dist)
    up <- cbind(up, sapply(p, function(i) i$up))
    down <- cbind(down, sapply(p, function(i) i$down))
    twoside <- cbind(twoside, sapply(p, function(i) i$twoside))

    g.up <- sapply(1:n.genes, function(i) p_combine(up[i,])$p)
    g.down <- sapply(1:n.genes, function(i) p_combine(down[i,])$p)
    g.twoside <- sapply(1:n.genes, function(i) p_combine(twoside[i,])$p)
    names(g.up) <- names(g.down) <- names(g.twoside) <- rownames(expr.matrix)
    
    zero.genes <- rowSums(is.finite(ref.expr.matrix)) == 0 | rowSums(is.finite(expr.matrix)) == 0
    n1 <- sum(zero.genes)
    n2 <- n.samples + n.refs
    if (n1 > 0)
    {
      m1 <- ref.expr.matrix[zero.genes, , drop = F]
      m1 <- ifelse(is.finite(m1), exp(m1), 0)
      m2 <- expr.matrix[zero.genes, , drop = F]
      m2 <- ifelse(is.finite(m2), exp(m2), 0)
      g.up[zero.genes] <- sapply(1:n1, function(i) p_combine(rep(stats::wilcox.test(m1[i, ], m2[i, ], exact = F, alternative = "less")$p.value, n2))$p)
      g.down[zero.genes] <- sapply(1:n1, function(i) p_combine(rep(stats::wilcox.test(m1[i, ], m2[i, ], exact = F, alternative = "greater")$p.value, n2))$p)
      g.twoside[zero.genes] <- sapply(1:n1, function(i) p_combine(rep(stats::wilcox.test(m1[i, ], m2[i, ], exact = F, alternative = "two.sided")$p.value, n2))$p)
    }
    
    results$gene <- list(up = g.up, down = g.down, twoside = g.twoside)
  }

  if ("sample" %in% summarize)
  {
    s.up <- sapply(1:n.samples, function(i) p_combine(up[,i])$p)
    s.down <- sapply(1:n.samples, function(i) p_combine(down[,i])$p)
    s.twoside <- sapply(1:n.samples, function(i) p_combine(twoside[,i])$p)
    names(s.up) <- names(s.down) <- names(s.twoside) <- colnames(expr.matrix)
    
    results$sample <- list(up = s.up, down = s.down, twoside = s.twoside)
  }

  return(results)
}

#' Calculate netDEG p-values
#' 
#' Perform the single or two side tests and calculate the p-values.
#' 
#' @param ref.ratio.dist The expression ratio distribution profile returned by \code{get_ratio_distribution} or \code{get_ratio_distribution2}.
#' @param expr.val Numeric vector of gene expression values in the sample.
#' @param log.expr Logical variable indicating whether the input expression vector is in logarithmic scale.
#'
#' @return This function will return a list with the following components:
#'   \item{up}{A numeric vector containing the p-values of up-regulation test.}
#'   \item{down}{A numeric vector containing the p-values of down-regulation test.}
#'   \item{twoside}{A numeric vector containing the p-values of twoside test.}
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
#' Calculate the lower and upper quantiles of expression ratios for each pair of genes,
#' and estimate the parameters of negative binomial distribution from reference expression data.
#' 
#' @param ref.expr.matrix The reference expression matrix. Each row represents a gene and each column represents a sample.
#' @param p.edge The expected probability of edges in the expression ratio network for a normal sample.
#' @param log.expr Logical variable indicating whether the input expression matrix is in logarithmic scale.
#' @param use.parallel Logical variable indicating to use the BiocParallel package to accelerate computation.
#' 
#' @return This function will return a list with the following components:
#'   \item{LB}{A numeric matrix with element [i,j] represents the lower quantile of expressioin ratios for gene pairs (i, j).}
#'   \item{NB}{A numeric vector with two elements: \code{size} and \code{mu}, 
#'   which are the estimated parameters of negative binomial distribution.}
#'   \item{p.edge}{The used input parameter \code{p.edge}.}
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
#' Calculate the lower and upper quantiles of expression ratios after trimming the extreme values, 
#' and estimate the parameters of negative binomial distribution from reference expression data.
#' 
#' @param ref.expr.matrix The reference expression matrix. Each row represents a gene and each column represents a sample.
#' @param p.edge The expected probability of edges in the expression ratio network for a normal sample.
#' @param p.trim The percentage of lower or upper extreme values to be trimmed from the expression ratios for each pair of genes.
#' @param log.expr Logical variable indicating whether the input expression matrix is in logarithmic scale.
#' @param use.parallel Logical variable indicating to use the BiocParallel package to accelerate computation.
#' 
#' @return This function will return a list with the following components:
#'   \item{LB}{A numeric matrix with element [i,j] represents the lower quantile of trimmed expressioin ratios for gene pairs (i, j).}
#'   \item{NB}{A numeric vector with two elements: \code{size} and \code{mu}, 
#'   which are the estimated parameters of negative binomial distribution.}
#'   \item{p.edge}{The used input parameter \code{p.edge}.}
#'   \item{p.trim}{The used input parameter \code{p.trim}.}
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
#' Construct the differential expression ratio network for a single sample.
#' 
#' @param ref.ratio.dist The expression ratio distribution profile returned by \code{get_ratio_distribution} or \code{get_ratio_distribution2}.
#' @param expr.val Numeric vector of gene expression values in the sample.
#' @param log.expr Logical variable indicating whether the input expression vector is in logarithmic scale.
#' 
#' @return This function will return a list with the following components:
#'   \item{net}{The binary adjacent matrix of differential expression ratio network.}
#'   \item{diff}{A numeric vector containing the adjusted degree differences of all genes.}
#'   \item{degree}{A list containing the raw degree differences and sums of all genes.}
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
#' Calculate the adjusted degree differences for all genes in the given network.
#' 
#' @param net The binary adjacent matrix of differential expression ratio network.
#' @param log.expr.val Numeric vector containing the logarithmic scale gene expression values.
#' @param p The parameter for calculating the adjusted degree differences.
#' 
#' @return This function will return a list with the following components:
#'   \item{diff}{A numeric vector containing the adjusted degree differences of all genes.}
#'   \item{degree}{A list containing the raw degree differences and sums of all genes.}
#' 
#' @export
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
#' The Fisher's method is used to combine the results from several independent tests.
#' 
#' @param p the numeric vector containing the p-values need to combine.
#' 
#' @return This function will return a list with the following components:
#'   \item{chisq}{The chi-squared statistic.}
#'   \item{df}{The degrees of freedom of chi-squared distribution.}
#'   \item{p}{The p-value of chi-squared statistic.}
#' 
#' @export
p_combine <- function(p)
{
  p <- p[!is.na(p)]
  p[p > 1] <- 1
  p[p < .Machine$double.xmin] <- .Machine$double.xmin
  chisq <- (-2) * sum(log(p))
  df <- 2 * length(p)
  p <- stats::pchisq(chisq, df, lower.tail = FALSE)
  if (p == 0) p <- .Machine$double.xmin / chisq
  list(chisq = chisq, df = df, p = p)
}
