#' netDEG: Differentially expressed gene identification method
#' 
#' Perform netDEG for two group samples.
#' 
#' @param ref.expr.matrix The reference expression matrix. Each row represents a gene and each column represents a sample.
#' @param expr.matrix The test expression matrix. Each row represents a gene and each column represents a sample.
#' @param p.edge The expected probability of edges in the expression ratio network for a normal sample.
#' @param summarize Character vector indicating how to summarize the results. Available methods are \code{c("gene", "sample")}.
#' @param summarize.method Character vector indicating the methods used to summarize the results. See \code{p_combine}.
#' @param summarize.shrink Numeric vector indicating the shrink parameter to summarize the results. See \code{p_combine}.
#' @param log.expr Logical variable indicating whether the input expression matrix is in logarithmic scale.
#' @param zero.as.dropout Logical variable indicating whether the zero expressions are regarded as dropouts.
#' @param scale.degree Logical variable indicating whether the degree values are scaled according to the dropout rate.
#' @param use.parallel Logical variable indicating to use the BiocParallel package to accelerate computation.
#' 
#' @return This function will return a list with the following components:
#'   \item{up}{A numeric matrix with same dimension as \code{expr.matrix}, containing the p-values of up-regulation test.}
#'   \item{down}{A numeric matrix with same dimension as \code{expr.matrix}, containing the p-values of down-regulation test.}
#'   \item{twoside}{A numeric matrix with same dimension as \code{expr.matrix}, containing the p-values of twoside test.}
#'   \item{rev}{A list containing the reverse comparison results, containing three components: \code{up}, \code{down}, 
#'   and \code{twoside}. Available if the gene method is specified in \code{summarize} argument.}
#'   \item{gene}{A list containing the gene-wise summaried results, containing three components: \code{up}, \code{down}, 
#'   and \code{twoside}. Available if the gene method is specified in \code{summarize} argument.}
#'   \item{sample}{A list containing the sample-wise summaried results, containing three components: \code{up}, \code{down},
#'   and \code{twoside}. Available if the sample method is specified in \code{summarize} argument.}
#' 
#' @export
netDEG <- function(ref.expr.matrix, expr.matrix, p.edge = 0.1,
                   summarize = c("gene", "sample"), summarize.method = c("sumlog", "sumlog"), summarize.shrink = c(Inf, Inf),
                   log.expr = FALSE, zero.as.dropout = TRUE, scale.degree = TRUE, use.parallel = FALSE)
{
  if (use.parallel && requireNamespace("BiocParallel")) {
    lapply <- BiocParallel::bplapply
  }

  if (!log.expr) {
    ref.expr.matrix <- log(ref.expr.matrix)
    expr.matrix <- log(expr.matrix)
  }
  if (!zero.as.dropout) {
    zero.expr = min(ref.expr.matrix[ref.expr.matrix != -Inf], 0, na.rm = TRUE) - log(10)
    ref.expr.matrix[ref.expr.matrix == -Inf] <- zero.expr
    expr.matrix[expr.matrix == -Inf] <- zero.expr
  }
  n.samples <- dim(expr.matrix)[2]
  
  message("Estimating the ratio distribution from reference samples")
  dist <- get_ratio_distribution(ref.expr.matrix, p.edge, log.expr = TRUE, scale.degree = scale.degree, use.parallel = use.parallel)

  message("Calculating the sample-specific p-values for test samples")
  p <- lapply(1:n.samples, function(i) netDEG_pvalue(dist, expr.matrix[,i], log.expr = TRUE, scale.degree = scale.degree))
  rm(dist)
  
  up <- sapply(p, function(i) i$up)
  down <- sapply(p, function(i) i$down)
  twoside <- sapply(p, function(i) i$twoside)
  dimnames(up) <- dimnames(down) <- dimnames(twoside) <- dimnames(expr.matrix)
  rm(p)

  results <- list(up = up, down = down, twoside = twoside)
  
  if ("gene" %in% summarize) {
    n.genes <- dim(expr.matrix)[1]
    n.refs <- dim(ref.expr.matrix)[2]

    message("Estimating the ratio distribution from test samples")
    dist <- get_ratio_distribution(expr.matrix, p.edge, log.expr = TRUE, scale.degree = scale.degree, use.parallel = use.parallel)

    message("Calculating the sample-specific p-values for reference samples")
    p <- lapply(1:n.refs, function(i) netDEG_pvalue(dist, ref.expr.matrix[,i], log.expr = TRUE, scale.degree = scale.degree))
    rm(dist)

    rev.up <- sapply(p, function(i) i$up)
    rev.down <- sapply(p, function(i) i$down)
    rev.twoside <- sapply(p, function(i) i$twoside)
    dimnames(rev.up) <- dimnames(rev.down) <- dimnames(rev.twoside) <- dimnames(ref.expr.matrix)
    rm(p)

    results$rev <- list(up = rev.up, down = rev.down, twoside = rev.twoside)

    message("Gene-wise summarization")
    
    method <- summarize.method[summarize == "gene"][1]
    if (is.na(method)) method <- "sumlog"
    shrink <- summarize.shrink[summarize == "gene"][1]
    if (is.na(shrink)) shrink <- Inf
    
    up <- cbind(up, rev.down)
    down <- cbind(down, rev.up)
    twoside <- cbind(twoside, rev.twoside)

    g.up <- sapply(1:n.genes, function(i) p_combine(up[i,], method, shrink)$p)
    g.down <- sapply(1:n.genes, function(i) p_combine(down[i,], method, shrink)$p)
    g.twoside <- sapply(1:n.genes, function(i) p_combine(twoside[i,], method, shrink)$p)
    names(g.up) <- names(g.down) <- names(g.twoside) <- rownames(expr.matrix)
    
    zero.genes <- rowSums(is.finite(ref.expr.matrix)) == 0 | rowSums(is.finite(expr.matrix)) == 0
    n1 <- sum(zero.genes)
    n2 <- n.samples + n.refs
    if (n1 > 0) {
      m1 <- ref.expr.matrix[zero.genes, , drop = FALSE]
      m1 <- ifelse(is.finite(m1), exp(m1), 0)
      m2 <- expr.matrix[zero.genes, , drop = FALSE]
      m2 <- ifelse(is.finite(m2), exp(m2), 0)
      g.up[zero.genes] <- sapply(1:n1, function(i) p_combine(rep(stats::wilcox.test(m1[i, ], m2[i, ], exact = FALSE, alternative = "less")$p.value, n2), method, shrink)$p)
      g.down[zero.genes] <- sapply(1:n1, function(i) p_combine(rep(stats::wilcox.test(m1[i, ], m2[i, ], exact = FALSE, alternative = "greater")$p.value, n2), method, shrink)$p)
      g.twoside[zero.genes] <- sapply(1:n1, function(i) p_combine(rep(stats::wilcox.test(m1[i, ], m2[i, ], exact = FALSE, alternative = "two.sided")$p.value, n2), method, shrink)$p)
    }
    
    results$gene <- list(up = g.up, down = g.down, twoside = g.twoside)
  }

  if ("sample" %in% summarize) {
    message("Sample-wise summarization")
    
    method <- summarize.method[summarize == "sample"][1]
    if (is.na(method)) method <- "sumlog"
    shrink <- summarize.shrink[summarize == "sample"][1]
    if (is.na(shrink)) shrink <- Inf
    
    s.up <- sapply(1:n.samples, function(i) p_combine(up[,i], method, shrink)$p)
    s.down <- sapply(1:n.samples, function(i) p_combine(down[,i], method, shrink)$p)
    s.twoside <- sapply(1:n.samples, function(i) p_combine(twoside[,i], method, shrink)$p)
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
#' @param scale.degree Logical variable indicating whether the degree values are scaled according to the dropout rate.
#'
#' @return This function will return a list with the following components:
#'   \item{up}{A numeric vector containing the p-values of up-regulation test.}
#'   \item{down}{A numeric vector containing the p-values of down-regulation test.}
#'   \item{twoside}{A numeric vector containing the p-values of twoside test.}
#' 
#' @export
netDEG_pvalue <- function(ref.ratio.dist, expr.val, log.expr = FALSE, scale.degree = FALSE)
{
  if (!log.expr) expr.val <- log(expr.val)
  score <- get_diff_ratio_net(ref.ratio.dist, expr.val, log.expr = TRUE, scale.degree = scale.degree)$diff
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
#' @param scale.degree Logical variable indicating whether the degree values are scaled according to the dropout rate.
#' @param use.parallel Logical variable indicating to use the BiocParallel package to accelerate computation.
#' 
#' @return This function will return a list with the following components:
#'   \item{LB}{A numeric matrix with element [i,j] represents the lower quantile of expressioin ratios for gene pairs (i, j).}
#'   \item{NB}{A numeric vector with two elements: \code{size} and \code{mu}, 
#'   which are the estimated parameters of negative binomial distribution.}
#'   \item{p.edge}{The used input parameter \code{p.edge}.}
#'   
#' @export
get_ratio_distribution <- function(ref.expr.matrix, p.edge = 0.1, log.expr = FALSE, scale.degree = FALSE, use.parallel = FALSE)
{
  if (use.parallel && requireNamespace("BiocParallel")) {
    lapply <- BiocParallel::bplapply
  }

  if (!log.expr) ref.expr.matrix <- log(ref.expr.matrix)
  if (use.parallel) {
    n.genes <- dim(ref.expr.matrix)[1]
    dist <- lapply(1:ceiling((n.genes-1)/2), function (i) .Call(ND_RatioDistributionParI, ref.expr.matrix, p.edge, i))
    dist <- .Call(ND_RatioDistributionParM, dist, n.genes)
    dist$p.edge <- p.edge
  }
  else {
    dist <- .Call(ND_RatioDistribution, ref.expr.matrix, p.edge)
  }
  diff <- unlist(lapply(1:dim(ref.expr.matrix)[2], function(i) get_diff_ratio_net(dist, ref.expr.matrix[,i], log.expr = TRUE, scale.degree = scale.degree)$diff))
  diff <- diff[!is.na(diff)]
  if (requireNamespace("fitdistrplus")) {
    dist$NB <- fitdistrplus::fitdist(abs(diff), "nbinom")$estimate
  }
  else {
    dist$NB <- MASS::fitdistr(abs(diff), "negative binomial")$estimate
  }
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
#' @param scale.degree Logical variable indicating whether the degree values are scaled according to the dropout rate.
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
get_ratio_distribution2 <- function(ref.expr.matrix, p.edge = 0.1, p.trim = 0.3, log.expr = FALSE, scale.degree = FALSE, use.parallel = FALSE)
{
  if (use.parallel && requireNamespace("BiocParallel")) {
    lapply <- BiocParallel::bplapply
  }
  
  if (!log.expr) ref.expr.matrix <- log(ref.expr.matrix)
  dist <- .Call(ND_RatioDistribution2, ref.expr.matrix, p.edge, p.trim)
  diff <- unlist(lapply(1:dim(ref.expr.matrix)[2], function(i) get_diff_ratio_net(dist, ref.expr.matrix[,i], log.expr = TRUE, scale.degree = scale.degree)$diff))
  diff <- diff[!is.na(diff)]
  if (requireNamespace("fitdistrplus")) {
    dist$NB <- fitdistrplus::fitdist(abs(diff), "nbinom")$estimate
  }
  else {
    dist$NB <- MASS::fitdistr(abs(diff), "negative binomial")$estimate
  }
  dist
}

#' Construct differential expression ratio network
#' 
#' Construct the differential expression ratio network for a single sample.
#' 
#' @param ref.ratio.dist The expression ratio distribution profile returned by \code{get_ratio_distribution} or \code{get_ratio_distribution2}.
#' @param expr.val Numeric vector of gene expression values in the sample.
#' @param log.expr Logical variable indicating whether the input expression vector is in logarithmic scale.
#' @param scale.degree Logical variable indicating whether the degree values are scaled according to the dropout rate.
#' 
#' @return This function will return a list with the following components:
#'   \item{net}{The binary adjacent matrix of differential expression ratio network.}
#'   \item{diff}{A numeric vector containing the adjusted degree differences of all genes.}
#'   \item{degree}{A list containing the raw degree differences and sums of all genes.}
#' 
#' @export
get_diff_ratio_net <- function(ref.ratio.dist, expr.val, log.expr = FALSE, scale.degree = FALSE)
{
  if (!log.expr) expr.val <- log(expr.val)
  edges <- .Call(ND_DiffRatioNet, ref.ratio.dist$LB, expr.val)
  n <- length(expr.val)
  net <- sparseMatrix(i = edges$i, j = edges$j, dims = c(n, n))
  d <- get_adjusted_deg_diff(net, expr.val, scale.degree)
  list(net = net, diff = d$diff, degree = d$degree)
}


#' Calculate adjusted degree differences for given network
#' 
#' Calculate the adjusted degree differences for all genes in the given network.
#' 
#' @param net The binary adjacent matrix of differential expression ratio network.
#' @param log.expr.val Numeric vector containing the logarithmic scale gene expression values.
#' @param scale.degree Logical variable indicating whether the degree values are scaled according to the dropout rate.
#' @param p The parameter for calculating the adjusted degree differences.
#' 
#' @return This function will return a list with the following components:
#'   \item{diff}{A numeric vector containing the adjusted degree differences of all genes.}
#'   \item{degree}{A list containing the raw degree differences and sums of all genes.}
#' 
#' @export
get_adjusted_deg_diff <- function(net, log.expr.val, scale.degree = FALSE, p = 0.5)
{
  g.NA <- !is.finite(log.expr.val)
  d.out <- rowSums(net)
  d.in <- colSums(net)
  if (scale.degree) {
    f <- length(g.NA) / sum(!g.NA)
    d.out <- round(d.out * f)
    d.in <- round(d.in * f)
  }
  d.out[g.NA] <- NA
  d.in[g.NA] <- NA
  d.sum <- d.out + d.in
  d.diff <- d.out - d.in
  adj.diff <- d.diff - ceiling(stats::median(d.diff[d.sum <= stats::quantile(d.sum, p, na.rm = TRUE)], na.rm = TRUE))
  list(diff = adj.diff, degree = list(diff = d.diff, sum = d.sum))
}


#' Calculate combined p-value
#' 
#' Combine the statistical significance results from several independent tests by using one of several methods.
#' 
#' @param p the numeric vector containing the p-values need to combine.
#' @param method the method use to combine the p-values, can be "sumlog" (Fisher's method), "sumz" (Stouffer’s method).
#' @param shrink the number of p-values used in calculation, which are uniform selected from original p-value vector.
#' 
#' @return This function will return a list with the following components:
#'   \item{p}{The combined p-value.}
#'   \item{v}{The value of statistic.}
#'   \item{}{Use "sumlog" method:}
#'   \item{chisq}{The value of chi-squared statistic.}
#'   \item{df}{The degrees of freedom of chi-squared distribution.}
#'   \item{}{Use "sumz" method:}
#'   \item{z}{The value of sum z statistic.}
#' 
#' @export
p_combine <- function(p, method = "sumlog", shrink = Inf)
{
  p <- p[!is.na(p)]
  p[p > 1] <- 1
  p[p < .Machine$double.xmin] <- .Machine$double.xmin
  if (shrink > 0 && length(p) > shrink) {
    p <- sort(p)[ceiling((2*(1:shrink)-1) / (2*shrink) * length(p))]
  }
  if (length(p) == 0) p <- 0.5
  if (method == "sumlog") {
    chisq <- (-2) * sum(log(p))
    df <- 2 * length(p)
    p <- stats::pchisq(chisq, df, lower.tail = FALSE)
    if (p == 0) p <- .Machine$double.xmin / chisq
    list(chisq = chisq, df = df, v = chisq, p = p)
  }
  else if (method == "sumz") {
    z <- sum(stats::qnorm(p, lower.tail = FALSE)) / sqrt(length(p))
    p <- stats::pnorm(z, lower.tail = FALSE)
    if (p == 0) p <- .Machine$double.xmin / z
    list(z = z, v = z, p = p)
  }
  else {
    stop("Unknown method in function p_combine: ", method)
  }
}
