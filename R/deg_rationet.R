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
#' @param adjust.p Logical variable indicating whether the individual-level p-values are adjusted by "BH" method.
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
netDEG <- function(ref.expr.matrix, expr.matrix, p.edge = 0.1, n.ref.genes = 10000,
                   summarize = c("gene", "sample"), summarize.method = c("sumlog", "sumlog"), summarize.shrink = c(Inf, Inf),
                   log.expr = FALSE, zero.as.dropout = TRUE, scale.degree = TRUE, adjust.p = TRUE, use.parallel = FALSE)
{
  if (use.parallel && requireNamespace("BiocParallel")) {
    lapply <- BiocParallel::bplapply
  }

  if (!log.expr) {
    ref.expr.matrix <- log(ref.expr.matrix)
    expr.matrix <- log(expr.matrix)
  }
  if (!zero.as.dropout) {
    zero.expr = min(ref.expr.matrix[is.finite(ref.expr.matrix)], 0, na.rm = TRUE) - log(2)
    ref.expr.matrix[!is.finite(ref.expr.matrix)] <- zero.expr
    expr.matrix[!is.finite(expr.matrix)] <- zero.expr
  }
  n.samples <- dim(expr.matrix)[2]
  
  message("Estimating the ratio distribution from reference samples")
  dist <- get_ratio_distribution(ref.expr.matrix, p.edge, n.ref.genes = n.ref.genes, log.expr = TRUE, scale.degree = scale.degree, use.parallel = use.parallel)

  message("Calculating the sample-specific p-values for test samples")
  p <- lapply(1:n.samples, function(i) netDEG_pvalue(dist, expr.matrix[,i], log.expr = TRUE, scale.degree = scale.degree))
  rm(dist)
  
  up <- sapply(p, function(i) i$up)
  down <- sapply(p, function(i) i$down)
  twoside <- sapply(p, function(i) i$twoside)
  dimnames(up) <- dimnames(down) <- dimnames(twoside) <- dimnames(expr.matrix)
  rm(p)

  if (adjust.p) {
    up <- t(apply(up, 1, function(x) p.adjust(x, "BH")))
    down <- t(apply(down, 1, function(x) p.adjust(x, "BH")))
    twoside <- t(apply(twoside, 1, function(x) p.adjust(x, "BH")))
  }

  results <- list(up = up, down = down, twoside = twoside)
  
  if ("gene" %in% summarize) {
    n.genes <- dim(expr.matrix)[1]
    n.refs <- dim(ref.expr.matrix)[2]

    message("Estimating the ratio distribution from test samples")
    dist <- get_ratio_distribution(expr.matrix, p.edge, n.ref.genes = n.ref.genes, log.expr = TRUE, scale.degree = scale.degree, use.parallel = use.parallel)

    message("Calculating the sample-specific p-values for reference samples")
    p <- lapply(1:n.refs, function(i) netDEG_pvalue(dist, ref.expr.matrix[,i], log.expr = TRUE, scale.degree = scale.degree))
    rm(dist)

    rev.up <- sapply(p, function(i) i$up)
    rev.down <- sapply(p, function(i) i$down)
    rev.twoside <- sapply(p, function(i) i$twoside)
    dimnames(rev.up) <- dimnames(rev.down) <- dimnames(rev.twoside) <- dimnames(ref.expr.matrix)
    rm(p)

    if (adjust.p) {
      rev.up <- t(apply(rev.up, 1, function(x) p.adjust(x, "BH")))
      rev.down <- t(apply(rev.down, 1, function(x) p.adjust(x, "BH")))
      rev.twoside <- t(apply(rev.twoside, 1, function(x) p.adjust(x, "BH")))
    }

    results$rev <- list(up = rev.up, down = rev.down, twoside = rev.twoside)

    message("Gene-wise summarization")
    
    method <- summarize.method[summarize == "gene"][1]
    if (is.na(method)) method <- "sumlog"
    shrink <- summarize.shrink[summarize == "gene"][1]
    if (is.na(shrink)) shrink <- Inf
    
    g.up <- cbind(up, rev.down)
    g.down <- cbind(down, rev.up)
    g.twoside <- cbind(twoside, rev.twoside)
    g.up <- sapply(1:n.genes, function(i) p_combine(g.up[i,], method, shrink)$p)
    g.down <- sapply(1:n.genes, function(i) p_combine(g.down[i,], method, shrink)$p)
    g.twoside <- sapply(1:n.genes, function(i) p_combine(g.twoside[i,], method, shrink)$p)
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
get_ratio_distribution <- function(ref.expr.matrix, p.edge = 0.1, n.ref.genes = 10000, log.expr = FALSE, scale.degree = FALSE, use.parallel = FALSE)
{
  if (use.parallel && requireNamespace("BiocParallel")) {
    lapply <- BiocParallel::bplapply
  }

  if (!log.expr) ref.expr.matrix <- log(ref.expr.matrix)
  n.non0 <- rowSums(is.finite(ref.expr.matrix))
  n.genes <- dim(ref.expr.matrix)[1]
  
  if (n.ref.genes <= 0) n.ref.genes <- n.genes
  n.ref.genes <- min(n.ref.genes, n.genes, sum(n.non0 > 0))
  ref.genes <- (1:n.genes) %in% sample.int(n.genes, n.ref.genes, replace = F, prob = n.non0)
  map.genes <- integer(n.genes)
  map.genes[ref.genes] <- cumsum(ref.genes)[ref.genes]
  map.genes[!ref.genes] <- cumsum(!ref.genes)[!ref.genes]
  dist <- list(n.ref.genes = n.ref.genes, ref.genes = ref.genes, map.genes = map.genes)

  ref.expr.matrix.A <- ref.expr.matrix[ref.genes, , drop = FALSE]
  ref.expr.matrix.B <- ref.expr.matrix[!ref.genes, , drop = FALSE]
  if (use.parallel) {
    dist$LB0 <- lapply(1:ceiling((n.ref.genes-1)/2), function (i) .Call(ND_RatioDistributionParI, ref.expr.matrix.A, p.edge, i))
    dist$LB0 <- .Call(ND_ParMerge, dist$LB0, n.ref.genes, -Inf, FALSE)
  }
  else {
    dist$LB0 <- .Call(ND_RatioDistribution, ref.expr.matrix.A, p.edge)
  }
  if (n.ref.genes < n.genes) {
    if (use.parallel) {
      n <- n.genes - n.ref.genes
      temp <- sapply(1:n.ref.genes, function(i) .Call(ND_RatioDistributionParAiB, ref.expr.matrix.A[i, ], ref.expr.matrix.B, p.edge))
      dist$LB1 <- t(temp[1:n, , drop = F])
      dist$LB2 <- temp[(n+1):(2*n), , drop = F]
    }
    else {
      temp <- .Call(ND_RatioDistributionAB, ref.expr.matrix.A, ref.expr.matrix.B, p.edge)
      dist$LB1 <- temp$LB1
      dist$LB2 <- temp$LB2
    }
  }
  else {
    dist$LB1 <- NULL
    dist$LB2 <- NULL
  }
  
  diff <- unlist(lapply(1:dim(ref.expr.matrix)[2], function(i) get_diff_ratio_net(dist, ref.expr.matrix[,i], log.expr = TRUE, scale.degree = scale.degree)$diff))
  diff <- diff[!is.na(diff)]
  if (requireNamespace("fitdistrplus")) {
    dist$NB <- fitdistrplus::fitdist(abs(diff), "nbinom")$estimate
  }
  else {
    dist$NB <- MASS::fitdistr(abs(diff), "negative binomial")$estimate
  }
  
  dist$p.edge <- p.edge
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
  edges <- .Call(ND_DiffRatioNet, ref.ratio.dist, expr.val)
  n <- length(expr.val)
  net <- sparseMatrix(i = edges$i, j = edges$j, dims = c(n, n))
  d <- get_adjusted_deg_diff(net, expr.val, ref.ratio.dist, scale.degree)
  list(net = net, diff = d$diff, degree = d$degree)
}


#' Calculate adjusted degree differences for given network
#' 
#' Calculate the adjusted degree differences for all genes in the given network.
#' 
#' @param net The binary adjacent matrix of differential expression ratio network.
#' @param log.expr.val Numeric vector containing the logarithmic scale gene expression values.
#' @param ref.ratio.dist The expression ratio distribution profile returned by \code{get_ratio_distribution} or \code{get_ratio_distribution2}.
#' @param scale.degree Logical variable indicating whether the degree values are scaled according to the dropout rate.
#' @param p The parameter for calculating the adjusted degree differences.
#' 
#' @return This function will return a list with the following components:
#'   \item{diff}{A numeric vector containing the adjusted degree differences of all genes.}
#'   \item{degree}{A list containing the raw degree differences and sums of all genes.}
#' 
#' @export
get_adjusted_deg_diff <- function(net, log.expr.val, ref.ratio.dist, scale.degree = FALSE, p = 0.5)
{
  g.NA <- !is.finite(log.expr.val)
  d.out <- rowSums(net[, ref.ratio.dist$ref.genes])
  d.in <- colSums(net[ref.ratio.dist$ref.genes, ])
  if (scale.degree) {
    f <- length(g.NA) / sum(!g.NA & ref.ratio.dist$ref.genes)
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
