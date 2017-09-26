
#' @export
make_DEG_pattern <- function(n.genes, n.samples, fold.change = 2, gene.rate = 0.3, sample.rate = 1.0, active.rate = 1.0, up.rate = 0.5)
{
  active.gene <- runif(n.genes) <= gene.rate
  active.sample <- runif(n.samples) <= sample.rate
  if (sum(active.gene) < 1) active.gene[sample.int(n.genes, 1)] <- TRUE
  if (sum(active.sample) < 1) active.sample[sample.int(n.samples, 1)] <- TRUE
  up.gene <- logical(n.genes)
  up.gene[active.gene] <- runif(sum(active.gene)) <= up.rate
  down.gene <- active.gene & !up.gene
  fc.sub <- matrix(fold.change, sum(active.gene), sum(active.sample))
  fc.sub[down.gene[active.gene], ] <- 1/fold.change
  fc.sub[runif(length(fc.sub)) > active.rate] <- 1
  fc <- matrix(1, n.genes, n.samples)
  fc[active.gene, active.sample] <- fc.sub
  return(list(FC=fc, gene=up.gene-down.gene, sample=as.integer(active.sample)))
}


#' @export
make_DEG_data <- function(n.genes, n.samples.A, n.samples.B, exp.mean = 8, exp.sd = 2, alpha = 0.2, size.factor.sd = 0.1, ...)
{
  # simulate DEG and heterogeneity
  deg <- make_DEG_pattern(n.genes, n.samples.B, ...)
  # simulate expression mean and dispersion
  mu0 <- 2^rnorm(n.genes, exp.mean, exp.sd)
  sd0 <- alpha * mu0
  # simulate group A
  sfA <- 2^rnorm(n.samples.A, sd = size.factor.sd)
  countsA <- matrix(rnorm(n.genes * n.samples.A, mean = mu0 %*% t(sfA), sd = sd0), nrow = n.genes)
  # simulate group B
  sfB <- 2^rnorm(n.samples.B, sd = size.factor.sd)
  countsB <- matrix(rnorm(n.genes * n.samples.B, mean = mu0 %*% t(sfB) * deg$FC, sd = sd0 * deg$FC), nrow = n.genes)
  countsA[countsA < 0] <- 0
  countsB[countsB < 0] <- 0
  return(list(DEG=deg, countsA=countsA, countsB=countsB))
}


#' @export
make_DEG_data2 <- function(n.genes, n.samples.A, n.samples.B, exp.mean = 8, exp.sd = 2, dispersion = 0.2, size.factor.sd = 0.1, ...)
{
# simulate DEG and heterogeneity
  deg <- make_DEG_pattern(n.genes, n.samples.B, ...)
# simulate expression mean
  mu0 <- 2^rnorm(n.genes, exp.mean, exp.sd)
  if (is.null(dispersion)) dispersion <- 4/mu0 + 0.1
# simulate group A
  sfA <- 2^rnorm(n.samples.A, sd = size.factor.sd)
  countsA <- matrix(rnbinom(n.genes * n.samples.A, mu = mu0 %*% t(sfA), size = 1/dispersion), nrow = n.genes)
# simulate group B
  sfB <- 2^rnorm(n.samples.B, sd = size.factor.sd)
  countsB <- matrix(rnbinom(n.genes * n.samples.B, mu = mu0 %*% t(sfB) * deg$FC, size = 1/dispersion), nrow = n.genes)
  return(list(DEG=deg, countsA=countsA, countsB=countsB))
}
