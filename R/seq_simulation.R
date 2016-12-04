
makeDiffExpSeqCount <- function(n.genes, n.samples.A, n.samples.B, fold.change = 2, deg.rate = 0.3, deg.homo = 1, dispersion = 0.2, size.factor.sd = 0.1)
{
# simulate DEG and heterogeneity
  mu0 <- rexp(n.genes, rate = 1/250)
  deg <- runif(n.genes) <= deg.rate
  fc <- matrix(1, n.genes, n.samples.B)
  fc[deg, ] <- fold.change
  fc[runif(length(fc)) > deg.homo] <- 1
# simulate group A
  sfA <- exp(rnorm(n.samples.A, sd = size.factor.sd))
  countsA <- matrix(rnbinom(n.genes * n.samples.A, mu = mu0 %*% t(sfA), size = 1/dispersion), nrow = n.genes)
# simulate group B
  sfB <- exp(rnorm(n.samples.B, sd = size.factor.sd))
  countsB <- matrix(rnbinom(n.genes * n.samples.B, mu = mu0 %*% t(sfB) * fc, size = 1/dispersion), nrow = n.genes)
  return(list(DEG = deg, FC = fc, countsA = countsA, countsB = countsB))
}


makeDiffExpSeqCount2 <- function(n.genes, n.samples.A, n.samples.B, fold.change = 2, deg.rate = 0.3, deg.homo = 1, exp.mean = 8, exp.sd = 2, size.factor.sd = 0.1)
{
# simulate expression mean and dispersion
  mu0 <- 2^rnorm(n.genes, exp.mean, exp.sd)
  dispersion <- 4/mu0 + 0.1
# simulate DEG and heterogeneity
  deg <- runif(n.genes) <= deg.rate
  fc <- matrix(1, n.genes, n.samples.B)
  fc[deg, ] <- fold.change
  fc[runif(length(fc)) > deg.homo] <- 1
# simulate group A
  sfA <- 2^rnorm(n.samples.A, sd = size.factor.sd)
  countsA <- matrix(rnbinom(n.genes * n.samples.A, mu = mu0 %*% t(sfA), size = 1/dispersion), nrow = n.genes)
# simulate group B
  sfB <- 2^rnorm(n.samples.B, sd = size.factor.sd)
  countsB <- matrix(rnbinom(n.genes * n.samples.B, mu = mu0 %*% t(sfB) * fc, size = 1/dispersion), nrow = n.genes)
  return(list(DEG = deg, FC = fc, countsA = countsA, countsB = countsB))
}
