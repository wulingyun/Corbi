makeDiffExpSeqCount <- function(n.genes, n.samples.A, n.samples.B, fold.change = 2, deg.rate = 0.3, homo.rate = 1, dispersion = 0.2, size.factor.sd = 0.1)
{
# simulate DEG and heterogeneity
  mu0 <- rexp(n.genes, rate = 1/250)
  deg <- runif(n.genes) <= deg.rate
  fc <- matrix(1, n.genes, n.samples.B)
  fc[deg, ] <- fold.change
  fc[runif(length(fc)) > homo.rate] <- 1
# simulate group A
  sfA <- exp(rnorm(n.samples.A, sd = size.factor.sd))
  countsA <- t(sapply(seq_len(n.genes), 
                      function(i) sapply(seq_len(n.samples.A), 
                                         function(j) rnbinom(1, mu = mu0[i] * sfA[j], size = 1/dispersion))))
# simulate group B
  sfB <- exp(rnorm(n.samples.B, sd = size.factor.sd))
  countsB <- t(sapply(seq_len(n.genes),
                      function(i) sapply(seq_len(n.samples.B),
                                         function(j) rnbinom(1, mu = mu0[i] * sfB[j] * fc[i, j], size = 1/dispersion))))
  return(list(DEG = deg, FC = fc, countsA = countsA, countsB = countsB))
}
