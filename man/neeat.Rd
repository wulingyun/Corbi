\name{neeat}
\alias{neeat}
\title{Network enhanced enrichment analysis tool}
\usage{
neeat(core.sets, gene.set = NULL, net = NULL, subnet = NULL,
  depths = NULL, method = "gene", rho = 0.5, max.depth = 10,
  n.perm = 10000, use.multinom = FALSE, z.threshold = 2,
  adjust.p = "BH", n.cpu = 1, batch.size = 5000)
}
\arguments{
  \item{core.sets}{Logical matrix indicated the core genes
  associated with specific functions or pathways.}

  \item{gene.set}{Logical vector indicated the gene set for
  evaluating.}

  \item{net}{The adjacent matrix of network.}

  \item{subnet}{The adjacent matrix of sub-network for
  evaluating in the NEEAT-subnet model.}

  \item{method}{A string indicated the NEEAT model,
  including "gene", "net" and "subnet". Use "hyper" for
  traditional hypergeometric test, in which the network
  information is ignored.}

  \item{rho}{The weight parameter for depths.}

  \item{max.depth}{Integer for the maximum depth considered
  in the NEEAT models.}

  \item{n.perm}{The number of permutations for calculating
  p-values.}

  \item{use.multinom}{Logical variable indicated whether
  use \code{\link{rmultinom}} to approximate
  \code{\link{rmultihyper}}.}

  \item{z.threshold}{The threshold for filtering out small
  Z-scores. The p-values are calculated only for the
  Z-scores \code{>= z.threshold}, otherwise p-values are
  set as 1.}

  \item{adjust.p}{The method for adjusting p-values for
  multiple testing. Use "none" for bypassing. See
  \code{\link{p.adjust}} for available methods.}

  \item{n.cpu}{The number of CPUs/cores used in the
  parallel computation.}

  \item{batch.size}{The desired size of batches in the
  parallel computation.}
}
\value{
This function will return a matrix of same columns as
\code{core.sets}, and each column containing the following
components for the correponding core gene set
\code{core.sets[,i]}: \item{\code{z.score}}{The Z-score for
\code{gene.set}} \item{\code{p.value}}{The statistic
significance for \code{gene.set} under specified NEEAT
model} \item{\code{raw.score}}{The raw score for
\code{gene.set} under specified NEEAT model}
\item{\code{avg.score}}{The average score for random
permutations of \code{gene.set}}
\item{\code{var.score}}{The variance of scores for random
permutations of \code{gene.set}}
}
\description{
Functional enrichment anslysis of gene sets by integrating
network information
}
\details{
This is a novel functional enrichment analysis tool based
on the gene set and network information.
}
\seealso{
\code{\link{get_core_sets}}
}

