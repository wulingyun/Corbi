\name{rmultihyper}
\alias{rmultihyper}
\title{The Multivariate Hypergeometric Distribution}
\usage{
rmultihyper(n, m, k)
}
\arguments{
  \item{n}{The number of observations}

  \item{m}{The vector containing the number of balls in
  each urn}

  \item{k}{The number of balls drawn from the urns}
}
\value{
This function will return a vector of the same length as
\code{m}, containing the number of balls drawn from each
urn
}
\description{
Generate random variables for the multivariate
hypergeometric distribution
}
\details{
This function generates random variables for the
multivariate hypergeometric distribution by iteratively
calling hypergeometric random variable generator
\code{\link{rhyper}}.
}
\seealso{
\code{\link{rhyper}}
}

