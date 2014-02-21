\name{rmultihyper}
\alias{rmultihyper}
\title{The Multivariate Hypergeometric Distribution}
\usage{
rmultihyper(n, k, m)
}
\arguments{
  \item{n}{The number of observations}

  \item{k}{The number of balls drawn from the urns}

  \item{m}{The integer vector containing the number of
  balls in each urn}
}
\value{
This function will return a matrix of \code{length(m)} rows
and \code{n} columns, containing the number of balls drawn
from each urn
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

