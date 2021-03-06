% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subnets.R
\name{extend_subnets}
\alias{extend_subnets}
\title{Extend subnetworks from smaller subnetworks}
\usage{
extend_subnets(subnet1, subnet2, size = 0)
}
\arguments{
\item{subnet1}{The matrix representing the first set of subnetworks}

\item{subnet2}{The matrix representing the second set of subnetworks}

\item{size}{The desired size of extended subnetworks}
}
\value{
A matrix represents the extended subnetworks, in which each row represents a subnetwork.
}
\description{
Extend subnetworks by pairwise overlapping two sets of smaller subnetworks.
}
\details{
Enumerate all possible subnetworks of desired size by pairwise overlapping two sets of 
subnetworks of size \code{s1} and \code{s2}. The desired size should be between 
\code{max(s1,s2)+1} and \code{s1+s2-1}. Invalid desired size will be replaced by the 
minimum allowed value \code{max(s1,s2)+1}.
}
\examples{

library(Corbi)
net <- matrix(FALSE, nrow=10, ncol=10)
net[sample.int(100, 20)] <- TRUE
net <- net | t(net)
subnets <- get_subnets(net, 3)
subnets[[4]] <- extend_subnets(subnets[[3]], subnets[[2]], 4)

}
