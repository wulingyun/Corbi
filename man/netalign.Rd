\name{net.align}
\alias{net.align}

\title{Network alignment method based on conditional random fields}

\description{
Find the maximal matching subnetworks from a target network for a query network
based on the conditional random fields (CRF) model.}

\usage{
net.align(query.net, target.net, node.sim, query.type=4,
	delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1,
	output="result.txt")
}

\arguments{
\item{query.net}{The input file name of the query network.}

\item{target.net}{The input file name of the target network.}

\item{node.sim}{The input file name of the node similarity scores between the query network and the target network.}

\item{query.type}{The querying network type: 1 - general, 2 - chain, 3 - tree, 4 - heuristic.}

\item{delta.d}{The parameter delta.d is a parameter for deletions.}

\item{delta.c}{The parameter delta.c is a parameter for consecutive deletions.}

\item{delta.e}{The parameter delta.e is a parameter for single deletion.}

\item{delta.s}{The parameter delta.s is a parameter for insertions.}

\item{output}{The suffix of output file name. The output contains two files in the working directory. One is the matching nodes and edges between query network and target network, the other is the unique matching node pairs.}
}

\details{
This is an approach for network alignment problem based on conditional random field (CRF)
model which uses the node similarity and structure information equally. This method is based on our network querying method \code{\link{net.query}}. This method uses an iterative  
strategy to get the one-to-one map between the query network and target netowrk.

More details can be seen in \code{\link{net.query}}.
}

\references{Qiang Huang, Ling-Yun Wu, and Xiang-Sun Zhang. 
CNetA: A new method for network alignment
combined biological and structural features. ISB2012, submitted.}


\examples{
\dontrun{
library(Corbi)

## An example: "querynet.txt", "targetnet.txt", "nodesim.txt" are 
## three input files in the working directory
net.align("querynet.txt", "targetnet.txt", "nodesim.txt")
}}
