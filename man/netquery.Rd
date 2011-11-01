\name{net.query}
\alias{net.query}
\alias{net.query.batch}

\title{Network querying method based on conditional random fields}

\description{
Find the best matching subnetworks from a large target network for small query networks
based on the conditional random fields (CRF) model.}

\usage{
net.query(query.net, target.net, node.sim, query.type=1,
	delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1,
	output="result.txt")
net.query.batch(query.nets, target.net, node.sim, query.type=1,
	delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1,
	output="result.txt")
}

\arguments{
\item{query.net}{The input file name of the query network.}

\item{query.nets}{The vector of input file names of the query networks.}

\item{target.net}{The input file name of the target network.}

\item{node.sim}{The input file name of the node similarity scores between the query network and the target network.}

\item{query.type}{The querying network type: 1 - general, 2 - chain, 3 - tree.}

\item{delta.d}{The parameter delta.d is a parameter for deletions.}

\item{delta.c}{The parameter delta.c is a parameter for consecutive deletions.}

\item{delta.e}{The parameter delta.e is a parameter for single deletion.}

\item{delta.s}{The parameter delta.s is a parameter for insertions.}

\item{output}{The suffix of output file name.}
}

\details{
This is an approach for network querying problem based on conditional random field (CRF)
model which can handle both undirected and directed networks, acyclic and cyclic networks,
and any number of insertions/deletions.

When querying several networks in the same target network, \code{\link{net.query.batch}} will save much time.

\itemize{
\item query.net: The query network file is written as follows:\cr
v1 v2 v3 v4 v5\cr
v3 v4 \cr
...   \cr
where v1, v2, v3, v4, v5 ... are the nodes' names and each line indicates there are edges between
the first node and other nodes in the line. For example, the first line denotes 4 edges:
(v1, v2), (v1, v3), (v1, v4), and (v1, v5).

\item target.net: The format of this file is the same as the query network file.

\item node.sim: This similarity file's format is as follows:\cr
v1 V1 s1 \cr
v1 V2 s2 \cr
...      \cr
v1 is the node from the query network, V1 is the node from the target network,
s1 is the similarity score between the node v1 and V1, and so on.

\item query.type: If the query is a general graph with loops, we set query.type = 1.
If the query is a chain, we set query.type = 2. If the query is a tree,
we set query.type = 3. The default value is 1.

\item delta.d:  The smaller delta.d is, the heavier penalty for deletions. 

\item delta.c: The smaller delta.c is, the heavier penalty for consecutive deletions.

\item delta.e: The smaller delta.e is, the heavier penalty for single deletion.

\item delta.s: The larger delta.s indicates heavier penalty for insertions.

}
}

\references{Qiang Huang, Ling-Yun Wu, and Xiang-Sun Zhang. 
An Efficient Network Querying Method Based on Conditional Random Fields. 
Bioinformatics, 2011, doi:10.1093/bioinformatics/BTR524.}


\examples{
\dontrun{
library(Corbi)

## An example: "querynet.txt", "targetnet.txt", "nodesim.txt" are 
## three input files in the working directory
net.query("querynet.txt", "targetnet.txt", "nodesim.txt", query.type=3)

## Batch example
net.query.batch(c("querynet.txt", "querynet2.txt"),
	"targetnet.txt", "nodesim.txt", query.type=3)
}}
