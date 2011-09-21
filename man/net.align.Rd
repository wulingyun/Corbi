\name{net.align}
\alias{net.align}

\title{network querying method based on conditional random fields}

\description{
Find the best matching subnetworks from a large target network for small querying networks based on the conditional random fields (CRFs) model. }

\usage{
net.align(quenet,net,nodefea,qtype=2, ...)
net.align(quenet,net,nodefea,qtype=2,delta.d=1e-10,delta.c=0.5,delta.e=1,delta.s=1,isnormalize=0,multitest=1)
}

\arguments{
\item{quenet}{The file name of the querying network.}

\item{net}{The file name of the target network. The format of this file is the same as the querying network file.}

\item{nodefea}{The file name of the node similarity scores between the querying network and the target network.}

\item{qtype}{The querying network type.}

\item{delta.d}{The parameter delta.d is a parameter for deletions.}

\item{delta.c}{The parameter delta.c is  a parameter for consecutive deletions.}

\item{delta.e}{The parameter delta.e is a parameter for single deletion.}

\item{delta.s}{The parameter delta.s is a parameter for insertions.}

\item{isnormalize}{Whether carry out the normalized operation to the node similarity scores or not? default not.}

\item{multitest}{The result files' number;}
}

\details{
This is an approach for network querying problem based on conditional random field (CRF) model which can handle both undirected and directed networks, acyclic and cyclic networks, and any number of insertions/deletions.

\itemize{
\item quenet: The querying file is written as follows:\cr
v1 v2 \cr
v3 v4 \cr
...   \cr
where v1,v2,v3,v4, ... are the nodes' names and each line contains an edge in the network.

\item net: The format of this file is the same as the querying network file.

\item nodefea: This files' format is as follows:\cr
v1 V1 s1 \cr
v1 V2 s2 \cr
...      \cr
v1 is the node from the querying network, 
V1 is the node from the target network, 
s1 is the similarity score between the node v1 and V1, and so on.

\item qtype: If the querying is a chain, we set qtype=0; If the querying is a tree, we set qtype=1; If the querying is a graph with loops, we set qtype=2; with the default is 2.

\item delta.d:  The smaller delta.d is, the heavier penalty for deletions. 

\item delta.c: The smaller delta.c is, the heavier penalty for consecutive deletions.

\item delta.e: The smaller delta.e is, the heavier penalty for single deletion.

\item delta.s: The larger delta.s indicates heavier penalty for insertions.

\item isnormalize: Whether carry out the normalized operation to the node similarity scores or not? Default not.

\item multitest: For multiple test examples, we need a number to distinguish the different result files.
}
}

\references{Qiang Huang, Ling-Yun Wu, and Xiang-Sun Zhang. 
An Efficient Network Querying Method Based on Conditional Random Fields. 
Bioinformatics, 2011, doi:10.1093/bioinformatics/BTR524.}


\examples{
library(CRF)
## Schematic of an example: "quenet.txt", "ynet.txt", "queblast.txt" are 
## the file names in the work directory
net.align("quenet.txt","ynet.txt","queblast.txt",qtype=1) 

}