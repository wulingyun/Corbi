#' Read network information from text file
#' 
#' Read the network information from a text file with specific format.
#' 
#' This function reads the network information from a text file with specific format:
#' each line contains two strings separated by spaces, which correspond to the 
#' names of two end points of one edge in the network.
#' 
#' @param file The name of text file
#' @return A list with the following components:
#'   \item{size}{The number of network nodes}
#'   \item{node}{The vector of network node names}
#'   \item{matrix}{The logical adjacency matrix}
#' 
#' @seealso \code{\link{write_net}}
#' 
#' @import Matrix
#' 
#' @export
read_net <- function(file)
{
  net.text <- as.matrix(utils::read.table(file, fill=TRUE, as.is=TRUE, col.names=1:max(utils::count.fields(file))))
  net.node <- unique(as.character(net.text))
  net.node <- net.node[net.node != ""]
  net.edge <- cbind(as.character(net.text[,1]), as.character(net.text[,-1]))
  net.edge <- net.edge[net.edge[,2] != "", ]
  net.size <- length(net.node)
  node.id <- seq_along(net.node)
  names(node.id) <- net.node
  net.matrix <- sparseMatrix(node.id[net.edge[,1]], node.id[net.edge[,2]], x=TRUE, dims=c(net.size, net.size), dimnames=list(net.node, net.node))
  list(size=net.size, node=net.node, matrix=net.matrix)
}

#' Write network information to text file
#' 
#' Write the network information to a text file with specific format.
#' 
#' This function writes the network information to a text file with specific format:
#' each line contains two strings separated by spaces, which correspond to the
#' names of two end points of one edge in the network.
#' 
#' @param net A list as returned by \code{\link{read_net}}
#' @param file The name of text file
#' 
#' @seealso \code{\link{read_net}}
#' 
#' @import Matrix
#' 
#' @export
write_net <- function(net, file)
{
  net.edge <- which(net$matrix != 0, arr.ind=1)
  net.edge <- matrix(net$node[net.edge], ncol=2)
  utils::write.table(net.edge, file, quote=FALSE, row.names=FALSE, col.names=FALSE)
}


#' Calculate shortest distances of unweighted network
#' 
#' Calculate all pairs of shortest distances of unweighted network
#' 
#' This function calculates all pairs of shortest distances of unweighted network
#' by using breadth-first-search (BFS) algorithm.
#' 
#' @param net.matrix Logical adjacency matrix of given unweighted network
#' @param source.nodes Logical vector to indicate the source nodes that 
#' need to calculate the shortest distances
#' @return This function will return the shortest distance matrix, where the element
#' \code{[i, j]} is the shortest distance between node i and j. Value -1 means unreachable.
#' If \code{source.nodes[i]} equals FALSE, the shortest distance from i to other nodes
#' will not be calculated and the row i will be all -1.
#' 
#' @import Matrix
#' 
#' @export
get_shortest_distances <- function(net.matrix, source.nodes = rep_len(TRUE, dim(net.matrix)[1]))
{
  edges <- which(net.matrix != 0, arr.ind = TRUE)
  edges <- edges[order(edges[,1]), ]
  index <- findInterval(0:dim(net.matrix)[1], edges[,1])
  .Call(NQ_ShortestDistances, edges, index, source.nodes)
}


read_sim <- function(file)
{
  sim.text <- utils::read.table(file, as.is=TRUE)
  sim.node1 <- unique(as.character(sim.text[,1]))
  sim.node2 <- unique(as.character(sim.text[,2]))
  sim.size1 <- length(sim.node1)
  sim.size2 <- length(sim.node2)
  sim.matrix <- matrix(0, sim.size1, sim.size2, dimnames=list(sim.node1, sim.node2))
  sim.matrix[cbind(as.character(sim.text[,1]),as.character(sim.text[,2]))] <- sim.text[,3]
  sim.matrix
}

