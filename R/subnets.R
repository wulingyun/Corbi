#' Get all sub-networks
#' 
#' @export
get_subnets <- function(net.matrix, max.size = 2)
{
  edges <- which(net.matrix != 0, arr.ind = T)
  edges <- edges[edges[,1] < edges[,2], ]
  edges <- edges[order(edges[,1]), ]
  .Call(BS_GetSubnets, edges, dim(net.matrix)[1], max.size)
}


#' Extend sub-networks by overlapping two sub-networks
#' 
#' @export
extend_subnets <- function(subnet1, subnet2, size = 0)
{
  s1 <- dim(subnet1)[2]
  s2 <- dim(subnet2)[2]
  if (s1 < 2 || s2 < 2) stop("The sizes of subnet1 and subnet2 must be at least 2!")
  min.size <- max(s1, s2) + 1
  max.size <- s1 + s2 - 1
  if (size < min.size || size > max.size) size <- min.size
  .Call(BS_ExtendSubnets, subnet1, subnet2, size)
}
