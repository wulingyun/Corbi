#' The best subnetworks
#' 
#' Search best subnetworks that maximize given objective functions.
#' 
#' Enumerate and search the best subnetworks that maximize given objective function. If the size of 
#' subnetworks <= \code{exhaust.size}, exact exhaustive searching is applied,
#' otherwise, heuristic searching algorithm is used.
#' 
#' @param func The objective function to maximize
#' @param net.matrix The adjacent matrix of network
#' @param max.size The maximal size of subnetworks
#' @param exhaust.size The maximal size of subnetworks that use exhaustive searching strategy
#' @param max.top The maiximal number of top candidates kept for evaluation of next size,
#' used in heuristic searching strategy
#' 
#' @return A list with the following two components:
#'   \item{subnets}{The list of top subnetworks in different sizes}
#'   \item{obj.values}{The list of objective values of corresponding subnetworks}
#'   
#' @seealso get_subnets, extend_subnets
#' 
#' @examples
#' 
#' library(Corbi)
#' net <- matrix(FALSE, nrow=10, ncol=10)
#' net[sample.int(100, 20)] <- TRUE
#' net <- net | t(net)
#' func <- function(subnet) max(subnet) - min(subnet)
#' result <- best_subnets(func, net, 5)
#' 
#' @export
best_subnets <- function(func, net.matrix, max.size = 10, exhaust.size = 5, max.top = 10000)
{
  exhaust.size <- min(max.size, exhaust.size)
  subnets <- get_subnets(net.matrix, exhaust.size)
  obj.values <- lapply(subnets, function (x) apply(x, 1, func))
  orders <- lapply(obj.values, function(x) order(x, decreasing = TRUE))
  subnets <- lapply(1:exhaust.size, function(x) matrix(subnets[[x]][orders[[x]], ], ncol = x))
  obj.values <- lapply(1:exhaust.size, function(x) obj.values[[x]][orders[[x]]])
  i <- exhaust.size + 1
  while (i <= max.size)
  {
    sub1 <- subnets[[i-1]]
    top <- min(max.top, dim(sub1)[1])
    sub1 <- matrix(sub1[1:top, ], ncol = i-1)
    if (top > 1) sub2 <- sub1
    else sub2 <- subnets[[2]]
    subnets[[i]] <- extend_subnets(sub1, sub2, i)
    obj.values[[i]] <- apply(subnets[[i]], 1, func)
    orders <- order(obj.values[[i]], decreasing = TRUE)
    subnets[[i]] <- matrix(subnets[[i]][orders, ], ncol = i)
    obj.values[[i]] <- obj.values[[i]][orders]    
    i <- i + 1
  }
  list(subnets=subnets, obj.values=obj.values)
}


#' All subnetworks of limited size
#' 
#' Enumerate all subnetworks of size <= \code{max.size} from given network.
#' 
#' @param net.matrix The adjacent matrix of network
#' @param max.size The maximal size of subnetworks
#' 
#' @return A list of generated subnetworks, with element $i$ corresponds the subnetworks
#' of size $i$. Each element is a matrix, in which each row represents a subnetwork.
#' 
#' @examples
#' 
#' library(Corbi)
#' net <- matrix(FALSE, nrow=10, ncol=10)
#' net[sample.int(100, 20)] <- TRUE
#' net <- net | t(net)
#' subnets <- get_subnets(net, 3)
#' 
#' @export
get_subnets <- function(net.matrix, max.size = 2)
{
  edges <- which(net.matrix != 0, arr.ind = TRUE)
  edges <- edges[edges[,1] < edges[,2], ]
  edges <- edges[order(edges[,1]), ]
  .Call(BS_GetSubnets, edges, dim(net.matrix)[1], max.size)
}


#' Extend subnetworks from smaller subnetworks
#'
#' Extend subnetworks by pairwise overlapping two sets of smaller subnetworks.
#'
#' Enumerate all possible subnetworks of desired size by pairwise overlapping two sets of 
#' subnetworks of size \code{s1} and \code{s2}. The desired size should be between 
#' \code{max(s1,s2)+1} and \code{s1+s2-1}. Invalid desired size will be replaced by the 
#' minimum allowed value \code{max(s1,s2)+1}.
#' 
#' @param subnet1 The matrix representing the first set of subnetworks
#' @param subnet2 The matrix representing the second set of subnetworks
#' @param size The desired size of extended subnetworks
#' 
#' @return A matrix represents the extended subnetworks, in which each row represents a subnetwork.
#' 
#' @examples
#' 
#' library(Corbi)
#' net <- matrix(FALSE, nrow=10, ncol=10)
#' net[sample.int(100, 20)] <- TRUE
#' net <- net | t(net)
#' subnets <- get_subnets(net, 3)
#' subnets[[4]] <- extend_subnets(subnets[[3]], subnets[[2]], 4)
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
