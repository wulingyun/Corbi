#' Find the optimal sub-networks
#' 
#' @export
best_subnets <- function(func, net.matrix, max.size = 10, exhaust.size = 5, max.top = 10000)
{
  exhaust.size <- min(max.size, exhaust.size)
  subnets <- get_subnets(net.matrix, exhaust.size)
  results <- lapply(subnets, function (x) apply(x, 1, func))
  orders <- lapply(results, function(x) order(x, decreasing = TRUE))
  subnets <- lapply(1:exhaust.size, function(x) matrix(subnets[[x]][orders[[x]], ], ncol = x))
  results <- lapply(1:exhaust.size, function(x) results[[x]][orders[[x]]])
  i <- exhaust.size + 1
  while (i <= max.size)
  {
    sub1 <- subnets[[i-1]]
    top <- min(max.top, dim(sub1)[1])
    sub1 <- matrix(sub1[1:top, ], ncol = i-1)
    if (top > 1) sub2 <- sub1
    else sub2 <- subnets[[2]]
    subnets[[i]] <- extend_subnets(sub1, sub2, i)
    results[[i]] <- apply(subnets[[i]], 1, func)
    orders <- order(results[[i]], decreasing = TRUE)
    subnets[[i]] <- matrix(subnets[[i]][orders, ], ncol = i)
    results[[i]] <- results[[i]][orders]    
    i <- i + 1
  }
  list(subnets=subnets, results=results)
}


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
