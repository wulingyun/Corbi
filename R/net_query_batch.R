#' @rdname net_query
#' @name net_query
#' @aliases net_query_batch
#' 
#' @param query.nets The vector of input file names of the query networks.
#' 
#' @export
#' 
net_query_batch <- function(query.nets, target.net, node.sim, query.type=4, delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1, output="result.txt")
{
  query.type <- as.numeric(query.type)
  delta <- lapply(list(d=delta.d, c=delta.c, e=delta.e, s=delta.s), as.numeric)
  target <- read.net(target.net)
  target$sim <- read.sim(node.sim)
  target$dist <- .Call(NQ_ShortestDistances, target$matrix, rep(T, target$size))
  
  for (query.net in query.nets)
  {
    query <- read.net(query.net)
    label <- simplify.target(query, target, delta)
    model <- build.model(query, label, delta)
    result <- solve.crf(model, query.type)
    write.result(query, label, model, result, paste(query.net, output, sep="_"))
  }
}
