#' @rdname net_query
#' @name net_query
#' @aliases net_query_batch
#' 
#' @param query.nets The vector of input file names of the query networks.
#' 
#' @examples
#' 
#' \dontrun{
#' ## Batch example
#' net_query_batch(c("querynet.txt", "querynet2.txt"),
#'   "targetnet.txt", "nodesim.txt", query.type=3)
#' }
#' 
#' @export
net_query_batch <- function(query.nets, target.net, node.sim, query.type=4, delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1, output="result.txt")
{
  query.type <- as.numeric(query.type)
  delta <- lapply(list(d=delta.d, c=delta.c, e=delta.e, s=delta.s), as.numeric)
  target <- read_net(target.net)
  target$sim <- read_sim(node.sim)
  target$dist <- get_shortest_distances(target$matrix)
  
  for (query.net in query.nets)
  {
    query <- read_net(query.net)
    label <- simplify_target(query, target, delta)
    model <- build_model(query, label, delta)
    result <- solve_crf(model, query.type)
    write_result(query, label, model, result, paste(query.net, output, sep="_"))
  }
}
