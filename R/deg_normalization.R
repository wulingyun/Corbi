


#' @export
get_ratio_variance <- function(expr.matrix, log.expr = FALSE)
{
  if (!log.expr) expr.matrix <- log(expr.matrix)
  .Call(ND_RatioVariance, expr.matrix)
}
