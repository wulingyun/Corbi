#' The Multivariate Hypergeometric Distribution
#'
#' Generate random variables for the multivariate hypergeometric distribution
#'
#' This function generates random variables for the multivariate hypergeometric
#' distribution by iteratively calling hypergeometric random variable generator
#' \code{\link{rhyper}}.
#' 
#' @param n The number of observations
#' @param k The number of balls drawn from the urns
#' @param m The integer vector containing the number of balls in each urn
#' @return This function will return a matrix of \code{length(m)} rows and \code{n} columns,
#' containing the number of balls drawn from each urn
#' 
#' @seealso \code{\link{rhyper}}
#'
#' @export
rmultihyper <- function(n, k, m)
{
  .Call(RMultiHyper, n, k, m)
}
