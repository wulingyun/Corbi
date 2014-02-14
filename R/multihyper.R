#' The Multivariate Hypergeometric Distribution
#'
#' Generate random variables for the multivariate hypergeometric distribution
#'
#' This function generates random variables for the multivariate hypergeometric
#' distribution by iteratively calling hypergeometric random variable generator
#' \code{\link{rhyper}}.
#' 
#' @param n The number of observations
#' @param m The vector containing the number of balls in each urn
#' @param k The number of balls drawn from the urns
#' @return This function will return a vector of the same length as \code{m},
#' containing the number of balls drawn from each urn
#' 
#' @seealso \code{\link{rhyper}}
#'
#' @export
rmultihyper <- function(n, m, k)
{
  .Call(RMultiHyper, n, m, k)
}
