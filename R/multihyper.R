#' The Multivariate Hypergeometric Distribution
#'
#' Generate random variables for the multivariate hypergeometric distribution
#'
#' This function generates random variables for the multivariate hypergeometric
#' distribution by iteratively calling hypergeometric random variable generator
#' \code{\link{rhyper}}.
#' 
#' @param n The number of observations.
#' @param k The total number of balls drawn from the urn.
#' @param m The integer vector containing the number of balls of each color in the urn.
#' Length of vector is the number of colors.
#' @return This function will return a matrix of \code{length(m)} rows and \code{n} columns,
#' and each column contains the number of balls of each color drawn from the urn.
#' 
#' @seealso \code{\link{rhyper}}
#'
#' @export
rmultihyper <- function(n, k, m)
  .Call(RMultiHyper, n, k, m)

pmultihyper <- function(x, k, m, w)
  .Call(PMultiHyper, x, k, m, w)

pmultinom <- function(x, k, m, w)
  .Call(PMultiNom, x, k, m, w)
