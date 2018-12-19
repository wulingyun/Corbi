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


#' The Multivariate Hypergeometric Distribution
#'
#' The distribution function for the weighted sums of multivariate hypergeometric distribution
#'
#' This function gives the distribution function for the weighted sums of multivariate hypergeometric
#' distribution by recursively calling the hypergeometric distribution density function
#' \code{\link{dhyper}}.
#' 
#' @param x The quantile of weighted sum.
#' @param k The total number of balls drawn from the urn.
#' @param m Integer non-negative vector of length N, containing the number of balls of each color in the urn.
#' N is the number of colors.
#' @param w Numeric non-negative vector of length N, specifying the weight of balls of each color.
#' 
#' @return This function will return the probablity of \eqn{P(X \le x)}.
#' 
#' @seealso \code{\link{dhyper}}
#'
#' @export
pmultihyper <- function(x, k, m, w)
  .Call(PMultiHyper, x, k, m, w)


#' The Multinomial Distribution
#'
#' The distribution function for the weighted sums of multinomial distribution
#'
#' This function gives the distribution function for the weighted sums of multinomial
#' distribution by recursively calling the binomial distribution density function
#' \code{\link{dbinom}}.
#' 
#' @param x The quantile of weighted sum.
#' @param k The total number of balls drawn from the urn.
#' @param m Numeric non-negative vector of length N, specifying the probability
#' for drawing the ball of each color; is internally normalized to sum 1. Infinite and missing values
#' are not allowed. N is the number of colors.
#' @param w Numeric non-negative vector of length N, specifying the weight of balls of each color.
#' 
#' @return This function will return the probablity of \eqn{P(X \le x)}.
#' 
#' @seealso \code{\link{dbinom}}, \code{\link{dmultinom}}, \code{\link{rmultinom}}
#'
#' @export
pmultinom <- function(x, k, m, w)
  .Call(PMultiNom, x, k, m, w)
