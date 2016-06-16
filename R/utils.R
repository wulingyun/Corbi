#' Extract a column from a matrix
#' 
#' Extract a specified column from a sparse matrix rapidly
#' 
#' This function implements faster column extraction algorithm for the 
#' \code{\link[=CsparseMatrix-class]{CsparseMatrix}} class in the package \pkg{Matrix}.
#' 
#' @param m The matrix
#' @param i The column index
#' 
#' @return This function will return the specified column as a vector of corresponding type.
#' 
#' @import Matrix
#' 
#' @export
column <- function(m, i)
{
  if (inherits(m, "CsparseMatrix")) {
    v <- vector(typeof(m@x), m@Dim[1])
    p <- (m@p[i]+1):m@p[i+1]
    if (p[1] <= p[length(p)])
      v[m@i[p]+1] <- m@x[p]
  }
  else
    v <- m[,i]
  v
}


#' Extract a submatrix from a matrix
#' 
#' Extract a specified submatrix from a sparse matrix rapidly
#' 
#' This function implements faster submatrix extraction algorithm for the 
#' \code{\link[=CsparseMatrix-class]{CsparseMatrix}} class in the package \pkg{Matrix}.
#' 
#' @param m The matrix
#' @param rows The integer vectors of row index(es)
#' @param cols The integer vectors of column index(es)
#' 
#' @return This function will return the specified submatrix as a matrix of corresponding type.
#' 
#' @export
submatrix <- function(m, rows, cols)
{
  sapply(cols, function(i) column(m, i)[rows])
}


#' The number of non-zero values of a submatrix
#' 
#' Retuen the number of non-zero values of the specified submatrix of a given sparse matrix rapidly
#' 
#' This function implements faster calculation algorithm for the 
#' \code{\link[=CsparseMatrix-class]{CsparseMatrix}} and \code{\link[=RsparseMatrix-class]{RsparseMatrix}}
#' class in the package \pkg{Matrix}.
#' 
#' @param m The matrix
#' @param rows The integer vector of row index(es) or logical vector indicated the selected rows
#' @param cols The integer vector of column index(es) or logical vector indicated the selected cols
#' 
#' @return This function will return the number of non-zero values in the specified submatrix.
#' 
#' @export
nnzero <- function(m, rows = 1:dim(m)[1], cols = 1:dim(m)[2])
{
  r <- logical(dim(m)[1])
  c <- logical(dim(m)[2])
  if (is.integer(rows)) {
    r[rows] <- TRUE
  }
  else {
    r <- rows
    rows <- which(r)
  }
  if (is.integer(cols)) {
    c[cols] <- TRUE
  }
  else {
    c <- cols
    cols <- which(c)
  }
  fun <- function(i) if (m@p[i] < m@p[i+1]) (m@p[i]+1):m@p[i+1] else NULL
  if (sum(r) == 0 || sum(c) == 0)
    0
  else if (inherits(m, "CsparseMatrix")) {
    p <- unlist(lapply(cols, fun))
    sum(m@x[p[r[m@i[p]+1]]] != 0)
  }
  else if (inherits(m, "RsparseMatrix")) {
    p <- unlist(lapply(rows, fun))
    sum(m@x[p[c[m@i[p]+1]]] != 0)
  }
  else
    Matrix::nnzero(m[rows, cols])
}


#' Cohen's kappa score
#' 
#' Calculate Cohen's kappa score for two vectors.
#' 
#' This function calculate Cohen's kappa score for two logical vectors.
#' 
#' @param x1 The first logical vector
#' @param x2 The second logical vector
#' 
#' @return The Cohen's kappa score
#' 
#' @export
kappa_score <- function(x1, x2)
{
  if (length(x1) != length(x2)) stop("x1 and x2 are not of same length!")
  t <- length(x1)
  p1 <- sum(x1)
  p2 <- sum(x2)
  p <- sum(x1 == x2) / t
  e <- (p1*p2 + (t-p1)*(t-p2)) / t^2
  (p-e) / (1-e)
}
