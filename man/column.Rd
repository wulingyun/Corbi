% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{column}
\alias{column}
\title{Extract a column from a matrix}
\usage{
column(m, i)
}
\arguments{
\item{m}{The matrix}

\item{i}{The column index}
}
\value{
This function will return the specified column as a vector of corresponding type.
}
\description{
Extract a specified column from a sparse matrix rapidly
}
\details{
This function implements faster column extraction algorithm for the 
\code{\link[=CsparseMatrix-class]{CsparseMatrix}} class in the package \pkg{Matrix}.
}
