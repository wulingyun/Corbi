#' Corbi - Collection of reliable biomedical informatics tools
#' 
#' This pakcage provides a bundle of analysis tools for biomedical research.
#'
#' This ia a collection of biomedical informatics tools developed by WuLab at AMSS, including:
#' 
#' Network comparison:
#' \itemize{
#'   \item \code{\link{net_query}} Network querying method based on conditional random fields
#'   \item \code{\link{net_query_batch}} Batch processing version of \code{\link{net_query}}
#'   \item \code{\link{net_align}} Network alignment method based on conditional random fields
#' }
#' 
#' Functional enrichment analysis:
#' \itemize{
#'   \item \code{\link{neeat}} Network enhanced functional enrichment analysis tool
#' }
#' 
#' 
#' 
#' @name Corbi-package
#' @aliases Corbi-package Corbi
#' @docType package
#' @keywords package
#' @author Ling-Yun Wu \email{wulingyun@@gmail.com}
#' 
#' @references Qiang Huang, Ling-Yun Wu, and Xiang-Sun Zhang. An Efficient
#' Network Querying Method Based on Conditional Random Fields. Bioinformatics,
#' 27(22):3173-3178, 2011.  
#' @references Qiang Huang, Ling-Yun Wu, and Xiang-Sun Zhang. Corbi: A new
#' R package for biological network alignment and querying. BMC Systems Biology,
#' 7(Suppl 2):S6, 2013.
#' 
#' @import CRF
#' @useDynLib Corbi, .registration = TRUE
#' 
NULL
