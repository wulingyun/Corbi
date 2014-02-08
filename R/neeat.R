#' Network enhanced enrichment analysis tool
#' 
#' Functional enrichment anslysis of gene sets by integrating network information
#' 
#' This is a novel functional enrichment analysis tool based on the gene set and network
#' information.
#' 
#' 
#'
#' @export
neeat <- function(w, r)
{
  .Call(NE_Depths, w, r)
}
