#' MarkRank
#' 
#' Prioritization of network biomarkers for complex diseases
#' 
#' MarkRank is a network biomarker identification method to prioritize disease genes 
#' by integrating multi-source information including the biological network, eg. 
#' protein-protein interaction (PPI) network, the prior information about related diseases,
#' and the discriminative power of cooperative gene combinations.
#' 
#' The random-walk based iteration formula is: 
#' 
#' R = alpha*[lambda*A1 + (1-lambda)*A2]*R + (1-alpha)*E_value.
#' 
#' 
#' @param dataset The microarray expression matrix of related disease. Each row represents
#' a sample and each column represents a gene.
#' @param label The 0-1 binary phenotype vector of dataset samples. The size of label must
#' accord with the sample number in dataset.
#' @param adj_matrix The 0-1 binary adjacent matrix of a connected biological network. 
#' Here the node set should be the same order as the gene set in expression matrix. 
#' @param alpha The convex combination coefficient of network effect and prior information vector E.
#' @param lambda The convex combination coefficient of two network effects.
#' @param eps The stop criteria for the iterative solution method.
#' @param E_value A vector containing the prior information about the importance of nodes. Default is 
#' absolute Pearson correlation coefficient (PCC).
#' @param trace Locaical variable indicated whether tracing information on the progress of the iterative
#' solution method is produced.
#' @param d Threshold for simplifying G_2 computation. Only the gene pairs whose shortest distances in PPI network are 
#' less than d participate in G_2 computation. Default is Inf.
#'
#' 
#' @return This function will return a list with components:
#'   \item{score}{The vector of final MarkRank scores for each gene.}
#'   \item{steps}{The final steps used by iterative method.}
#'   \item{NET2}{The weighted adjacent matrix of discriminative potential network.}
#'   \item{initial_pars}{The initial parameter values used in computation.}
#'   \item{dis}{The pairwise distance matrix of input network. Null if d=Inf}
#' 
#' @references Duanchen Sun, Xianwen Ren, Eszter Ari, Tamas Korcsmaros, Peter Csermely,
#' Ling-Yun Wu. Prioritization of network biomarkers for complex diseases via MarkRank.
#' Manuscript, 2016.
#' 
#' @import Matrix
#' 
#' @export
markrank <- function(dataset, label, adj_matrix, alpha, lambda, eps=1e-10, E_value=NULL, trace=TRUE, d=Inf)
{
  
  m <- nrow(dataset)											# Sample number.
  n <- ncol(dataset)											# Gene number.
  if (length(label) != m) {
    stop("The size of label must accord with the sample number in dataset.")
  }
  if (!setequal(rownames(adj_matrix), colnames(dataset))) {
    stop("The mapping for genes is not performed.")
  }
  if (!all(rownames(adj_matrix) == colnames(dataset))) {
    stop("The order of genes is not coincident in expression matrix and adjacent matrix.")
  }
  NET1 <- adj_matrix
  degs <- rowSums(NET1)
  if (sum(degs == 0) > 0) {
    stop("The input biological network should be a connected network.")
  }	

# The coefficient matrix A = lambda*A1 + (1-lambda)*A2.
# The coefficient matrix A1 and A2 are derived from NET1 and NET2 (the adjacent matrix of 
# original biological network and the discriminative potential network) respectively.
# The matrix equation is : (I-alpha*A)*R = (1-alpha)*E, where the I is n*n unit matrix.
# Solving the final matrix equation by iterative method or matrix inversion method. 
# First one is preferred.
  
  D1 <- diag(1/degs)
  A1 <- t(NET1)%*%D1
 
  if (d != Inf){
    dis <- get_shortest_distances(adj_matrix)
  }else{
    dis <- NULL
  }
  if (trace) print("Computing discriminative potential network ...")
  D2 <- Matrix(0, n, n, sparse=TRUE, dimnames=list(colnames(dataset), colnames(dataset)))
  system.time(NET2 <- .markrank.compute_net2(dataset, label, dis, d, trace=trace))
  diag(D2) <- 1/rowSums(NET2)
  A2 <- t(NET2)%*%D2
  A  <- lambda*A1 + (1-lambda)*A2
  
  if (class(E_value) == "NULL"){
	  PCC <- NULL														
      for (i in 1:n){
        PCC[i] <- stats::cor(dataset[,i], label, method="pearson")
      }
    E_value <- abs(PCC)											# Use the absolute value of Pearson correlation coefficient(PCC) as the prior information.
  }
  R1 <- as.matrix(stats::runif(n))
  R2 <- as.matrix(stats::runif(n))
  tm <- 1														# Iteration steps.
  while (sum(abs(R1-R2)) >= eps){
    tm <- tm + 1
    R1 <- R2
    R2 <- alpha*A%*%R1 + (1-alpha)*E_value
  }
  R <- as.vector(R2/sum(R2))
  names(R) <- colnames(dataset)
  
  return(list(score = R,
              steps = tm,
              NET2  = NET2,
              initial_pars = list(alpha = alpha, lambda = lambda, eps = eps),
		  	      dis   = dis
  ))
}


#' @import Matrix mpmi
#' 
.markrank.compute_net2 <- function(dataset, label, dis=NULL, d=Inf, trace=FALSE)
{
  l <- ncol(dataset)
  label <- as.matrix(as.numeric(label))
  MI1 <- mminjk(dataset, label, level=0L, na.rm=FALSE)
  MI2 <- Matrix(0, l, l, dimnames=list(colnames(dataset),colnames(dataset)), sparse=TRUE)
  
  
  if (d == Inf){
	  for (i in 1:(l-1)) {
		  if (trace == TRUE && i%%10 == 0) print(i)
		  dataset_tmp <- (dataset[,(i+1):l, drop=FALSE] + dataset[,i])/sqrt(2)
		  MI2[(i+1):l, i] <- mminjk(dataset_tmp, label, level=0L, na.rm=FALSE)
	  }			
  }else{
	  for (i in 1:(l-1)) {
		  if (trace == TRUE && i%%10 == 0) print(i)
		
		  index <- which(dis[i,(i+1):l] <= d)
		  if (length(index) > 0){
	      inds  <- ((i+1):l)[index]
		    dataset_tmp <- (dataset[,inds, drop=FALSE] + dataset[,i])/sqrt(2)
		    MI2[inds, i] <- mminjk(dataset_tmp, label, level=0L, na.rm=FALSE)
		  }
	  }
  }
  NET2 <- MI2 + t(MI2) - as.vector(MI1)
  NET2[which(NET2 <= 0)] <- 0
  diag(NET2) <- 0
  NET2 <- Matrix(NET2, sparse=TRUE)
  
  return(NET2)
}



