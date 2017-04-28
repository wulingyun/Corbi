#' MarkRank
#' 
#' MarkRank is a novel proposed network-based model, which can identify the cooperative 
#' biomarkers for heterogeneous complex disease diagnoses.
#' 
#' MarkRank is a network-based biomarker identification method to prioritize disease genes 
#' by integrating multi-source information including the biological network, e.g 
#' protein-protein interaction (PPI) network, the prior information about related diseases,
#' and the discriminative power of cooperative gene combinations. MarkRank shows that
#' explicit modeling of gene cooperative effects can greatly improve biomarker identification
#' for complex disease, especially for diseases with high heterogeneity.
#' 
#' MarkRank algorithm contains mainly two steps: 1) The construction of gene cooperation network G_2 
#' and 2) a random walk based iteration procedure. The following descriptions will help the users to
#' using \code{markrank} more convenient:
#' 
#' 1) As for the construction of the gene cooperation network, we  suggest the user to set
#' \code{trace=TRUE} to output the G_2 computation process. The G_2 construction step finished
#' if the output number is identical to the gene number of the input expression matrix. The parameter \code{d}
#' introduced the structure information of used biological network to facilitate the construction
#' of G_2, only the gene pairs whose shortest distances in network are less than \code{d} participate
#' the G_2 computation. We suggest \code{d=Inf}, the default value, to fully use the information of expression
#' matrix. If the user given a preset \code{d}, the distance matrix of input network \code{dis} 
#' will be returned.
#' 
#' 2) MarkRank uses a random-walk based iteration procedure to score each gene. The detailed formula is: 
#' 
#' \code{score} = \code{alpha}*[\code{lambda}*A1 + (1-\code{lambda})*A2]*\code{score} + (1-\code{alpha})*\code{E_value}.
#' 
#' The users could set an appropriate parameter settings in their pracitical application.
#' Our suggested value is \code{alpha}=0.8 and \code{lambda}=0.2. The model input parameter combinations and iteration steps will
#' be returned in output components \code{initial_pars} and \code{steps}, respectively. Because the iteration step is separate with
#' the cooperation network construction, the user can use the parameter \code{Given_NET2} to tune
#' the model parameters. In detail, the user could set 
#' 
#' \code{Given_NET2 = result$NET2} 
#' 
#' in \code{markrank} input to avoid the repeated computation of G_2, where the object \code{result}
#' is the returned variable of \code{markrank} function.
#' 
#' 3) The final MarkRank score for each gene is in output \code{score}. The users could sort
#' this result and use the top ranked genes for further analysis.
#' 
#' 
#' @param dataset The microarray expression matrix of related disease. Each row represents
#' a sample and each column represents a gene.
#' @param label The 0-1 binary phenotype vector of dataset samples. The size of label must
#' accord with the sample number in dataset.
#' @param adj_matrix The 0-1 binary adjacent matrix of a connected biological network. 
#' Here the node set should be the same order as the gene set in expression matrix. 
#' @param alpha The convex combination coefficient of network effect and prior information vector \code{E_value}.
#' The range of alpha is in \code{[0,1]}. A larger alpha will lay more emphasis on the 
#' network information. The default value is 0.8.
#' @param lambda In the random walk-based iteration, matrix A1 reflects the stucture information of the 
#' biological network, whereas matrix A2 reflects the cooperative effect of gene combinations.
#' Parameter lambda is the convex combination coefficient of two network effects. The range of lambda is
#' in \code{[0,1]}. A larger lambda will lay more emphasis on the A1. The default value is 0.2.
#' @param eps The stop criteria for the iterative solution method. The default value is 1e-10.
#' @param E_value A vector containing the prior information about the importance of nodes. Default is the 
#' absolute Pearson correlation coefficient (PCC).
#' @param trace Locaical variable indicated whether tracing information on the progress of the gene cooperation
#' network construction is produced.
#' @param d Threshold for simplifying the G_2 computation. Only the gene pairs whose shortest distances in PPI network are 
#' less than d participate in the G_2 computation. The default value is Inf.
#' @param Given_NET2 Whether a computed cooperation network is given for tuning parameter. See Details
#' for a more specific description.
#' 
#' @return This function will return a list with the following components:
#'   \item{score}{The vector of final MarkRank scores for each gene.}
#'   \item{steps}{The final iteration steps in random walk based scoring procedure.}
#'   \item{NET2}{The weighted adjacent matrix of gene cooperation network.}
#'   \item{initial_pars}{The initial/input parameter values used in MarkRank.}
#'   \item{dis}{The pairwise distance matrix of input network. This variable will be \code{Null} if input d=Inf.}
#' 
#' @references Duanchen Sun, Xianwen Ren, Eszter Ari, Tamas Korcsmaros, Peter Csermely,
#' Ling-Yun Wu. Discovering cooperative biomarkers for heterogeneous complex disease diagnoses.
#' Manuscript, 2017.
#' 
#' @import Matrix
#' 
#' @export
markrank <- function(dataset, label, adj_matrix, alpha=0.8, lambda=0.2, eps=1e-10, E_value=NULL, trace=TRUE, d=Inf, Given_NET2=NULL)
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
  if (class(Given_NET2) == "NULL"){
  	system.time(NET2 <- .markrank.compute_net2(dataset, label, dis, d, trace=trace))
  }
  D2 <- Matrix(0, n, n, sparse=TRUE, dimnames=list(colnames(dataset), colnames(dataset)))
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



