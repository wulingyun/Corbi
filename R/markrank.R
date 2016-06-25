#' @import Matrix
#' 
#' @export
markrank <- function(dataset, label, PPI_adj, alpha, lambda, eps=1e-10, E_value=NULL, trace=TRUE)
{
  # INPUT DESCRIPTION: 
  # dataset : The microarray expression matrix of related disease. Each row represents a sample and each column represents a gene.
  # label   : The 0-1 binary phenotype vector of dataset samples. The size of label must accordance with the sample number in dataset.
  # PPI_adj : The 0-1 binary adjacent matrix of a connected biological network. Here the node set should be the same order as the gene set in expression matrix. 
  
  # The coefficient matrix A = lambda*A1 + (1-lambda)*A2. 
  # The coefficient matrix A1 and A2 are derived from NET1 and NET2 (the adjacent matrix of original PPI network and the discriminative potential network) respectively.
  # The matrix equation is : (I-alpha*A)*R = (1-alpha)*E, where the I is n*n unit matrix.
  # Solving the final matrix equation by iterative method or matrix inversion method. First one is preferred.
  
  # E_value : Contain the prior information about the importance of nodes. Default is absolute Pearson correlation coefficient (PCC).
  # alpha   : The convex combination coefficient of network effect and prior information vector E.
  # lambda  : The convex combination coefficient of two network effects.
  # eps     : The cut-off value of the iteration solving method.

  m <- nrow(dataset)												# Sample number.
  n <- ncol(dataset)												# Gene number.
  if (length(label) != m){
    stop("The size of label must accordance with the sample number in dataset.")
  }
  if (!setequal(rownames(PPI_adj), colnames(dataset))){
    stop("The mapping for genes is not performed")
  }
  NET1 <- PPI_adj
  degs <- rowSums(NET1)
  isos <- which(degs == 0)
  if (length(isos) > 0){
    stop("The input biological network should be a connected network.")
  }	
  
  D1 <- diag(1/degs)
  A1 <- t(NET1)%*%D1
  print("Discriminative potential network G2 computing...")
  D2 <- Matrix(0, n, n, sparse=TRUE, dimnames=list(colnames(dataset), colnames(dataset)))
  system.time(NET2 <- .markrank.compute_net2(dataset, label, trace=trace))
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
  tm <- 1															# Iteration steps.
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
              initial_pars = list(alpha = alpha, lambda = lambda, eps = eps)
  ))
}


#' @import Matrix mpmi
#' 
.markrank.compute_net2 <- function(dataset, label, trace=FALSE)
{
  l <- ncol(dataset)
  label <- as.matrix(as.numeric(label))
  MI1 <- mminjk(dataset, label, level=0L, na.rm=FALSE)
  MI2 <- Matrix(0, l, l, dimnames=list(colnames(dataset),colnames(dataset)), sparse=TRUE)
  for (i in 1:(l-1)){
    if (trace == TRUE && i%%10 == 0) print(i)
    dataset_tmp <- (dataset[,(i+1):l, drop=FALSE] + dataset[,i])/sqrt(2)
    MI2[(i+1):l, i] <- mminjk(dataset_tmp, label, level=0L, na.rm=FALSE)
  }
  NET2 <- MI2 + t(MI2) - as.vector(MI1)
  NET2[which(NET2 <= 0)] <- 0
  diag(NET2) <- 0
  NET2 <- Matrix(NET2, sparse=TRUE)
  
  return(NET2)
}
