library(Corbi)
library(MASS)
library(matrixcalc)
setwd("markrank")
source("search_net.R")
options(scipen=999)

sub_info <- read_net("subnetwork.txt")			# "size" "node" "adj_matrix"	
size     <- 10															# The number of pre-set differential expression node.
sam_num  <- 50															# Sample number in simulation dataset.
half <- sam_num/2

loadin <- TRUE
if (loadin == TRUE){
  load('dataset.RData')
  load('de_list.RData')
  load('answer.RData')
}else{                                      # Generate a new dataset.
  preset  <- search_net(sub_info, node_size = size, ori_name = TRUE)			# Extract a connected subnetwork from the original network. 	
  de_list <- as.character(unique(as.vector(preset)))
  l <- sub_info$size[1]											# The number of network node.
  p <- length(de_list)											# The number of pre-set differential expression node.
  
  dataset <- matrix(0, sam_num, l, dimnames = list(paste("sample", 1:sam_num, sep=""), sub_info$node))	
  vars  <- 1
  sigma <- matrix(0)
  while(!is.positive.definite(sigma)){										# The covariance matrix should be a positive definite matrix.
    vars  <- vars + 1 
    sigma <- sub_info$adj_matrix											# Construct the correlation matrix of the distribution.
    sigma[which(sigma == 1)] <- rnorm(length(which(sigma == 1)), 4, 1)		# The correlation matrix have a higher correlation within adjacent nodes.
    sigma[which(sigma == 0)] <- rnorm(length(which(sigma == 0)), 2, 1)		# The correlation matrix have a lower  correlation within non-adjacent nodes.
    diag(sigma) <- rnorm(l, vars, 1)
    sigma <- (sigma + t(sigma))/2
  }
  multi_mean <- rnorm(l, 5, 1)
  dataset <- mvrnorm(sam_num, multi_mean, sigma)								# Create simulated dataset with corresponding dimension. 
  dataset[1:half, de_list] <- dataset[1:half, de_list] * rnorm(half*p, 2, 0.1)
  
  save(dataset, file='dataset.RData')
  save(de_list, file='de_list.RData')
  answer <- NULL
}
label <- c(rep(0, half), rep(1, half))

# Prioritize disease genes using MarkRank.
adj_matrix <- as.matrix(sub_info$matrix)
adj_matrix <- adj_matrix[colnames(dataset), colnames(dataset)]
adj_matrix <- adj_matrix + t(adj_matrix)

print(system.time(result1 <- markrank(dataset, label, adj_matrix, alpha=0.8, lambda=0.2, eps=1e-10, trace=F, d=Inf)))
if (!is.null(answer) && sum(abs(result1$score - answer$result1$score)) > 1e-10){
  stop('Computation error!!')
}
s1 <- sort(result1$score, decreasing=TRUE)
print('The score of pre-set differential expression genes.')
print(result1$score[de_list])				      										
print("False discovery genes are")
print(setdiff(names(s1[1:10]), de_list))												

# Set different d for simplifying G_2 computation.
d <- 2
print(system.time(result2 <- markrank(dataset, label, adj_matrix, alpha=0.8, lambda=0.2, eps=1e-10, trace=F, d=d)))
if (!is.null(answer) && sum(abs(result2$score - answer$result2$score)) > 1e-10){
  stop('Computation error!!')
}
s2 <- sort(result2$score, decreasing=TRUE)
print('The score of pre-set differential expression genes.')
print(result2$score[de_list])				      										
print("False discovery genes are")
print(setdiff(names(s2[1:10]), de_list))												
