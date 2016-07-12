search_net <- function(PPI_net, node_size=150, ori_name=FALSE){
	# The input variable PPI_net is the result of function read_net, which include the following items : "size"	"node"	"edge"	"adj_matrix"

	obj <- PPI_net$adj_matrix 								# The adjacent matrix of original network.
	n   <- nrow(obj)		  								# The number of network node.

	all_search <- NULL
	all_neighs <- NULL
	
	start_node <- PPI_net$node[sample(n, 1)]
	neigh_node <- names(which(obj[start_node,]==1))
	
	all_search <- c(all_search, start_node)
	all_neighs <- c(all_neighs, neigh_node)
	
	for (i in 1:(node_size-1)){
		num <- length(all_neighs)
		tmp <- all_neighs[sample(num, 1)]
		neigh_node <- setdiff(names(which(obj[tmp,]==1)), all_search)
	
		all_search <- c(all_search, tmp)
		all_neighs <- setdiff(all_neighs, tmp)
		all_neighs <- unique(c(all_neighs, neigh_node))
	}

	if (ori_name == FALSE){
		subinfo <- obj[all_search, all_search]
		subinfo[lower.tri(subinfo)] <- 0 
		subnetwork <- which(subinfo==1, arr.ind=T)
		rownames(subnetwork) <- NULL
	}
	if (ori_name == TRUE){
		subinfo <- obj[all_search, all_search]
		subinfo[lower.tri(subinfo)] <- 0 
		subnetwork <- which(subinfo==1, arr.ind=T)
		
		subnetwork[,1] <- as.numeric(rownames(subinfo)[subnetwork[,1]])
		subnetwork[,2] <- as.numeric(colnames(subinfo)[subnetwork[,2]])
		rownames(subnetwork) <- NULL
	}
	
	logical_result <- NULL							# Testify the fact that the edge in sub-network is unique.
	for (i in 1:(nrow(subnetwork)-1)){
		for (j in (i+1):nrow(subnetwork)){
			logical_result <- c(logical_result, setequal(subnetwork[i,], subnetwork[j,]))
		}
	}	
	return(subnetwork)
}
