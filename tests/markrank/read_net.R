read_net <- function(filename, directed=FALSE, usename=FALSE){
	# This function can read the edge file (input file of cytoscape) and obtain the related information.
	# Note that the original edge file should in matrix form. e.g. the edge are derived from "which" function. Therefore the original file is numeric (no quote).
	
	if (usename == FALSE){
		if (class(filename) == "character"){
			net_text <- as.matrix(read.table(filename, fill=T, as.is=T, col.names=1:max(count.fields(filename))))
		}
		if (class(filename) == "matrix"){
			net_text <- filename
		}
		net_node <- as.character(sort(unique(c(net_text[,1], net_text[,2])))) 		# The name list of network node, which can mapping to the original node name list.
		net_size <- length(net_node)												# The number of network node.
		
		# The adjacent matrix of the object network, which is not a symmetric matrix if the network is a directed graph and vice versa.
		net_matrix <- matrix(0, net_size, net_size, dimnames=list(net_node, net_node))
		if (directed == FALSE){
			net_matrix[net_text] <- 1	
			net_matrix <- net_matrix + t(net_matrix)
		}else{
			net_matrix[net_text] <- 1		
		}
		list(size=c(net_size,nrow(net_text)), node=net_node, edge=net_text, adj_matrix=net_matrix)
	}else{
		if (class(filename) == "character"){
			net_text <- as.matrix(read.table(filename, fill=T, as.is=T, col.names=1:max(count.fields(filename))))
		}
		if (class(filename) == "matrix"){
			net_text <- filename
		}
		net_text[,1] <- as.character(net_text[,1])
		net_node <- unique(c(net_text[,1], net_text[,2]))							# The name list of network node, which can mapping to the original node name list.
		net_size <- length(net_node)												# The number of network node.
		
		# The adjacent matrix of the object network, which is not a symmetric matrix if the network is a directed graph and vice versa.
		net_matrix <- matrix(0, net_size, net_size, dimnames=list(net_node, net_node))
		if (directed == FALSE){
			net_matrix[net_text] <- 1	
			net_matrix <- net_matrix + t(net_matrix)
		}else{
			net_matrix[net_text] <- 1		
		}
		list(size=c(net_size, sum(net_matrix)/2), node=net_node, edge=net_text, adj_matrix=net_matrix)
	}
}
