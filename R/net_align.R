#' Network alignment method based on conditional random fields
#' 
#' Find the maximal matching subnetworks from a target network for a query
#' network based on the conditional random fields (CRF) model.
#' 
#' This is an approach for network alignment problem based on conditional
#' random field (CRF) model which uses the node similarity and structure
#' information equally. This method is based on our network querying method
#' \code{\link{net_query}}. This method uses an iterative strategy to get the
#' one-to-one map between the query network and target netowrk.
#' 
#' More details can be seen in \code{\link{net_query}}.
#' 
#' @param query.net The input file name of the query network.
#' @param target.net The input file name of the target network.
#' @param node.sim The input file name of the node similarity scores between
#' the query network and the target network.
#' @param query.type The querying network type: 1 - general, 2 - chain, 3 -
#' tree, 4 - heuristic.
#' @param delta.d The parameter delta.d is a parameter for deletions.
#' @param delta.c The parameter delta.c is a parameter for consecutive
#' deletions.
#' @param delta.e The parameter delta.e is a parameter for single deletion.
#' @param delta.s The parameter delta.s is a parameter for insertions.
#' @param output The suffix of output file name. The output contains two files
#' in the working directory. One is the matching nodes and edges between query
#' network and target network, the other is the unique matching node pairs.
#' 
#' @references Qiang Huang, Ling-Yun Wu, and Xiang-Sun Zhang. CNetA: Network
#' alignment by combining biological and topological features. In Proceedings
#' of 2012 IEEE International Conference on Systems Biology (ISB), 220-225,
#' IEEE, 2012.
#' @references Qiang Huang, Ling-Yun Wu, and Xiang-Sun Zhang. Corbi: A new
#' R package for biological network alignment and querying. BMC Systems Biology,
#' 7(Suppl 2):S6, 2013.
#' 
#' @examples
#' 
#' \dontrun{
#' library(Corbi)
#' 
#' ## An example: "querynet.txt", "targetnet.txt", "nodesim.txt" are
#' ## three input files in the working directory
#' net_align("querynet.txt", "targetnet.txt", "nodesim.txt")
#' }
#' 
#' @import CRF
#' 
#' @export
net_align <- function(query.net, target.net, node.sim, query.type=4, delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1, output="result.txt")
{
# options
#   query.net: input file name of query network
#   target.net: input file name of target network
#   node.sim: input file name of node similarity
#   query.type: the type of query network, 1 - general, 2 - chain, 3 - tree, 4 - heuristic
#   delta: the parameters \Delta_d, \Delta_c, \Delta_e, \Delta_s
#   output: the output filename

# read the input files
#   target$node: the node names of target network
#   target$matrix: the adjacency matrix of target network
#   target$sim: the node similarity matrix
#   query$node: the node names of query network
#   query$matrix: the adjacency matrix of query network

	query.type <- as.numeric(query.type)
	delta <- lapply(list(d=delta.d, c=delta.c, e=delta.e, s=delta.s), as.numeric)
	target <- read_net(target.net)
	target$sim <- read_sim(node.sim)
	target0 <- target
 
	query <- read_net(query.net)

	# iter times
	k=0;

	while (TRUE){

	## query --> target

		# compute the shortest distance matrix for the target network
		# and simplify the target network
		# build and solve CRF model

		label <- simplify_target(query, target, delta)
		model <- build_model(query, label, delta)
		result <- solve_crf(model, query.type)
	
		# change the target subnet as query network 

		modelchange <- target_part(label,result,target,query)
		targetpart <- modelchange$target
		query$sim <- modelchange$sim
		 

	## target part --> query

		label1 <- simplify_target( targetpart, query, delta)
		model1 <- build_model(targetpart, label1, delta)
		result1 <- solve_crf(model1, query.type)

		# update the common match pairs
		hitcom <- match_com(label,result,label1,result1)
		

	## given the stop rule

		if(k==0){ 
			hitcom0 <- hitcom 
		}else{
			if(all(dim(hitcom0)==dim(hitcom))){
				if(all(hitcom0==hitcom)){
					break;
				}else{
					hitcom0 <- hitcom
				}
			}else{
				hitcom0 <- hitcom 
			}
		}

	## update
		# the node similarity matrix
		newtarget <- update_sim(hitcom,target,query)
		target$sim <- newtarget$sim

		k=k+1;

	}
	
	# write unique corresponding list
	write_result1(query, target0, label, model, result, paste(query.net, output, sep="_"))  

}

target_part <- function(label,result,target0,query){

	sub <- result<=label$size
	node <- label$node[unique(result[sub])]
	sim <- t(label$sim[,node])
	# dimnames(sim) <- list(node,rownames(label$sim))

	target <- target0
	target$size <- length(node)
	target$node <- node

	# strategy 1, use the target sub network
	target$matrix <- target0$matrix[node,node]

	list(target=target,sim=sim)
}

match_com <- function(label,result,label1,result1)
{
	nodetarget <- label$node[result[result<=label$size]]
	nodequery <- rownames(label$sim)[result<=label$size]
	label1.name <- c(label1$node,"gap")
	N <- length(nodequery)
	sub <- matrix(FALSE,N,1)
	for (i in 1:N){
		if(nodequery[i]==label1.name[result1[nodetarget[i]==rownames(label1$sim)]]){
			sub[i] <- TRUE
		}	
	}

	hitcom <- cbind(nodequery[sub],nodetarget[sub])
	hitcom
}

update_sim <- function (hitcom,target,query)
{
	# update sim
	sim <- target$sim
	sim[hitcom[,1],] <- 0
	sim[,hitcom[,2]] <- 0
 	sim[hitcom] <-1 

	list(sim=sim)
}

uni <- function(pairs)
{
	# first, remove the gap
	pairs <- pairs[pairs[,2]!="gap",]
	# the duplicated elements sub
	subdup <- duplicated(pairs[,2])
	# the duplicated elements
	dup <- unique(pairs[subdup,2])
	# the unique matching sub
	subdup <- (1 - as.character(pairs[,2]) %in% dup)>0
	# the unique mathcing pairs
	unipair <- pairs[subdup,]

	unipair
}

write_result1 <- function(query, target, label, model, result, filename="result.txt")
{
	query.name <- query$node
	label.name <- c(label$node, "gap")

	#
	select <- matrix(TRUE,target$size,1)
	net.dist <- matrix(,target$size,target$size,dimnames=dimnames(target$matrix))
	net.dist[select,select] <- get_shortest_distances(target$matrix, select)[select, select]
	net.dist[net.dist == -1] <- Inf
	net.dist[cbind(1:target$size, 1:target$size)] <- 1

	label.dist <- matrix("gap", nrow=label$size+1, ncol=label$size+1)
	label.dist[1:label$size, 1:label$size] <- net.dist[label$node,label$node]


	con <- file(as.character(filename), "w")
	writeLines("node match:", con, sep="\n")
	for (i in 1:query$size)
	{
		writeLines(paste(query.name[i], "\t", label.name[result[i]]), con, sep="\n");
	}
	writeLines("", con, sep="\n")

	direction <- c("--->", "<---", "----")
	writeLines("edge match:", con, sep="\n")
	for (i in 1:model$n.edges)
	{
		x1 <- model$edges[i,1]
		x2 <- model$edges[i,2]
		y1 <- result[x1]
		y2 <- result[x2]
		distance <- c(label.dist[y1,y2], label.dist[y2,y1], min(label.dist[y1,y2], label.dist[y2,y1]))
		d <- (query$matrix[x1,x2] > 0) + 2 * (query$matrix[x2,x1] > 0)
		writeLines(paste(query.name[x1], direction[d], query.name[x2], "\t", label.name[y1], direction[d], label.name[y2], "\t", distance[d]), con, sep="\n")
	}
	close(con)

	#==============================================================
	# write unique corresponding list
	# query --> target; target --> query
	
	pairs1 <- cbind(query.name,label.name[result])
	x1table <- uni(pairs1)
	x1 <- c("query-->target",dim(x1table)[1])
	ftable <- rbind(x1,x1table)
	utils::write.table(ftable,paste("list",filename),row.names=FALSE,col.names=FALSE,quote=FALSE)
}
