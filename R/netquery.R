#' Network querying method based on conditional random fields
#' 
#' Find the best matching subnetworks from a large target network for small
#' query networks based on the conditional random fields (CRF) model.
#' 
#' This is an approach for network querying problem based on conditional random
#' field (CRF) model which can handle both undirected and directed networks,
#' acyclic and cyclic networks, and any number of insertions/deletions.
#' 
#' When querying several networks in the same target network,
#' \code{\link{net.query.batch}} will save much time.
#' 
#' \itemize{ \item query.net: The query network file is written as follows:\cr
#' v1 v2 v3 v4 v5\cr v3 v4 \cr ...  \cr where v1, v2, v3, v4, v5 ... are the
#' nodes' names and each line indicates there are edges between the first node
#' and other nodes in the line. For example, the first line denotes 4 edges:
#' (v1, v2), (v1, v3), (v1, v4), and (v1, v5).
#' 
#' \item target.net: The format of this file is the same as the query network
#' file.
#' 
#' \item node.sim: This similarity file's format is as follows:\cr v1 V1 s1 \cr
#' v1 V2 s2 \cr ...  \cr v1 is the node from the query network, V1 is the node
#' from the target network, s1 is the similarity score between the node v1 and
#' V1, and so on.
#' 
#' \item query.type: If query.type = 1, the loopy belief propagation (LBP)
#' algorithm will be applied, which is an approximate algorithm for a general
#' graph with loops. If the query is a chain or tree, there are exact
#' algorithms. Set query.type = 2 when the query is a chain, and query.type = 3
#' when the query is a tree. The heuristic algorithm will be used when
#' query.type = 4, which will try the exact algorithm (junction tree algorithm)
#' first and resort to LBP algorithm when the exact algorithm failed. The
#' default value is 4.
#' 
#' \item delta.d: The smaller delta.d is, the heavier penalty for deletions.
#' 
#' \item delta.c: The smaller delta.c is, the heavier penalty for consecutive
#' deletions.
#' 
#' \item delta.e: The smaller delta.e is, the heavier penalty for single
#' deletion.
#' 
#' \item delta.s: The larger delta.s indicates heavier penalty for insertions.
#' 
#' }
#' 
#' @aliases net.query net.query.batch
#' @param query.net The input file name of the query network.
#' @param query.nets The vector of input file names of the query networks.
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
#' @param output The suffix of output file name.
#' @references Qiang Huang, Ling-Yun Wu, and Xiang-Sun Zhang. An Efficient
#' Network Querying Method Based on Conditional Random Fields. Bioinformatics,
#' 27(22):3173–3178, 2011.
#' @examples
#' 
#' \dontrun{
#' library(Corbi)
#' 
#' ## An example: "querynet.txt", "targetnet.txt", "nodesim.txt" are
#' ## three input files in the working directory
#' net.query("querynet.txt", "targetnet.txt", "nodesim.txt", query.type=3)
#' 
#' ## Batch example
#' net.query.batch(c("querynet.txt", "querynet2.txt"),
#' 	"targetnet.txt", "nodesim.txt", query.type=3)
#' }
#' 
#' @export net.query
net.query <- function(query.net, target.net, node.sim, query.type=4, delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1, output="result.txt")
{
# options
#   query.net: input file name of query network
#   target.net: input file name of target network
#   node.sim: input file name of node similarity
#   query.type: the type of query network, 1 - general, 2 - chain, 3 - tree, 4 - heuristic, 5 - ilp
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
	target <- read.net(target.net)
	target$sim <- read.sim(node.sim)

	query <- read.net(query.net)

# compute the shortest distance matrix for the target network
# and simplify the target network

	label <- simplify.target(query, target, delta)

# build and solve CRF model

	model <- build.model(query, label, delta)
	result <- solve.crf(model, query.type)

# write result to output file

	write.result(query, label, model, result, paste(query.net, output, sep="_"))
}

net.query.batch <- function(query.nets, target.net, node.sim, query.type=4, delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1, output="result.txt")
{
	query.type <- as.numeric(query.type)
	delta <- lapply(list(d=delta.d, c=delta.c, e=delta.e, s=delta.s), as.numeric)
	target <- read.net(target.net)
	target$sim <- read.sim(node.sim)
	target$dist <- .Call(NQ_ShortestDistances, target$matrix, rep(T, target$size))

	for (query.net in query.nets)
	{
		query <- read.net(query.net)
		label <- simplify.target(query, target, delta)
		model <- build.model(query, label, delta)
		result <- solve.crf(model, query.type)
		write.result(query, label, model, result, paste(query.net, output, sep="_"))
	}
}

read.net <- function(file)
{
	net.text <- as.matrix(read.table(file, fill=T, as.is=T, col.names=1:max(count.fields(file))))
	net.node <- unique(as.character(net.text))
	net.node <- net.node[net.node != ""]
	net.size <- length(net.node)
	net.edge <- cbind(as.character(net.text[,1]), as.character(net.text[,-1]))
	net.edge <- net.edge[net.edge[,2] != "", ]
	net.matrix <- matrix(0, net.size, net.size, dimnames=list(net.node, net.node))
	net.matrix[net.edge] <- 1
	list(size=net.size, node=net.node, matrix=net.matrix)
}

write.net <- function(net, file)
{
	net.edge <- which(net$matrix == 1, arr.ind=1)
	net.edge <- matrix(net$node[net.edge], ncol=2)
	write.table(net.edge, file, quote=F, row.names=F, col.names=F)
}

read.sim <- function(file)
{
	sim.text <- read.table(file, as.is=T)
	sim.node1 <- unique(as.character(sim.text[,1]))
	sim.node2 <- unique(as.character(sim.text[,2]))
	sim.size1 <- length(sim.node1)
	sim.size2 <- length(sim.node2)
	sim.matrix <- matrix(0, sim.size1, sim.size2, dimnames=list(sim.node1, sim.node2))
	sim.matrix[cbind(as.character(sim.text[,1]),as.character(sim.text[,2]))] <- sim.text[,3]
	sim.matrix
}

simplify.target <- function(query, target, delta)
{
	net.sim <- matrix(0, query$size, target$size, dimnames=list(query$node, target$node))
	n1 <- query$node[query$node %in% rownames(target$sim)]
	n2 <- target$node[target$node %in% colnames(target$sim)]
	net.sim[n1, n2] <- target$sim[n1, n2]
	select <- colSums(net.sim > delta$d) > 0
	net.size <- sum(select)
	if (net.size <= 0)
	{
		stop("The simplified target network is empty, you should give a smaller cut values!")
	}
	net.node <- target$node[select]
	net.matrix <- target$matrix[select, select]
	net.sim <- net.sim[, select]
	net.sim[net.sim <= delta$d] <- 0
	if (is.null(target$dist))
	{
		net.dist <- .Call(NQ_ShortestDistances, target$matrix, select)[select, select]
	}
	else
	{
		net.dist <- target$dist[select, select]
	}
	net.dist[net.dist == -1] <- 10000
	net.dist[cbind(1:net.size, 1:net.size)] <- 1
	list(size=net.size, node=net.node, matrix=net.matrix, sim=net.sim, dist=net.dist)
}

build.model <- function(query, label, delta)
{
	query.size <- query$size
	query.net <- query$matrix
	query.net[query.net != 0] <- 1
	label.size <- label$size
	label.gap <- label.size + 1
	n.labels <- label.gap
	label.weight <- 1 / label$dist ^ delta$s

	S <- cbind(label$sim, delta$d)

	W <- matrix(nrow=n.labels, ncol=n.labels)
	W[1:label.size, 1:label.size] <- label.weight
	W[label.gap,] <- delta$e
	W[,label.gap] <- delta$e
	W[label.gap, label.gap] <- delta$c

	crf.net <- query.net + t(query.net)
	crf.net[crf.net != 0] <- 1
	crf <- make.crf(crf.net, rowSums(S > 0))

	crf$state.map <- matrix(label.gap, nrow=crf$n.nodes, ncol=crf$max.state)
	for (i in 1:crf$n.nodes)
	{
		crf$state.map[i, 1:crf$n.states[i]] <- which(S[i,] > 0)
		crf$node.pot[i,] <- S[i, crf$state.map[i,]]
	}

	for (e in 1:crf$n.edges)
	{
		n1 <- crf$edges[e, 1]
		n2 <- crf$edges[e, 2]
		m1 <- 1:crf$n.states[n1]
		m2 <- 1:crf$n.states[n2]
		S1 <- matrix(crf$node.pot[n1, m1], crf$n.states[n1], crf$n.states[n2])
		S2 <- matrix(crf$node.pot[n2, m2], crf$n.states[n1], crf$n.states[n2], byrow=T)
		W1 <- W[crf$state.map[n1, m1], crf$state.map[n2, m2]]
		W2 <- W[crf$state.map[n2, m2], crf$state.map[n1, m1]]
		crf$edge.pot[[e]] <- (S1 + S2) * pmax(W1 * query.net[n1, n2], t(W2) * query.net[n2, n1]) / 2
	}
	crf
}

decode.heuristic <- function(crf)
{
	result <- try(decode.junction(crf), T)
	if (class(result) == "try-error")
	{
		result <- decode.lbp(crf)
	}
	result
}

solve.crf <- function(model, query.type)
{
	decode.method <- list(decode.lbp, decode.chain, decode.tree, decode.heuristic, decode.ilp)
	if (!is.numeric(query.type) || query.type > length(decode.method) || query.type < 1) query.type = 4
	result <- decode.method[[query.type]](model)
	result <- model$state.map[cbind(1:model$n.nodes, result)]
}

write.result <- function(query, label, model, result, filename="result.txt")
{
	query.name <- query$node
	label.name <- c(label$node, "gap")
	label.dist <- matrix("gap", nrow=label$size+1, ncol=label$size+1)
	label.dist[1:label$size, 1:label$size] <- label$dist

	con <- file(as.character(filename), "w")
	writeLines("node match:", con, sep="\n")
	for (i in 1:query$size)
	{
		writeLines(paste(query.name[i], "  ", label.name[result[i]]), con, sep="\n");
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
}
