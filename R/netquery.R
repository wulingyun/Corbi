net.query <- function(query.net, target.net, node.sim, query.type=2, delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1, output="result.txt")
{
# options
#   query.net: input file name of query network
#   target.net: input file name of target network
#   node.sim: input file name of node similarity
#   query.type: the type of query network, 0 - chain, 1 - tree, 2 - general
#   delta: the parameters \Delta_d, \Delta_c, \Delta_e, \Delta_s
#   output: the output filename

# read the input files
#   target$node: the node names of target network
#   target$matrix: the adjacency matrix of target network
#   target$sim: the node similarity matrix
#   query$node: the node names of query network
#   query$matrix: the adjacency matrix of query network

	delta <- list(d=delta.d, c=delta.c, e=delta.e, s=delta.s)
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

net.query.batch <- function(query.nets, target.net, node.sim, query.type=2, delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1, output="result.txt")
{
	delta <- list(d=delta.d, c=delta.c, e=delta.e, s=delta.s)
	target <- read.net(target.net)
	target$sim <- read.sim(node.sim)
	target$dist <- .Call("NQ_ShortestDistances", target$matrix, rep(T, target$size))

	for (query.net in query.nets)
	{
		query <- read.net(query.net)
		label <- simplify.target(query, target, delta)
		model <- build.model(query, label, delta)
		result <- solve.crf(model, query.type)
		write.result(query, label, model, result, paste(query.net, output, sep="_"))
	}
}

read.net <- function(net)
{
	net.text <- as.matrix(read.table(net, fill=T, as.is=T, col.names=1:max(count.fields(net))))
	net.node <- unique(as.vector(net.text))
	net.node <- net.node[net.node != ""]
	net.size <- length(net.node)
	net.edge <- cbind(net.text[,1], as.vector(net.text[,-1]))
	net.edge <- net.edge[net.edge[,2] != "", ]
	net.matrix <- matrix(0, net.size, net.size, dimnames=list(net.node, net.node))
	net.matrix[net.edge] <- 1
	list(size=net.size, node=net.node, matrix=net.matrix)
}

read.sim <- function(sim)
{
	sim.text <- read.table(sim, as.is=T)
	sim.node1 <- unique(as.vector(sim.text[,1]))
	sim.node2 <- unique(as.vector(sim.text[,2]))
	sim.size1 <- length(sim.node1)
	sim.size2 <- length(sim.node2)
	sim.matrix <- matrix(0, sim.size1, sim.size2, dimnames=list(sim.node1, sim.node2))
	sim.matrix[as.matrix(sim.text[,1:2])] <- sim.text[,3]
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
		net.dist <- .Call("NQ_ShortestDistances", target$matrix, select)[select, select]
	}
	else
	{
		net.dist <- target$dist[select, select]
	}
	net.dist[net.dist == -1] <- Inf
	net.dist[cbind(1:net.size, 1:net.size)] <- 1
	list(size=net.size, node=net.node, matrix=net.matrix, sim=net.sim, dist=net.dist)
}

build.model <- function(query, label, delta)
{
	node.pot <- S <- cbind(label$sim, delta$d)
	
	query.size <- query$size
	query.net <- query$matrix
	query.net[query.net != 0] <- 1
	label.size <- label$size
	label.gap <- label.size + 1
	label.weight <- 1 / label$dist ^ delta$s

	W <- matrix(nrow=label.gap, ncol=label.gap)
	W[1:label.size, 1:label.size] <- label.weight
	W[label.gap,] <- delta$e
	W[,label.gap] <- delta$e
	W[label.gap, label.gap] <- delta$c

	crf.net <- query.net + t(query.net)
	crf.net[crf.net != 0] <- 1
	edge.size <- sum(crf.net[upper.tri(crf.net)])
	edge.pot <- array(0, dim=c(label.gap, label.gap, edge.size))

	k <- 1
	for (i in 1:(query.size-1))
	{
		for (j in (i+1):query.size)
		{
			if (crf.net[i,j])
			{
				S1 <- matrix(S[i,], label.gap, label.gap)
				S2 <- matrix(S[j,], label.gap, label.gap, byrow=T)
				edge.pot[,,k] <- (S1 + S2) * pmax(W * query.net[i,j], t(W) * query.net[j,i]) / 2
				k <- k + 1
			}
	       }
	}
	crf <- make.crf(crf.net, label.gap)
	crf$node.pot <- node.pot
	crf$edge.pot <- edge.pot
	crf
}

solve.crf <- function(model, query.type)
{
	if (query.type == 0)
	{
		result <- decode.chain(model)
	}
	if (query.type == 1)
	{
		result <- decode.tree(model)
	}
	if (query.type == 2)
	{
		result <- decode.lbp(model)
	}
	result
}

write.result <- function(query, label, model, result, filename="result.txt")
{
	query.name <- query$node
	label.name <- c(label$node, "gap")

	con <- file(as.character(filename), "w")
	writeLines("node match:", con, sep="\n")
	for (i in 1:query$size)
	{
		writeLines(paste(query.name[i], "  ", label.name[result[i]]), con, sep="\n");
	}
	writeLines("", con, sep="\n")

	writeLines("edge match:", con, sep="\n")
	for (i in 1:model$n.edges)
	{
		x1 <- model$edges[i,1]
		x2 <- model$edges[i,2]
		y1 <- result[x1]
		y2 <- result[x2]
		if (y1 > label$size || y2 > label$size)
		{
			writeLines(paste(query.name[x1], "--", query.name[x2], "\t", label.name[y1], "--", label.name[y2], "\tgap"), con, sep="\n")
		}
		else
		{
			writeLines(paste(query.name[x1], "--", query.name[x2], "\t", label.name[y1], "--", label.name[y2], "\t", label$dist[y1, y2]), con, sep="\n")
		}
	}
	close(con)
}
