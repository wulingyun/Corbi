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

	data <- read.data(query.net, target.net, node.sim)

# data:
#   query$node: the node names of query network
#   query$matrix: the adjacency matrix of query network
#   target$node: the node names of target network
#   target$matrix: the adjacency matrix of target network
#   node.sim: the node similarity matrix

	data$delta <- list(d=delta.d, c=delta.c, e=delta.e, s=delta.s)

# compute the shortest distance matrix for the target network
# and simplify the target network

	data$labels <- simplify.net(data)

# build and solve CRF model

	data$model <- build.model(data)
	data$result <- solve.crf(data, query.type)

# write result to output file

	write.result(data, output)
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

read.node.sim <- function(query, target, node.sim)
{
	sim.text <- read.table(node.sim, as.is=T)
	sim.text <- sim.text[sim.text[,1] %in% query$node, ]
	sim.text <- sim.text[sim.text[,2] %in% target$node, ]
	sim.matrix <- matrix(0, query$size, target$size, dimnames=list(query$node, target$node))
	sim.matrix[as.matrix(sim.text[,1:2])] <- sim.text[,3]
	sim.matrix
}

read.data <- function(query.net, target.net, node.sim)
{
	query <- read.net(query.net)
	target <- read.net(target.net)
	similarity <- read.node.sim(query, target, node.sim)
	list(query=query, target=target, node.sim=similarity)
}

simplify.net <- function(data)
{
	delta.d <- data$delta$d
	select <- colSums(data$node.sim > delta.d) > 0
	net.size <- sum(select)
	if (net.size <= 0)
	{
		stop("The simplified target network is empty, you should give a smaller cut values!")
	}
	net.node <- data$target$node[select]
	net.matrix <- data$target$matrix[select, select]
	node.sim <- data$node.sim[, select]
	node.sim[node.sim <= delta.d] <- 0
	net.dist <- .Call("NQ_ShortestDistances", data$target$matrix, select)[select, select]
	net.dist[net.dist == -1] <- Inf
	net.dist[cbind(1:net.size, 1:net.size)] <- 1
	list(size=net.size, node=net.node, matrix=net.matrix, dist=net.dist, node.sim=node.sim)
}

build.model <- function(data)
{
	S <- cbind(data$labels$node.sim, data$delta$d)
	node.pot <- S
	
	query.size <- data$query$size
	query.net = data$query$matrix
	query.net[query.net != 0] <- 1
	label.size <- data$labels$size
	label.gap <- label.size + 1
	label.weight <- 1 / data$labels$dist ^ data$delta$s

	W <- matrix(nrow=label.gap, ncol=label.gap)
	W[1:label.size, 1:label.size] <- label.weight
	W[label.gap,] <- data$delta$e
	W[,label.gap] <- data$delta$e
	W[label.gap, label.gap] <- data$delta$c

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

solve.crf <- function(data, query.type)
{
	if (query.type == 0)
	{
		result <- decode.chain(data$model)
	}
	if (query.type == 1)
	{
		result <- decode.tree(data$model)
	}
	if (query.type == 2)
	{
		result <- decode.lbp(data$model)
	}
	result
}

write.result <- function(data, filename="result.txt")
{
	query.name <- data$query$node
	target.name <- c(data$labels$node, "gap")
	label.size <- data$labels$size
	label.dist <- data$labels$dist

	con <- file(as.character(filename), "w")
	writeLines("node match:", con, sep="\n")
	for (i in 1:data$query$size)
	{
		writeLines(paste(query.name[i], "  ", target.name[data$result[i]]), con, sep="\n");
	}
	writeLines("", con, sep="\n")

	writeLines("edge match:", con, sep="\n")
	for (i in 1:data$model$n.edges)
	{
		x1 <- data$model$edges[i,1]
		x2 <- data$model$edges[i,2]
		y1 <- data$result[x1]
		y2 <- data$result[x2]
		if (y1 > label.size || y2 > label.size)
		{
			writeLines(paste(query.name[x1], "--", query.name[x2], "\t", target.name[y1], "--", target.name[y2], "\tgap"), con, sep="\n")
		}
		else
		{
			writeLines(paste(query.name[x1], "--", query.name[x2], "\t", target.name[y1], "--", target.name[y2], "\t", label.dist[y1, y2]), con, sep="\n")
		}
	}

	close(con)
}
