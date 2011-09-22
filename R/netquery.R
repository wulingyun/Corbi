net.query <- function(query.net, target.net, node.sim, query.type=2, delta.d=1e-10, delta.c=0.5, delta.e=1, delta.s=1, normalize=F, multitest=1)
{
	# options
	#   query.net: input file name of query network
	#   target.net: input file name of target network
	#   node.sim: input file name of node similarity
	#   query.type: the type of query network, 0 - chain, 1 - tree, 2 - loopy
	#   delta: the parameters \Delta_d, \Delta_c, \Delta_e, \Delta_s
	#   isnormalize: whether carry out the normalized operation to the node feature values or not? default not
	#   multitest: the result files' number; only using in many experiments

	# read the input files
	data <- read.data(query.net, target.net, node.sim)

	# data:
	#   query$node: the node names of query network
	#   query$matrix: the adjacency matrix of query network
	#   target$node: the node names of target network
	#   target$matrix: the adjacency matrix of target network
	#   node.sim: the node similarity matrix

	# normalize the node similarity matrix
	if (normalize || max(data$node.sim) <= 0)
	{
		a <- min(data$node.sim)
		b <- max(data$node.sim)
		data$node.sim <- (data$node.sim - a) / (b - a)
	}

	# compute the shortest distance matrix for the target network
	data$target$dist <- .Call("NA_ShortestDistances", data$target$matrix)

	# simplify the target network
	data <- simplify.net(data, delta.d);

	emp <- GEM(delta.d, data$node.sim);

	# compute the transition probability
	D_X = data$query$matrix;
	s<-D_X-t(D_X);
	s1<-sum(s[s>0]);
	s2<-sum(s[s<0]);
	if(s1==0 & s2==0){F1<-tranf1(data$target$dist,emp,D_X,delta.s,delta.e,delta.c);}
	else{F1<-tranf2(data$target$dist,emp,D_X,delta.s,delta.e,delta.c);}

	D_X<-D_X+t(D_X);
	D_X[D_X>1]<-1;

	Y<-dyf(data$query$node,emp,F1$F1,data$target$node,data$target$dist,D_X,query.type);

	resultxt(data$query$node,Y,D_X,data$target$node,data$target$dist,multitest)
}

read.net <- function(net)
{
	n.col <- max(sapply(strsplit(as.matrix(read.table(net, as.is=T, sep="\n")), " "), length))
	net.text <- as.matrix(read.table(net, fill=T, as.is=T, col.names=1:n.col))
	net.node <- unique(as.vector(net.text))
	net.node <- net.node[net.node != ""]
	net.size <- length(net.node)
	net.edge <- cbind(net.text[,1], as.vector(net.text[,-1]))
	net.edge <- net.edge[net.edge[,2] != "", ]
	net.matrix <- matrix(0, net.size, net.size)
	rownames(net.matrix) <- colnames(net.matrix) <- net.node
	net.matrix[net.edge] <- 1
	list(size=net.size, node=net.node, matrix=net.matrix)
}

read.node.sim <- function(query, target, node.sim)
{
	sim.text <- read.table(node.sim, as.is=T)
	sim.text <- sim.text[sim.text[,1] %in% query$node, ]
	sim.text <- sim.text[sim.text[,2] %in% target$node, ]
	sim.matrix <- matrix(0, nrow=query$size, ncol=target$size)
	rownames(sim.matrix) <- query$node
	colnames(sim.matrix) <- target$node
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

simplify.net <- function(data, delta.d)
{
	s <- colSums(data$node.sim > delta.d) > 0
	data$target$size <- sum(s)
	if (data$target$size > 0)
	{
		data$target$node <- data$target$node[s]
		data$target$matrix <- data$target$matrix[s, s]
		data$node.sim <- data$node.sim[,s]
		data$node.sim[data$node.sim <= delta.d] <- 0
		d <- data$target$dist[s, s]
		d[d == -1] <- Inf
		for (i in 1:dim(d)[1]) d[i, i] <- 1
		data$target$dist <- d
	}
	else
	{
		stop("The simplified target network is empty, you should give a smaller cut values!")
	}
	data
}
