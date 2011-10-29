
PA.scores <- function(D1, D2, afp.length = 5)
{
	.Call("PA_Scores", D1, D2, dim(D1)[1], dim(D2)[1], afp.length)
}

tri2mat <- function(v, n)
{
	m <- array(0, dim=c(n,n))
	m[upper.tri(m)] <- v
	m[lower.tri(m)] <- t(m)[lower.tri(m)]
	m
}

mat2tri <- function(m)
{
	m[upper.tri(m)]
}

pro.align <- function(D1, D2, gap.penalty1 = 0.5, gap.penalty2 = 0.5, afp.length = 5, cutoff = 0.1)
{
	n.d1 <- dim(D1)[1]
	n.d2 <- dim(D2)[1]
	n.nodes <- n.d1
	n.states <- n.d2 * 2

	adj <- matrix(0, nrow=n.nodes, ncol=n.nodes)
	for (i in 2:n.nodes) adj[i-1, i] <- 1
	crf <- make.crf(adj, n.states)

	scores <- PA.scores(D1, D2, afp.length)
	scores$node.score <- exp(-scores$node.score)
	scores$node.score[scores$node.score <= exp(-cutoff)] <- 0
	scores$edge.score <- exp(-scores$edge.score)
	scores$edge.score[scores$edge.score <= exp(-cutoff)] <- 0

	crf$node.pot[,] <- 1
	crf$node.pot[,1:n.d2] <- scores$node.score

	ep <- matrix(0, nrow=n.states, ncol=n.states)
	ep[cbind(1:n.d2, (n.d2+1):n.states)] <- gap.penalty1
	ep[cbind((n.d2+1):n.states, (n.d2+1):n.states)] <- gap.penalty2
	temp <- matrix(gap.penalty1, nrow=n.d2, ncol=n.d2)
	temp[lower.tri(temp, diag=T)] <- 0
	ep[(n.d2+1):n.states, 1:n.d2] <- temp

	for (e in 1:crf$n.edges)
	{
		crf$edge.pot[[e]] <- ep
		crf$edge.pot[[e]][1:n.d2, 1:n.d2] <- scores$edge.score[,, crf$edges[e,1]]
	}

	label <- decode.chain(crf)
	label[label > n.d2] <- 0
	label
}
