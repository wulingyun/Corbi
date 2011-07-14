afp.score <- function(D1, D2, i1, i2, afp.length = 5)
{
	c1 <- i1:(i1+afp.length-1)
	c2 <- i2:(i2+afp.length-1)
	d <- D1[c1, c1] - D2[c2, c2]
	sqrt(mean((d[upper.tri(d)])^2))
}

afp.match <- function(D1, D2, afp.length = 5)
{
	n1 <- dim(D1)[1] - afp.length + 1
	n2 <- dim(D2)[1] - afp.length + 1
	m.afp <- matrix(0, nrow=n1, ncol=n2)
	for (i in 1:n1)
		m.afp[i, 1:n2] <- sapply(1:n2, function(x) afp.score(D1, D2, i, x, afp.length))
	m.afp
}

afp.dist <- function(D1, D2, n.nodes, n.afp, afp.length = 5)
{
	.Call("AFP_Distance", D1, D2, n.nodes, n.afp, afp.length)
}

pro.align <- function(D1, D2)
{
	afp.cutoff <- 1
	afp.length <- 5
	gap.penalty <- 0.5

	m.afp <- afp.match(D1, D2, afp.length)
	m.afp[m.afp >= afp.cutoff] <- 1
	m.afp <- 1 - m.afp

	n.d1 <- dim(D1)[1]
	n.d2 <- dim(D2)[1]

	n.nodes <- n.d1
	n.afp <- n.d2 - afp.length + 1
	n.states <- n.afp * afp.length + 2

	adj <- matrix(0, nrow=n.nodes, ncol=n.nodes)
	adj[cbind(1:(n.nodes-1), 2:n.nodes)] <- 1

	crf <- make.crf(adj, n.states)

	crf$node.pot[,] <- 0
	crf$node.pot[,1:2] <- gap.penalty
	for (i in 1:afp.length)
		crf$node.pot[i:(i+n.nodes-afp.length),(n.afp*(i-1)+3):(n.afp*i+2)] <- m.afp

	ep <- matrix(0, nrow=n.states, ncol=n.states)
	ep[1,1:(n.afp+2)] <- 1
	ep[,2] <- 1
	ep[cbind(3:(n.afp*(afp.length-1)+2), (n.afp+3):(n.afp*afp.length+2))] <- 1

	m.dist <- afp.dist(D1, D2, n.nodes, n.afp, afp.length)
	m.dist[m.dist >= afp.cutoff] <- 1
	m.dist <- 1 - m.dist

	for (e in 1:crf$n.edges)
	{
		crf$edge.pot[,,e] <- ep
		if (e <= n.nodes-afp.length)
			for (i in 1:min(afp.length, e))
				crf$edge.pot[(n.afp*(i-1)+3):(n.afp*i+2), 3:(n.afp+2), e] <- m.dist[, , i, e]
	}

	dec <- decode.chain(crf)
	dec[dec == 1 | dec == 2] <- 0
	dec[dec != 0] <- dec[dec != 0] - 2
	dec <- dec %% n.afp + dec %/% n.afp
	dec
}
