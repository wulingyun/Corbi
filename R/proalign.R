
afp.score <- function(D1, D2, afp.length = 5)
{
	.Call("AFP_Score", D1, D2, dim(D1)[1], dim(D2)[1], afp.length)
}

afp.dist <- function(D1, D2, n.nodes, n.afp, afp.length = 5)
{
	.Call("AFP_Distance", D1, D2, n.nodes, n.afp, afp.length)
}

pro.align <- function(D1, D2, gap.penalty = 0.1, allow.mid.gap = T, afp.length = 5, afp.cutoff = 1)
{
	n.d1 <- dim(D1)[1]
	n.d2 <- dim(D2)[1]
	n.nodes <- n.d1
	n.afp <- n.d2 - afp.length + 1
	n.states <- n.afp * afp.length + 2

	adj <- matrix(0, nrow=n.nodes, ncol=n.nodes)
	adj[cbind(1:(n.nodes-1), 2:n.nodes)] <- 1
	crf <- make.crf(adj, n.states)

	m.afp <- afp.score(D1, D2, afp.length)
	m.afp[m.afp >= afp.cutoff] <- 1
	m.afp <- 1 - m.afp

	crf$node.pot[,] <- 0
	crf$node.pot[,1:2] <- gap.penalty
	for (i in 1:afp.length)
		crf$node.pot[i:(i+n.nodes-afp.length),(n.afp*(i-1)+3):(n.afp*i+2)] <- m.afp

	ep <- matrix(0, nrow=n.states, ncol=n.states)
	ep[1,1:(n.afp+2)] <- 1
	ep[c(2,(n.afp*(afp.length-1)+3):(n.afp*afp.length+2)),2] <- 1
	if (allow.mid.gap) ep[2,2:(n.afp+2)] <- 1
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

	core <- decode.chain(crf)
	core[core == 1 | core == 2] <- 0
	core[core != 0] <- core[core != 0] - 2
	core <- core %% n.afp + core %/% n.afp
	core
}
