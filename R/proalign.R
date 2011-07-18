
afp.score <- function(D1, D2, afp.length = 5)
{
	.Call("AFP_Score", D1, D2, dim(D1)[1], dim(D2)[1], afp.length)
}

afp.dist <- function(D1, D2, n.nodes, n.afp, afp.length = 5)
{
	.Call("AFP_Distance", D1, D2, n.nodes, n.afp, afp.length)
}

get.seed <- function(core)
{
	seed <- matrix(0, nrow=0, ncol=2)
	in.seed <- F
	for (i in 1:length(core))
	{
		if (!in.seed && core[i] != 0)
		{
			in.seed <- T
			start <- i
		}
		if (in.seed && core[i] == 0)
		{
			in.seed <- F
			end <- i-1
			seed <- rbind(seed, c(start, end))
		}
	}
	if (in.seed)
	{
		end <- length(core)
		seed <- rbind(seed, c(start, end))
	}
	seed
}

tri2mat <- function(v, n)
{
	m <- array(0, dim=c(n,n))
	m[upper.tri(m, diag=T)] <- v
	m[lower.tri(m)] <- t(m)[lower.tri(m)]
	m
}

mat2tri <- function(m)
{
	m[upper.tri(m, diag=T)]
}

get.alignment <- function(D1, D2, seed, core, gap.penalty = 0.5)
{
	n.d1 <- dim(D1)[1]
	n.d2 <- dim(D2)[1]
	n.seed <- dim(seed)[1]
	if (n.seed <= 0) return()

	match.score <- sapply(mat2tri(D1), function(x) exp(-(x - mat2tri(D2))^2))
	map.d1 <- tri2mat(1:dim(match.score)[2], n.d1)
	map.d2 <- tri2mat(1:dim(match.score)[1], n.d2)

	n.states <- n.d2 * 2
	gap.states <- 1:n.d2
	normal.states <- (n.d2+1):n.states
	ep <- matrix(0, nrow=n.states, ncol=n.states)
	diag(ep[gap.states, gap.states]) <- 1
	diag(ep[normal.states, gap.states]) <- 1
	temp <- matrix(0, nrow=n.d2, ncol=n.d2)
	temp[upper.tri(temp)] <- 1
	ep[gap.states, normal.states] <- temp

	alignment <- matrix(0, nrow=n.seed, ncol=n.d1)
	for (i in 1:n.seed)
	{
		fixed <- rep(F, n.d1)
		fixed[seed[i,1]:seed[i,2]] <- T
		n.fixed <- sum(fixed)
		n.nodes <- n.d1 - n.fixed
		if (n.fixed >= n.d1)
		{
			alignment[i,] <- core
			break
		}
		fixed.id <- (1:n.d1)[fixed]
		node.id <- (1:n.d1)[!fixed]

		adj <- matrix(0, nrow=n.nodes, ncol=n.nodes)
		for (j in 2:n.nodes) adj[j-1, j] <- 1

		crf <- make.crf(adj, n.states)

		crf$node.pot[,] <- 0
		crf$node.pot[,gap.states] <- gap.penalty

		u1 <- core[seed[i,1]] - 1
		u2 <- core[seed[i,2]] + 1
		for (j in 1:n.nodes)
		{
			if (node.id[j] < seed[i,1] && u1 >= 1)
			{
				crf$node.pot[j,n.d2+(1:u1)] <- 1
				for (k in 1:n.fixed)
				{
					temp <- match.score[map.d2[1:u1, core[fixed.id[k]]], map.d1[node.id[j], fixed.id[k]]]
					crf$node.pot[j,n.d2+(1:u1)] <- crf$node.pot[j,n.d2+(1:u1)] * temp
				}
			}
			if (node.id[j] > seed[i,2] && u2 <= n.d2)
			{
				crf$node.pot[j,n.d2+(u2:n.d2)] <- 1
				for (k in 1:n.fixed)
				{
					temp <- match.score[map.d2[core[fixed.id[k]], u2:n.d2], map.d1[fixed.id[k], node.id[j]]]
					crf$node.pot[j,n.d2+(u2:n.d2)] <- crf$node.pot[j,n.d2+(u2:n.d2)] * temp
				}
			}
		}

		crf$edge.pot[,,] <- ep
		for (j in 1:crf$n.edges)
		{
			temp <- matrix(match.score[map.d2[1:n.d2, 1:n.d2], map.d1[node.id[crf$edges[j,1]], node.id[crf$edges[j,2]]]], nrow=n.d2, ncol=n.d2)
			temp[lower.tri(temp, diag=T)] <- 0
			crf$edge.pot[normal.states, normal.states, j] <- temp
		}

		dec <- decode.chain(crf) - n.d2
		dec[dec < 0] <- 0
		alignment[i, node.id] <- dec
		alignment[i, fixed.id] <- core[fixed]
	}
	alignment
}

pro.align <- function(D1, D2, gap.penalty = 0.5, allow.mid.gap = T, afp.length = 5, afp.cutoff = 1)
{
	n.d1 <- dim(D1)[1]
	n.d2 <- dim(D2)[1]
	n.nodes <- n.d1
	n.afp <- n.d2 - afp.length + 1
	n.states <- n.afp * afp.length + 2

	adj <- matrix(0, nrow=n.nodes, ncol=n.nodes)
	adj[cbind(1:(n.nodes-1), 2:n.nodes)] <- 1
	crf <- make.crf(adj, n.states)

	m.afp <- exp(-afp.score(D1, D2, afp.length))
	m.afp[m.afp <= exp(-afp.cutoff)] <- 0

	crf$node.pot[,] <- 0
	crf$node.pot[,1:2] <- gap.penalty
	for (i in 1:afp.length)
		crf$node.pot[i:(i+n.nodes-afp.length),(n.afp*(i-1)+3):(n.afp*i+2)] <- m.afp

	ep <- matrix(0, nrow=n.states, ncol=n.states)
	ep[1,1:(n.afp+2)] <- 1
	ep[c(2,(n.afp*(afp.length-1)+3):(n.afp*afp.length+2)),2] <- 1
	if (allow.mid.gap) ep[2,2:(n.afp+2)] <- 1
	ep[cbind(3:(n.afp*(afp.length-1)+2), (n.afp+3):(n.afp*afp.length+2))] <- 1

	m.dist <- exp(-afp.dist(D1, D2, n.nodes, n.afp, afp.length))
	m.dist[m.dist <= exp(-afp.cutoff)] <- 0

	for (e in 1:crf$n.edges)
	{
		crf$edge.pot[,,e] <- ep
		if (e <= n.nodes-afp.length)
			for (i in 1:min(afp.length, e))
				crf$edge.pot[(n.afp*(i-1)+3):(n.afp*i+2), 3:(n.afp+2), e] <- m.dist[, , i, e]
	}

	core <- decode.chain(crf)
	core.gap <- core == 1 | core == 2
	core <- core - 3
	core <- core %% n.afp + core %/% n.afp + 1
	core[core.gap] <- 0

	result <- list()
	result$core <- core
	result$seed <- get.seed(core)
	result$alignment <- get.alignment(D1, D2, result$seed, core, gap.penalty)
	result
}
