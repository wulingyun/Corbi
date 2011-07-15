
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

get.match.score <- function(d1, d2)
{
	(1 - abs((d1 - d2) / (d1 + d2)))
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
	core.gap <- core == 1 | core == 2
	core <- core - 3
	core <- core %% n.afp + core %/% n.afp + 1
	core[core.gap] <- 0

	seed <- get.seed(core)
	n.seed <- dim(seed)[1]

	match.score <- sapply(mat2tri(D1), function(x) get.match.score(x, mat2tri(D2)))
	map.d1 <- tri2mat(1:dim(match.score)[2], n.d1)
	map.d2 <- tri2mat(1:dim(match.score)[1], n.d2)

	alignment <- matrix(0, nrow=n.seed, ncol=n.d1)
	for (i in 1:n.seed)
	{
		fixed <- rep(F, n.nodes)
		fixed[seed[i,1]:seed[i,2]] <- T
		if (sum(fixed) >= n.d1)
		{
			alignment[i,] <- core
			break
		}

		normal.nodes <- (1:n.nodes)[!fixed]
		adj <- matrix(0, nrow=n.nodes, ncol=n.nodes)
		for (j in 2:length(normal.nodes)) adj[normal.nodes[j-1], normal.nodes[j]] <- 1
		adj[fixed,] <- 1
		diag(adj) <- 0

		n.states <- n.d2 + 1

		crf <- make.crf(adj, n.states)

		crf$node.pot[,] <- 1
		crf$node.pot[,1] <- gap.penalty

		crf$edge.pot[1, , ] <- 1
		crf$edge.pot[, 1, ] <- 1

		for (j in 1:crf$n.edges)
		{
			m <- matrix(match.score[map.d2[1:n.d2, 1:n.d2], map.d1[crf$edges[j,1], crf$edges[j,2]]], nrow=n.d2, ncol=n.d2)
			m[lower.tri(m, diag=T)] <- 0
			crf$edge.pot[2:(n.d2+1), 2:(n.d2+1), j] <- m
		}

		clamp = rep(0, n.nodes)
		clamp[fixed] <- core[fixed] + 1
		print(clamp)
		alignment[i,] <- decode.conditional(crf, clamp, decode.chain) - 1
	}

	result <- list()
	result$alignment <- alignment
	result$core <- core
	result$seed <- seed
	result
}
