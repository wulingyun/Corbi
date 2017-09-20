#' @export
getRatioDistribution <- function(expr) {
  n <- dim(expr)[1]
  ratio_mean <- matrix(0, n, n)
  ratio_sd <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      r <- (expr[i,] - expr[j,]) / (expr[i,] + expr[j,])
      ratio_mean[i,j] <- mean(r, na.rm = T)
      ratio_sd[i,j] <- sd(r, na.rm = T)
    }
  }
  ratio_mean <- ratio_mean - t(ratio_mean)
  ratio_sd <- ratio_sd - t(ratio_sd)
  return(list(mean = ratio_mean, sd = ratio_sd))
}
      

#' @export
getRatioNet <- function(ratio.dist, expr, cutoff = 1.0)
{
  n <- length(expr)
  z_matrix <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      r <- (expr[i] - expr[j]) / (expr[i] + expr[j])
      z <- (r - ratio.dist$mean[i,j]) / ratio.dist$sd[i,j]
      if (!is.na(z)) {
        if (z > 0) z_matrix[i,j] <- z
        else z_matrix[j,i] <- -z
      }
    }
  }
  adj_matrix <- ifelse(z_matrix > cutoff, 1, 0)
  return(list(adj = adj_matrix, z = z_matrix))
}



inferByDegree <- function(adj.matrix) {
  in_degree <- colSums(adj.matrix)
  out_degree <- rowSums(adj.matrix)
  up_score <- out_degree - in_degree
  down_score <- -up_score
  up_score[up_score < 0] <- 0
  down_score[down_score < 0] <- 0
  return(list(up = up_score, down = down_score, all = up_score-down_score))
}


inferByVertexCover_lpSolve <- function(adj.matrix, z.matrix, coverage = 1.0, binary = FALSE) {
  require(lpSolve)
  
  # LP variables:
  #  X indicate node is up
  #  Y indicate node is down  
  #  Z indicate edge is covered
  
  n_nodes <- dim(adj.matrix)[1]
  n_edges <- sum(adj.matrix)
  n_variables <- n_nodes * 2 + n_edges
  
  if (n_edges <= 0) {
    scores <- list(up = numeric(n_nodes), down = numeric(n_nodes))
  }
  else {
    edges <- which(adj.matrix == 1, arr.ind = TRUE)
    
    # objective function
    obj_fun <- c(rep(1, n_nodes * 2), rep(0, n_edges))
    
    # 1. coefficient matrix of X, Y, Z for Xi + Yi <= 1, for any node i in V(G)
    coeff1 <- matrix(0, n_nodes, n_variables)
    coeff1[cbind(1:n_nodes, 1:n_nodes)] <- 1
    coeff1[cbind(1:n_nodes, (n_nodes+1):(n_nodes*2))] <- 1
    
    # 2. coefficient matrix of X, Y, Z for Xi + Yj >= Zij, for any edge (i, j) in E(G)
    coeff2 <- matrix(0, n_edges, n_variables)
    for (e in 1:n_edges) {
      coeff2[e, edges[e,1]] <- 1
      coeff2[e, n_nodes+edges[e,2]] <- 1
      coeff2[e, n_nodes*2+e] <- -1
    }
    
    # 3. coefficient vector of X, Y, Z for sum(Z) >= k                         #  k is edge.num.covered
    coeff3 <- matrix(0, 1, n_variables)
    coeff3[1, (n_nodes*2+1):n_variables] <- 1
    
    constraint_mat <- rbind(coeff1, coeff2, coeff3)
    constraint_rhs <- c(rep(1, n_nodes), rep(0, n_edges), ceiling(n_edges*coverage))
    constraint_dir <- c(rep("<=", n_nodes), rep(">=", n_edges), ">=")
    
    # variable in interval [0, 1]
    variable_mat <- diag(x = 1, n_variables, n_variables)
    variable_min <- rep(0, n_variables)
    variable_max <- rep(1, n_variables)
    
    # final constraints
    constr_mat <- rbind(constraint_mat, variable_mat, variable_mat)
    constr_dir <- c(constraint_dir, rep( ">=", n_variables), rep( "<=", n_variables))
    constr_rhs <- c(constraint_rhs, variable_min, variable_max)
    
    # LP solution
    optimum <- lp(direction = "min",
                  objective.in = obj_fun,
                  const.mat = constr_mat,
                  const.dir = constr_dir,
                  const.rhs = constr_rhs,
                  transpose.constraints = TRUE,
                  all.bin = binary)
    
    x <- optimum$solution[1:n_nodes]
    y <- optimum$solution[(n_nodes+1):(n_nodes*2)]
    
    scores <- calcVertexCoverScore(x, y, edges, z.matrix, binary)
  }
  return(scores)
}


inferByVertexCover <- function(adj.matrix, z.matrix, coverage = 1.0, binary = FALSE) {
  require(Matrix)
  require(Rglpk)
  
  # LP variables:
  #  X indicate node is up
  #  Y indicate node is down  
  #  Z indicate edge is covered
  
  n_nodes <- dim(adj.matrix)[1]
  n_edges <- sum(adj.matrix)
  n_variables <- n_nodes * 2 + n_edges
  
  if (n_edges <= 0) {
    scores <- list(up = numeric(n_nodes), down = numeric(n_nodes))
  }
  else {
    edges <- which(adj.matrix == 1, arr.ind = TRUE)
    
    # objective function
    obj_fun <- c(rep(1, n_nodes * 2), rep(0, n_edges))
    
    # 1. coefficient matrix of X, Y, Z for Xi + Yi <= 1, for any node i in V(G)
    vi <- rep(seq(n_nodes), 2)
    vj <- seq(n_nodes*2)
    vx <- rep(1, n_nodes*2)

    # 2. coefficient matrix of X, Y, Z for Xi + Yj >= Zij, for any edge (i, j) in E(G)
    vi <- c(vi, rep(seq(n_edges), 3)+n_nodes)
    vj <- c(vj, edges[,1], edges[,2]+n_nodes, seq(n_edges)+n_nodes*2)
    vx <- c(vx, rep(1, n_edges*2), rep(-1, n_edges))

    # 3. coefficient vector of X, Y, Z for sum(Z) >= k
    vi <- c(vi, rep(n_nodes+n_edges+1, n_edges))
    vj <- c(vj, seq(n_edges)+n_nodes*2)
    vx <- c(vx, rep(1, n_edges))

    constraint_mat <- sparseMatrix(vi, vj, x = vx, dims = c(n_nodes+n_edges+1, n_variables))
    constraint_rhs <- c(rep(1, n_nodes), rep(0, n_edges), ceiling(n_edges*coverage))
    constraint_dir <- c(rep("<=", n_nodes), rep(">=", n_edges), ">=")
    
    # variable in interval [0, 1]
    variable_bounds <- list(lower = list(ind = 1:n_variables, val = rep(0, n_variables)),
                            upper = list(ind = 1:n_variables, val = rep(1, n_variables)))
    
    # LP solution
    optimum <- Rglpk_solve_LP(obj = obj_fun,
                              mat = constraint_mat,
                              dir = constraint_dir,
                              rhs = constraint_rhs,
                              bounds = variable_bounds,
                              types = ifelse(binary, "B", "C"),
                              max = FALSE)
    
    x <- optimum$solution[1:n_nodes]
    y <- optimum$solution[(n_nodes+1):(n_nodes*2)]
    
    scores <- calcVertexCoverScore(x, y, edges, z.matrix, binary)
  }
  return(scores)
}


calcVertexCoverScore <- function(x, y, edges, z.matrix, binary) {
  n_nodes <- length(x)
  if (binary) {
    up_score <- x
    down_score <- y
  }
  else {
    z <- x[edges[,1]] + y[edges[,2]]
    z1 <- x[edges[,1]] / z
    z2 <- y[edges[,2]] / z
    z1[is.na(z1) | is.infinite(z1)] <- 0
    z2[is.na(z2) | is.infinite(z2)] <- 0
    up_score <- sapply(1:n_nodes, function(i) sum(z1[edges[,1] == i]))
    down_score <- sapply(1:n_nodes, function(i) sum(z2[edges[,2] == i]))
    z1w <- (z1 + 0.5) * z.matrix[edges]
    z2w <- (z2 + 0.5) * z.matrix[edges]
    up_score_w <- sapply(1:n_nodes, function(i) sum(z1w[edges[,1] == i]))
    down_score_w <- sapply(1:n_nodes, function(i) sum(z2w[edges[,2] == i]))
  }
  return(list(up = up_score, down = down_score, all = up_score-down_score, allw = up_score_w-down_score_w))
}


#' @export
inferByRankComp <- function (expr.ctrl, expr.case, cutoff = 0.99) {
  n_genes <- dim(expr.ctrl)[1]
  n_ctrls <- dim(expr.ctrl)[2]
  n_cases <- dim(expr.case)[2]
  
  up_score <- matrix(0, n_genes, n_cases)
  down_score <- matrix(0, n_genes, n_cases)
  for (k in 1:n_genes) {
    diff <- expr.ctrl[k,] - t(expr.ctrl)
    sp_up <- which((colSums(diff > 0)/n_ctrls) > cutoff)
    sp_down <- which((colSums(diff < 0)/n_ctrls) > cutoff)
    n_up <- length(sp_up)
    n_down <- length(sp_down)
    
    n_up_rev <- numeric(n_cases)
    n_down_rev <- numeric(n_cases)
    if (n_up > 0) {
      diff <- expr.case[k,] - t(expr.case[sp_up, , drop = F])
      n_up_rev <- rowSums(diff < 0)
    }
    if (n_down > 0) {
      diff <- expr.case[k,] - t(expr.case[sp_down, , drop = F])
      n_down_rev <- rowSums(diff > 0)
    }
    
    sp_counts <- cbind(n_up, n_up-n_up_rev+n_down_rev, n_down, n_down-n_down_rev+n_up_rev)
    p_value <- apply(sp_counts, 1, function(x) fisher.test(matrix(x, 2, 2))$p.value)
    
    up_score[k,] <- 1 - p_value
    up_score[k, n_up_rev > n_down_rev] <- 0
    down_score[k,] <- 1 - p_value
    down_score[k, n_up_rev < n_down_rev] <- 0
  }
  return(list(up = up_score, down = down_score))
}
