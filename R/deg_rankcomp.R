#' @export
rankcomp <- function (expr.ctrl, expr.case, cutoff = 0.99) {
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
    
    up_score[k,] <- -log(p_value)
    up_score[k, n_up_rev > n_down_rev] <- 0
    down_score[k,] <- -log(p_value)
    down_score[k, n_up_rev < n_down_rev] <- 0
  }
  return(list(up = up_score, down = down_score))
}
