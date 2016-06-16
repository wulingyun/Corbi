library(Corbi)
library(tools)

setwd("net_query")

net_query("query1.txt", "target.txt", "sim.txt")
if (Rdiff("answer1.txt", "query1.txt_result.txt") != 0) 
  stop("Failed in net_query.")

query.nets <- c("query1.txt", "query2.txt")
net_query_batch(query.nets, "target.txt", "sim.txt")
if (Rdiff("answer1.txt", "query1.txt_result.txt") != 0)
  stop("Failed in net_query_batch.")
if (Rdiff("answer2.txt", "query2.txt_result.txt") != 0)
  stop("Failed in net_query_batch.")
