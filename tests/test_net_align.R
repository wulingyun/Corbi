library(Corbi)
library(tools)

setwd("net_align")

net_align("net1.txt", "net3.txt", "sim.txt")
if (Rdiff("answer1.txt", "net1.txt_result.txt") != 0) 
  stop("Failed in net_align.")
if (Rdiff("answer1_uni.txt", "list net1.txt_result.txt") != 0) 
  stop("Failed in net_align.")
