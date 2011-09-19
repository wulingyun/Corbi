dyf<-function(X,emp,F1,name,g,D_X,t)
{
 

C<-dim(emp);


nnodes<-C[1];
nstates<-C[2];

crf<-make.crf(D_X,nstates);

for (i in 1:nnodes) {crf$node.pot[i,]<-emp[i,]}

for (i in 1:crf$n.edges) {crf$edge.pot[,,i]<-F1[,,i]}


if(t==0){Y<-decode.chain(crf)}
if(t==1){Y<-decode.tree(crf)}
if(t==2){Y<-decode.lbp(crf)}

return(Y)
#list(Y,Y1,D)
}