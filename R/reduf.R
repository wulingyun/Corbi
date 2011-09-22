reduf<-function(E,N,G,deltad)

{
E[E<=deltad]<-0;
s<-colSums(E);
b<-s>0;

g=G[b,b];

if(length(g)>0){
g[g==-1]=Inf; ## this may be change as the distance matrix method

name<-c("");
name=N[b];
g[cbind(1:dim(g)[1],1:dim(g)[1])]=1;

E1<-E[,b];}
else{print("The simplification network is empty, you should give a smaller cut values!")}

list(g=g,name=name,E1=E1)

}