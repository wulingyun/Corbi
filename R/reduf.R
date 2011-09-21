reduf<-function(E,N,G,deltad)

{
E[E<=deltad]<-0;
s<-colSums(E);
b<-s>0;

g=G[b,b];

if(length(g)>0){
g[g==0]=Inf; ## this may be change as the distance matrix method

name<-c("");
name=N[b];
for (j in 1:dim(g)[1]) {if(g[j,j]==Inf){g[j,j]=2;}}

E1<-E[,b];}
else{print("The simplification network is empty, you should give a smaller cut values!")}

list(g=g,name=name,E1=E1)

}