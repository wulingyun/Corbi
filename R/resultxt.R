resultxt<-function(X,Y,D_X,name,g,multitest){

# add the gap node
N<-dim(g)[1];
name[N+1]<-"gap";
g1<-matrix(1,1,N);g2<-matrix(1,N+1,1);
g<-rbind(g,g1);
g<-cbind(g,g2);

filename<-paste("result",as.character(multitest),sep="");
con<- file(filename,"w");
writeLines("node match:",con,sep = "\n");
for (i in 1:length(X)){
writeLines(paste(X[i],"  ",name[Y[i]]),con,sep = "\n");
}
writeLines("",con,sep = "\n");
writeLines("edge match:",con,sep = "\n")
K<-sum(D_X[upper.tri(D_X)]);

for (i in 1:(dim(D_X)[1]-1)){
     for (j in (i+1):dim(D_X)[1]){
            if(D_X[i,j]==1){
                 if(g[Y[i],Y[j]]>1){
                 writeLines(paste(X[i],"--",X[j],"\t",name[Y[i]],"--",name[Y[j]],"   gap"),con,sep = "\n");}
                 else{
                 writeLines(paste(X[i],"--",X[j],"\t",name[Y[i]],"--",name[Y[j]]),con,sep = "\n");}
                           }
     }
}# end for i


close(con)
}