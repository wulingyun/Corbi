tranf1<-function(D,emp,D_X,s,e,c)
{

D<-1/D^s;C<-dim(D);
C1<-dim(D_X);C2<-dim(emp);
E=emp[,1:(C2[2]-1)];


k<-1;
K<-sum(D_X[upper.tri(D_X)]);
F1<-array(0,dim=c(C[1]+1,C[1]+1,K))

for (i in 1:(C1[1]-1))
    { for (j in (i+1):C1[1])
         { if(D_X[i,j]==1)
            {E1<-matrix(E[i,],C[2],C[2]);
             E2<-matrix(t(E[j,]),C[2],C[2],byrow=TRUE);
             D_1<-(E1+E2)*D/2;
             D1=e*matrix(1,nrow=1,ncol=C[1]);
             D0=cbind(D_1,t(D1));D01=cbind(D1,c);
             D0=rbind(D0,D01);
             F1[,,k]<-D0;k<-k+1;    }
          }
     }

list(F1=F1);

}