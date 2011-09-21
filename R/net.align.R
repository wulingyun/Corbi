net.align<-function(quenet,net,nodefea,qtype=2,delta.d=1e-10,delta.c=0.5,delta.e=1,delta.s=1,isnormalize=0,multitest=1)
{
# paramters explain
# quenet: query network file name
# net: target network file name
# nodefea: node feature function values
# qtype: the query network type
# delta.d: the parameter \Delta_d
# delta.c: the parameter \Delta_c
# delta.e: the parameter \Delta_e
# delta.s: the parameter \Delta_s
# isnormalize: whether carry out the normalized operation to the node feature values or not? default not
# multitest: the result files' number; only using in many experiments

# read the files
netalign<-trans(quenet,net,nodefea);
# netalign: 
# netalign$X query subnetwork name; netalign$D_X the adjacency matrix of query network; netalign$N target network names; netalign$G the adjacency matrix of the target network; netalign$E the node similar matrix

# is or not normalize the node feature function
if(isnormalize>0){a=min(netalign$E);b=max(netalign$E);netalign$E=(netalign$E-a)/(b-a);}
if(max(netalign$E)<=0){a=min(netalign$E);b=max(netalign$E);netalign$E=(netalign$E-a)/(b-a);}

# compute the distance matrix for the target network
D_G<-dis(netalign$G,dim(netalign$G)[1]);

ded=delta.d; # the similar value to the deletion node Vd
# simplification the target network
netalign1=reduf(netalign$E,netalign$N,D_G,delta.d);



emp=GEM(ded,netalign1$E1);

# compute the transition probability
D_X=netalign$D_X;
s<-D_X-t(D_X);
s1<-sum(s[s>0]);
s2<-sum(s[s<0]);
if(s1==0 & s2==0){F1<-tranf1(netalign1$g,emp,D_X,delta.s,delta.e,delta.c);}
else{F1<-tranf2(netalign1$g,emp,D_X,delta.s,delta.e,delta.c);}

D_X<-D_X+t(D_X);
D_X[D_X>1]<-1;

Y<-dyf(netalign$X,emp,F1$F1,netalign1$name,netalign1$g,D_X,qtype);

resultxt(netalign$X,Y,D_X,netalign1$name,netalign1$g,multitest)

}