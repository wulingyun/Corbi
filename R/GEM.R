GEM<-function(ded,E)

{
emp0<-E;
C<-dim(E);
emp1<-ded*matrix(1,nrow=C[1],ncol=1);
emp<-cbind(emp0,emp1);

return(emp)

}