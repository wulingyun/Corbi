trans<-function(query,target,feature){
# change the text to CRF method
# query: query network text, two nodes at least
# target: target network text, two nodes at least
# feature: the node feature score text

#=============================================================

  queryname='';k=0;I=c();J=c();S<-1;
  con <- file(query, "r")
  line=readLines(con,n=1);
  while( length(line) != 0 ) {            
        subc<-strsplit(line," ");
          ik0=match(subc[[1]][1],queryname,nomatch=0);
          if(ik0){ik=ik0;}
          else{
              k=k+1; queryname[k]=subc[[1]][1]; ik=k;}

        for (i in 2:length(subc[[1]])){
          J[S]=match(subc[[1]][i],queryname,nomatch=0);
          if(J[S]==0){
              k=k+1; queryname[k]=subc[[1]][i]; J[S]=k;}

          I[S]<-ik;S=S+1;
         }
        
        line=readLines(con,n=1); 
  } # end while

  close(con)
  dquery<-matrix(0,k,k);dquery[cbind(I,J)]<-1;
  
#=============================================================

  targetname='';k=0;I=c();J=c();S<-1;
  con <- file(target, "r")
  line=readLines(con,n=1);
  while( length(line) != 0 ) {        
        subc<-strsplit(line," ");
          ik0=match(subc[[1]][1],targetname,nomatch=0);
          if(ik0){ik=ik0;}
          else{
              k=k+1; targetname[k]=subc[[1]][1]; ik=k;}

        for (i in 2:length(subc[[1]])){
          J[S]=match(subc[[1]][i],targetname,nomatch=0);
          if(J[S]==0){
          k=k+1; targetname[k]=subc[[1]][i]; J[S]=k;}

          I[S]<-ik;S=S+1;
         }

          line=readLines(con,n=1); 
  } # end while

  close(con)
  dtarget<-matrix(0,k,k);dtarget[cbind(I,J)]<-1;

#=============================================================

  fea=matrix(0,dim(dquery)[1],dim(dtarget)[1]);
  #I<-c();J<-c();K<-c();S<-1;
  con <- file(feature, "r")
  line=readLines(con,n=1);
  while( length(line) != 0 ) {
        line<-gsub("\t"," ",line);            
        subc<-strsplit(line," ");
        fea[match(subc[[1]][1],queryname,nomatch=0),match(subc[[1]][2],targetname,nomatch=0)]<-as.numeric(subc[[1]][3]); 
        #k=1;nameone='';
        #for (i in 1:length(subc[[1]])){if(sum(match(subc[[1]][i],"",nomatch=0))==0){nameone[k]=subc[[1]][i];k=k+1;}}
        #K[S]=as.numeric(nameone[3]); 
        #J[S]=match(nameone[2],targetname,nomatch=0); I[S]=match(nameone[1],queryname,nomatch=0);
        #fea[I[S],J[S]]<-K[S];
        #S=S+1;
        
        line=readLines(con,n=1); 
  } # end while
  
  
  close(con)

#=============================================================


  list(X=queryname,D_X=dquery,N=targetname,G=dtarget,E=fea)
}