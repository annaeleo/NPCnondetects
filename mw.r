# date: August 2018
# name procedure: Mann-Withney test statistics
# @author: Eleonora Carrozzo (ref. paper: Multivariate permutation tests for two sample testing in presence of nondetects with application to microarray data, Arboretti, Bathke, Carrozzo, Pesarin, Salmaso, SMMR)
mw<-function(x,label, type="midrank", alt){
  nt=length(label)
  n = table(label)
 
  f = table(x)
  f.min = cumsum(f)-f
  
  k=dim(f)[1]
  
  F = (f/2 + f.min)
  
  if(type=="midrank"){
  x1 <- numeric(length(x))
  for(jj in names(F)){
    x1[x==names(F[jj])] <- F[jj]
  }
  
  R1 = x1[1:n[1]]
  R2 = x1[(n[1]+1):nt]
  
  if(alt=="greater") T_mw = sum(R1)
  if(alt=="less") T_mw = sum(R2)
  if(alt=="two.sided") T_mw = max(sum(R1)/n[1],sum(R2)/n[2]) #max(sum(F2/(F*(nt-F))),sum(F1/(F*(nt-F))))
  }
  
  if(type=="fmidrank"){
     W <- table(x,label)  
     F1.mw <- F*W[,1]
     F2.mw <- F*W[,2]
    
     if(alt=="greater") T_mw = sum(F1.mw)
     if(alt=="less") T_mw = sum(F2.mw)
     if(alt=="two.sided") T_mw = max(F1.mw, F2.mw)
  }
  return(T_mw)
}