# date: August 2018
# name procedure: Anderson-Darling test statistics with Normalized distribution
# @author: Eleonora Carrozzo (ref. paper: Multivariate permutation tests for two sample testing in presence of nondetects with application to microarray data, Arboretti, Bathke, Carrozzo, Pesarin, Salmaso, SMMR)
ad.ruym<-function(x,label,alt){     
  nt=length(label)
  n = table(label)
  #Tad=array(0,dim=c(B+1))
  #Pv.ad=array(0,dim=c(B+1))
  
  W=table(x,label)
  k=dim(W)[1]
  
  f1 = W[,1]  
  f1.min = cumsum(f1)-f1
  f2 = W[,2]
  f2.min = cumsum(f2)-f2
  
  f=f1+f2
  f.min = cumsum(f)-f
  
  F1 = (f1/2 + f1.min)/n[1]
  
  F2 = (f2/2 + f2.min)/n[2]
 
  F = (f/2 + f.min)/nt
 
  #F =round(F,1)
  
  if(alt=="greater") T_ad =  sum((F2-F1)/(F*(1-F))^(1/2))# 
  if(alt=="less") T_ad = sum((F1-F2)/(F*(1-F))^(1/2))# 
  if(alt=="two.sided") T_ad = sum((F2-F1)^2/(F*(1-F))) 
  return(T_ad)
}