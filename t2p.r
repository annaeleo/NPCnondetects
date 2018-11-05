# date: August 2018
# name procedure: t2p (from statistics to p-value)
# @author: Eleonora Carrozzo (ref. paper: Multivariate permutation tests for two sample testing in presence of nondetects with application to microarray data, Arboretti, Bathke, Carrozzo, Pesarin, Salmaso, SMMR)
# ref. Pesarin & Salmaso (2010)
t2p<-function(T){

if(is.null(dim(T))){T<-array(T,dim=c(length(T),1))}
oth<-seq(1:length(dim(T)))[-1]

B<-dim(T)[1]-1							
p<-dim(T)[2]	
if(length(dim(T))==3){C<-dim(T)[3]}

rango<-function(x){
r=1-rank(x[-1],ties.method="min")/B+1/B
#return(c(mean(x[-1]>=x[1]),r))
return(c(mean(x>=x[1]),r))
}

P=apply(T,oth,rango)
return(P)
}
