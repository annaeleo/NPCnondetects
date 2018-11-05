# date: August 2018
# name procedure: Permutation tests for Nondetects
# @author: Eleonora Carrozzo (ref. paper: Multivariate permutation tests for two sample testing in presence of nondetects with application to microarray data, Arboretti, Bathke, Carrozzo, Pesarin, Salmaso, SMMR)

# The function performs a k-variate permutation test for two independent samples.

# 03.05.2017
# @author: Eleonora Carrozzo
# @param
# data : an array of (k + 1) columns. The first contain the groups and other the observations for each variable.
# stat : specifies the test statistics: "dm" for difference of means; "ad" for Anderson-Darling(*); "ad.Ruym" for Anderson-Darling with Ruy
# alt :  specifies if we want a two-sided test ("two.sided"), a one tailed test first_group > second_group ("greater") or a one tailed test first_group < second_group ("less)
# f.comb : combination function. On default is "F" for "Fisher" combination function

# @Output
# it prints on screen:
#  - the p-value of the test;
#  - who is the first and the second sample;
#  - the alternative;
# it returns a list of objects: 
#  - T.distribution is a vector or array of (B+1)-elements i.e. the permutation distribution of the test statistics. Element in the first row is the observed test statistics. 
#  - P.distribution is a vector or array of (B+1)-elements i.e. the permutation distribution of the p-value statistics. Element in the first row is the observed p-value statistics.

# *Note that Anderson-Darling test statistic is on distribution, so that the alternative X1<X2 corresponds to F(X1)>F(X2). The alternative to specify is that reffered to groups: i.e. if we have the alternative as F(X1)>F(X2) we have to specify alt="less".

perm.data<-function(data, stat="dm", B=1000, alt, f.comb = "F"){
  source("t2p.r")
  source("ad.r")
  source("ad_Ruym.r")
  source("mw.r")
  source("comb.r")
  
  #set.seed(1234)
  n=table(data[,1])  
  N = n[1] + n[2]
  
  k=ncol(data)-1
  
  T<-array(0,dim=c((B+1),k))
  
  Group1 = names(n)[1]; Group2 = names(n)[2]

  if(stat=="dm"){
    contr = rep(1/n,n)
    contr[-c(1:n[1])]=-contr[-c(1:n[1])]
    
    T[1,] = array(contr, dim=c(1,N)) %*% as.matrix(data[,-1])
    
    for(bb in 2:(B+1)){
      u <- sample(1:N,N)
      X.star <- data[u,-1]
      T[bb,] = array(contr, dim=c(1,N)) %*% as.matrix(X.star)
    }
    if(alt=="greater") {T=T ; alternative = paste("Alternative: differences between two groups are greater than 0")}
    if(alt=="less") {T=-T; alternative = paste("Alternative: differences between two groups are less than 0")}
    if(alt=="two.sided") {T= abs(T); alternative = paste("Alternative: differences between two groups are different from 0")}
  }#dm
  
  if(stat=="ad"){
    if(k>1) T[1,] <- apply(data[,-1], 2, function(z){ad(z,data[,1],alt)})
    if(k==1) T[1] <- ad(data[,-1],data[,1],alt)
    
    for(bb in 2:(B+1)){
      u <- sample(1:N,N)
      X.star <- data[u,-1]
      if(k>1) T[bb,] <- apply(X.star, 2, function(z){ad(z,data[,1],alt)}) 
      if(k==1) T[bb] <- ad(X.star,data[,1],alt)
    }
    if(alt=="greater") {alternative = paste("Alternative: distribution of ",Group1," is less than ", Group2)}
    if(alt=="less") {alternative = paste("Alternative: distribution of ",Group1," is greater than ", Group2)}
    if(alt=="two.sided") {alternative = paste("Alternative: distribution of ",Group1," and ", Group2, "differ")}
  }#ad
  
  if(stat=="ad_Ruym"){
    if(k>1) T[1,] <- apply(data[,-1], 2, function(z){ad.ruym(z,data[,1],alt)}) #Per il multivariato
    if(k==1) T[1] <- ad.ruym(data[,-1],data[,1],alt) #solo univariato
    
    for(bb in 2:(B+1)){
      u <- sample(1:N,N)
      X.star <- data[u,-1]
      if(k>1) T[bb,] <- apply(X.star, 2, function(z){ad.ruym(z,data[,1],alt)}) 
      if(k==1) T[bb] <- ad.ruym(X.star,data[,1],alt)
    }
    if(alt=="greater") {alternative = paste("Alternative: distribution of ",Group1," is less than ", Group2)}
    if(alt=="less") {alternative = paste("Alternative: distribution of ",Group1," is greater than ", Group2)}
    if(alt=="two.sided") {alternative = paste("Alternative: distribution of ",Group1," and ", Group2, "differ")}
  }
  
  if(stat=="mw"){
    if(k>1) T[1,] <- apply(data[,-1], 2, function(z){mw(z,label=data[,1],alt=alt)})
    if(k==1) T[1] <- mw(data[,-1],data[,1],alt=alt)
    
    for(bb in 2:(B+1)){
      u <- sample(1:N,N)
      X.star <- data[u,-1]
      if(k>1) T[bb,] <- apply(X.star, 2, function(z){mw(z,label=data[,1],alt=alt)})
      if(k==1) T[bb] <- mw(X.star,label=data[,1],alt=alt)
    }
    if(alt=="greater") {alternative = paste("Alternative: distribution of ",Group1," is less than ", Group2)}
    if(alt=="less") {alternative = paste("Alternative: distribution of ",Group1," is greater than ", Group2)}
    if(alt=="two.sided") {alternative = paste("Alternative: distribution of ",Group1," and ", Group2, "differ")}
  }
  
  P.partial = t2p(T)
  
  if(k>1){
    T.comb = comb(P.partial, f.comb)
    P.Glob = t2p(T.comb)
  }
  cat("### Permutation Two-Sample Test #####\n")
  cat("#####################################\n")
  cat("   \n")
  cat("Partial p.values =", P.partial[1,], "\n")
  cat("Global p.value =", P.Glob[1], "\n")
  cat("First sample: ", Group1,", Second sample: ", Group2, "\n")
  cat(alternative, "\n")
  
  if(k>1) return(list(T.distribution = T, P.partial = P.partial[1,], P.Glob = P.Glob[1]))
  if(k==1) return(list(P.partial = P.partial, T.partial = T))
}


                                                                # @author: Eleonora Carrozzo