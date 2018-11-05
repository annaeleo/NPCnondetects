# June 2017
# date: August 2018
# name procedure: Anderson-Darling test statistics
# @author: Eleonora Carrozzo (ref. paper: Multivariate permutation tests for two sample testing in presence of nondetects with application to microarray data, Arboretti, Bathke, Carrozzo, Pesarin, Salmaso, SMMR)
# ref. Pesarin & Salmaso (2010)
# @param
# x : a vector of scores.
# label : a vector of labels specifing groups. 
# alt :  specifies if we want a two-sided test ("two.sided"), a one tailed test first_group > second_group ("greater") or a one tailed test first_group < second_group ("less). (*)

# @Output
# it prints on screen:
#  - the p-value of the test;
#  - who are the first and the second samples;
#  - the alternative;
# it returns a list of objects: 
#  - T.distribution is a vector or array of (B+1)-elements i.e. the permutation distribution of the test statistics. Element in the first row is the observed test statistics. 
#  - P.distribution is a vector or array of (B+1)-elements i.e. the permutation distribution of the p-value statistics. Element in the first row is the observed p-value statistics.

# *Note that Anderson-Darling test statistics is on distribution, so that the alternative X1<X2 corresponds to F(X1)>F(X2) (and viceversa). 
#  The alternative to specify is that reffered to groups: i.e. if we have the alternative as F(X1)>F(X2) we have to specify alt="less".

ad<-function(x,label, alt){     # 21Dic12 aggiornata con l'opzione ipotesi
nt=length(label)   
#Tad=array(0,dim=c(B+1))
#Pv.ad=array(0,dim=c(B+1))

W=table(x,label)
k=dim(W)[1]
f1=W[,1]
f2=W[,2]
f=f1+f2

F1=cumsum(f1)[1:(k-1)]
F2=cumsum(f2)[1:(k-1)]
F=F1+F2


if(alt=="greater") T_ad =  sum((F2-F1)/(F*(1-F))^(1/2))# 
if(alt=="less") T_ad = sum((F1-F2)/(F*(1-F))^(1/2))# 
if(alt=="two.sided") T_ad = sum((F2-F1)^2/(F*(1-F))) 
return(T_ad)
}

