# date: August 2018
# name procedure: NPC 
# @author: Eleonora Carrozzo (ref. paper: Multivariate permutation tests for two sample testing in presence of nondetects with application to microarray data, Arboretti, Bathke, Carrozzo, Pesarin, Salmaso, SMMR)
#ref. Pesarin & Salmaso (2010)
comb<-function(pv,fcomb){
	if(fcomb=="F"){
		#cat('Combination Function: Fisher \n')
		fi<-apply(pv,1,function(x) -2*log(prod(x,na.rm=TRUE)))##Fisher
		}
	if(fcomb=="L"){
	#	cat('Combination Function: Liptak \n')
		fi<-apply(pv,1,function(x) sum(qnorm(1-x),na.rm=TRUE))##Liptak
		}
	if(fcomb=="T"){
	#	cat('Combination Function: Tippet \n')
		fi<-apply(pv,1,function(x) max((1-x),na.rm=TRUE)) ##Tippet
		}
		return(fi)
}