####################################################################
## Author: Xue Lin, Jian Li.
## Maintainer: Xue Lin <xue.lin@njmu.edu.cn>
## License: Artistic 2.0
## Part of the Kataegis package
## Reference: To be published
####################################################################

##Input:
###the segmentations of the distances and call the kataegis

##Required by:
### kata()

##Requires:
### pcf()

callKataegis<-function(data,kmin=2,gamma=25,assembly="hg19",len=1000,nmut=6){ #the kmin for PCF is 2 and gamma is 25 which is tranned on a set of breast cancer data. And the len is the length of segment and the nmut is the number of mutations in the region.
	tmp<-pcf(data,kmin=kmin,gamma=gamma,assembly=assembly,verbose=FALSE)

	if(len=="Inf"){
		tmp<-tmp[,2:length(tmp[1,])]
		colnames(tmp)<-c("chr","arm","start.pos","end.pos","n.mut","dist.mean")
		return(tmp)
	}
	else if(len%%1==0){
		#tmp<-tmp[which((tmp[,5]-tmp[,4])<len && tmp[,6]>=nmut),2:length(tmp[1,])]
	  tmp<-tmp[which(tmp[,7]<=len & tmp[,6]>=nmut),2:length(tmp[1,])]
		colnames(tmp)<-c("chr","arm","start.pos","end.pos","n.mut","dist.mean")
		return(tmp)
	}
	else{
		cat("\nError!Please set the length of the segments with an integer  to identify the kataegis.\n")
		stop()
	}
}
