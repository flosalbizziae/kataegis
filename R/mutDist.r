####################################################################
## Author: Xue Lin, Jian Li.
## Maintainer: Xue Lin <xue.lin@njmu.edu.cn>
## License: Artistic 2.0
## Part of the Kataegis package
## Reference: To be published
####################################################################

##Input:
### the output of the readMAF() or readVCF()

##Required by:
### interMutDist()

##Requires:
### none

#
mutDist<-function(data){ #data is the read in vcf file
	data<-data[order(data[,1],as.numeric(data[,2])),] #sort by chromosome position
	dist<-data.frame(matrix(ncol=3,nrow=length(data[,1])))
	j=0
	for(i in 2:length(data[,1])){
		j=j+1
		if(data[i,1]==data[i-1,1]){	
			dist[j,1]=gsub("chr","",data[i,1])
			x=as.numeric(data[i-1,2])
			y=as.numeric(data[i,2])
			dist[j,2]=x+((y-x)+1)/2
			dist[j,3]=as.numeric(data[i,2])-as.numeric(data[i-1,2])	
		}
		else{
			j=j-1
		}
	}
	colnames(dist)<-c("CHR","POS","DISTANCE")
	dist<-dist[complete.cases(dist),]
	dist<-dist[order(dist[,1],dist[,2],dist[,3]),]
	return(dist)
}