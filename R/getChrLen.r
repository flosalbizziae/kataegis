####################################################################
## Author: Xue Lin, Jian Li.
## Maintainer: Xue Lin <xue.lin@njmu.edu.cn>
## License: Artistic 2.0
## Part of the Kataegis package
## Reference: To be published
####################################################################

##Input:
### the embeded genome ideograms

##Required by:
### ptest()

##Requires:
### None

#to calculate the chromosome length and return a dataframe
getChrLen<-function(assembly="hg19"){ #defualt genome assembly is hg19.
	if(assembly %in% c("hg38","hg19","hg18","hg17","hg16","mm7","mm8","mm9","mm10")){
		chrlen<-c()
		chrname<-c()
		#for( i in 1:length(unique(get(assembly)[,1]))){# changed by XL on 2020-11-19
		for( i in seq_len(length(unique(get(assembly)[,1])))){
			chr<-unique(get(assembly)[,1])[i]
			len<-max(get(assembly)[which(get(assembly)[,1]==chr),3])
			chrname[i]=as.character(chr)
			chrlen[i]=len
		}
		names(chrlen)<-chrname
		return(chrlen)
	}
	else{
		cat("\nError! The assembly should be hg38, hg19, hg18, hg17, hg16, mm10, mm9, mm8, or mm7.\n")
		stop()
	}
}
