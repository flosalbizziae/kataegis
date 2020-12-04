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
### None

##Requires:
### mutDist()

#
interMutDist<-function(data){
	if(is.matrix(data)){
		#the input is only a simple matrix which means the input is of one whole sample
		x<-mutDist(data)
		return(x) #return a data frame
	}
	else{
		#the input is a list holding multiple samples
		if(is.list(data)){ #the input is a list
			y<-list()
			#for(i in 1:length(data)){# changed by XL on 20202-11-19
			for(i in seq_len(length(data))){
				if(is.matrix(data[[i]]) && length(data[[i]])>0){ #which has more than 1 line
					x<-mutDist(data[[i]])
					if(length(x[,1])){	#the samples without inter-mutational distances are abandoned
						y[[names(data)[i]]]<-x
					}
				}
			}
			if(length(y)<length(data)){
				#the samples without the inter-mutational distances are abandoned.
				difference<-setdiff(names(data),names(y))
				cat("\nWarning! The following samples have too few locus for calculating inter-mutational distances, which will be abandonded in the subsequent analysis.\n")
				#for(j in 1:length(difference)){# changed by XL on 2020-11-19
				for(j in seq_len(length(difference))){
					cat(difference[j])
					cat("\n")
				}
				cat("\n")
			}
			return(y) #return a list
		}
		else{
			cat("\n\nError! The data is not a list or a matrix!\n\n")
			stop()
		}
	}
}
