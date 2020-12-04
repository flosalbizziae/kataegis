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
### callKataegis()

kata<-function(data,kmin=2,gamma=25,assembly="hg19",len=1000,nmut=6,verbose=TRUE){ #the kmin for PCF is 2 and gamma is 25 which is tranned on a set of breast cancer data. And the len is the length of segment and the nmut is the number of mutations in the region.
	if(is.data.frame(data)){
		#the input is only a simple matrix which means the input is of one whole sample
		y<-callKataegis(data,kmin=kmin,gamma=gamma,assembly=assembly,len=len,nmut=nmut)
		if(nrow(y)!=0){	#The y should not be empty
			return(y)
		}
		else{
			cat(sprintf("\nNo Kataegis called with parameters kmin=%i, gamma=%i, assembly=%s, len=%i, nmut=%i\n",kmin,gamma,assembly,len,nmut))
			stop()
		}


	}
	else{
		if(is.list(data)){
			#the input is a list holding multiple samples
			y<-list()
			#for(i in 1:length(data)){ #changed by XL on 2020-11-19
			for(i in seq_len(length(data))){
				#iterate among samples
				if(verbose){
					cat("Sample ")
					cat(names(data)[i])
					cat(" is being processed.\n")
				}
				x<-callKataegis(data[[i]],kmin=kmin,gamma=gamma,assembly=assembly,len=len,nmut=nmut)
				if(nrow(x)!=0){
					y[[names(data)[i]]]<-x
				}
				else{
					if(verbose){
						cat("Sample ")
						cat(names(data)[i])
						cat(" has no kataegis called, which will be excluded in the final result.\n")
					}
				}
			}
			return(y)
		}
		else{
			cat("\nError! Your input is not an appropriate type for kataegis calling, either a data frame or a list of data frame is acceptable.\n")
			stop()
		}
	}
}
