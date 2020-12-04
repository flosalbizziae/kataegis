####################################################################
## Author: Xue Lin, Jian Li.
## Maintainer: Xue Lin <xue.lin@njmu.edu.cn>
## License: Artistic 2.0
## Part of the Kataegis package
## Reference: To be published
####################################################################

##Input:
###the kataegis regions called by kata() and the readVCF()/readMAF() output

##Required by:
### fcobarseplot()

##Requires:
### None

#To annotate the kataegis region with the mutation type for further visulization
#annoKata<-function(x, y, len=10000, n=10){
annoKata<-function(x, y, n=10){
  #x is the read in data in the first step with readVCF
  #y is the kataegis data
  #n is a integer, which is the bin-size, i.e. number of the division of the
    #flanking region totally, default is 10, which means that the upstream is
    #devided into 5 and the downstream is devided into 5, the merge controls
    #whether to merge all the dataframe in a list for plot.

  #check the parameters
	#if(!is.integer(len) & len>10e12){
	#	cat("\n\nError! The len sould be an integer which should be less than 10E+12!\n\n")
	#	stop()
	#}

	if(!is.integer(n) & n>20){
		cat("\n\nError! The n should be an interger which should be less than 20
		    and exactly divided by 2!\n\n")
		stop()
	}

  #if(0){#bloack all old codes
	#annotation start
	if(is.data.frame(y) && is.matrix(x)){
	  #foci divison method -2
	  #captallize the characters of neucleotides
	  x[,3]<-toupper(x[,3])
	  x[,4]<-toupper(x[,4])

	  n_pieces=floor(n/2)*2+1
	  tile<-matrix(rep(0,4*n_pieces),nrow=4,ncol=n_pieces,byrow=TRUE)
	  #same order TGCA
	  row.names(tile)<-c("T","G","C","A")
	  colid<-c()
	  for(i in seq_len(n_pieces)){#changed by XL on 2020-11-19
	    colid[i]=-(floor(n/2)+1)+i
	  }
	  colnames(tile)<-colid

	  for(i in seq_len(length(y[,1]))){#changed by XL on 2020-11-19
	    step_size=(as.numeric(y[i,4])-as.numeric(y[i,3]))/n_pieces
	    for(j in seq_len(n_pieces)){#changed by XL on 2020-11-19
	      all<-c()
	      all<-x[which(x[,1]==y[i,1] & as.numeric(x[,2])>as.numeric(y[i,3])+
	               (j-1)*step_size & as.numeric(x[,2])<=as.numeric(y[i,3])+j*step_size),]
	      if(is.vector(all)){
	        #only one mutation is returned
	        if(all[4]=="T"){tile[1,j]=tile[1,j]+1}
	        if(all[4]=="G"){tile[2,j]=tile[2,j]+1}
	        if(all[4]=="C"){tile[3,j]=tile[3,j]+1}
	        if(all[4]=="A"){tile[4,j]=tile[4,j]+1}
	      }else{
	        #more than one mutaiton returned or less than one mutation returned a data frame
	        tile[1,j]=tile[1,j]+length(all[which(all[,4]=="T"),1])
	        tile[2,j]=tile[2,j]+length(all[which(all[,4]=="G"),1])
	        tile[3,j]=tile[3,j]+length(all[which(all[,4]=="C"),1])
	        tile[4,j]=tile[4,j]+length(all[which(all[,4]=="A"),1])
	      }
	    }
	  }
	  tile<-sweep(tile,2,colSums(tile),"/")
	  return(tile)
	  ######################
	}else{
		cat("\nError! x and y should be data frame and matrix!\n")
		stop()
	}
}
