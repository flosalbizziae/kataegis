####################################################################
## Author: Xue Lin, Jian Li.
## Maintainer: Xue Lin <xue.lin@njmu.edu.cn>
## License: Artistic 2.0
## Part of the Kataegis package
## Reference: To be published
####################################################################

##Input:
###the output of the readVCF()/readMAF functions

##Required by:
### kataplot()

##Requires:
### None

#To generate a barplot of the conversion types of all the muations.
conbarplot<-function(x, name="all",cols=6,type="png"){
  #x: a matrix for plot, which is the format of the readVCF()/readMAF output.
  #name: the sample name for the plot output, default is all.
  #type: the plot output format, default is png, other are jpeg, tiff, and eps

  if(!is.character(name)){
    cat("\nThe name should be a string as a file name to save your plot.\n")
    stop()
  }

  if(!is.character(type)){
    cat("\nThe type should be a character.\n")
    stop()
  }else{
    if(type=="png"|type=="jpeg"|type=="tiff"){
      if(!is.matrix(x)){
        cat("\nPlease input a matrix to plot.\n")
        stop()
      }
      else{
        tg=0
        tc=0
        ta=0
        ct=0
        cg=0
        ca=0
        #for(i in 1:length(x[,1])){#changed by XL on 2020-11-19
        for(i in seq_len(length(x[,1]))){
          if((x[i,3]=="T" && x[i,4]=="G") | (x[i,3]=="A" && x[i,4]=="C")){
            tg=tg+1
          }
          if((x[i,3]=="T" && x[i,4]=="C") | (x[i,3]=="A" && x[i,4]=="G")){
            tc=tc+1
          }
          if((x[i,3]=="T" && x[i,4]=="A") | (x[i,3]=="A" && x[i,4]=="T")){
            ta=ta+1
          }
          if((x[i,3]=="C" && x[i,4]=="T") | (x[i,3]=="G" && x[i,4]=="A")){
            ct=ct+1
          }
          if((x[i,3]=="C" && x[i,4]=="G") | (x[i,3]=="G" && x[i,4]=="C")){
            cg=cg+1
          }
          if((x[i,3]=="C" && x[i,4]=="A") | (x[i,3]=="G" && x[i,4]=="T")){
            ca=ca+1
          }
        }
        tmp<-c(tg,tc,ta,ct,cg,ca)
        names(tmp)<-c('T>G','T>C','T>A','C>T','C>G','C>A')
        if(type=="png"){
          png(paste(name,"_content_bar.png",sep=''))
        }
        if(type=="jpeg"){
          jpeg(paste(name,"_content_bar.jpg",sep=''))
        }
        if(type=="tiff"){
          tiff(paste(name,"_content_bar.tiff",sep=''))
        }
       #if(type=="eps"){
      #   eps(paste(name,"_content_bar.eps",sep=''))
      # }
        barplot(tmp,col=rainbow(cols+1),horiz=TRUE,border=NA,las=1,xpd=FALSE)
        graphics::box()
        #box()
        dev.off()
      }

    }else{
      #cat("\nThe type could be one of the png, jpeg, tiff or eps.\n")
      cat("\nThe type could be one of the png, jpeg, or tiff.\n")
      stop()
    }
  }


}

