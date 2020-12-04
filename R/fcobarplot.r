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
### annoKata()

#To generate a barplot of the conversion types of all the muations.
fcobarplot<-function(x, y, name="all", col=6, type="png", n=20, width=800, height=400){
  #x: a matrix for plot, which is the format of the readVCF()/readMAF output.
  #y: the kataegis data
  #name: the sample name for the plot output, default is all.
  #type: the plot output format, default is png, other are jpeg, tiff, and eps
  #col: the color palatte number of different types of nucleotides
  #n: a integer, which is the bin-size, i.e. number of the division of the flanking region totally, default is 10, which means that the upstream is devided into 5 and the downstream is devided into 5, the merge controls whether to merge all the dataframe in a list for plot.
  #width: the width of the output image, default is 200.
  #height: the height of the output image, default is 100.
  if(!is.integer(width) & width>3000){
    cat("\nThe width should be an integer no larger than 3000!\n")
    stop()
  }
  if(!is.integer(height) & height>3000){
    cat("\nThe height should be an integer no larger than 3000!\n")
    stop()
  }
  if(!is.character(name)){
    cat("\nThe name should be a string as a file name to save your plot.\n")
    stop()
  }
  if(!is.integer(n) & n>50){
    cat("\n\nError! The n should be an interger which should be less than 50 and exactly divided by 2!\n\n")
    stop()
  }
  if(!is.numeric(col)){
    cat("\n\nThe col should be a integer!\n\n")
    stop()
  }
  if(!is.character(type)){
    cat("\nThe type should be a character.\n")
    stop()
  }else{
    if(type=="png"|type=="jpeg"|type=="tiff"){
      if(!is.matrix(x)){
        cat("\nPlease input a matrix to x.\n")
        stop()
      }
      if(!is.data.frame(y)){
        cat("\nPlease input a dataframe to y.\n")
        stop()
      }else{
        #plot cames here
        result<-annoKata(x,y,n=n)
        colpool<-rainbow(col+1)
        ntcol<-c()
        ntype<-row.names(result)
        #for(i in 1:length(ntype)){#changed by XL on 2020-11-19
        for(i in seq_len(length(ntype))){
          if(ntype[i]=="T"){ntcol[i]=colpool[i]}
          if(ntype[i]=="G"){ntcol[i]=colpool[i]}
          if(ntype[i]=="C"){ntcol[i]=colpool[i]}
          if(ntype[i]=="A"){ntcol[i]=colpool[i]}
        }
        if(type=="png"){
          png(paste(name,"_spectra.png",sep=''),width=width,height=height)
        }
        if(type=="jpeg"){
          jpeg(paste(name,"_spectra.jpg",sep=''),width=width,height=height)
        }
        if(type=="tiff"){
          tiff(paste(name,"_spectra.tiff",sep=''),width=width,height=height)
        }
        #if(type=="eps"){
         # eps(paste(name,"_spectra.eps",sep=''),width=width,height=height)
        #}
        #par(xpd=T, mar=par()$mar+c(0,0,1,0))
        par(xpd=TRUE, mar=c(4,4,6,4))
        barplot(result,space=0,col=ntcol)
        legend("top",ntype,fill=ntcol,horiz=TRUE,bty='n',xpd=TRUE,border=NA,inset=-0.1)
        dev.off()
      }
    }else{
      #cat("\nAccepted figure output types are: png, jpeg, tiff and eps.\n")
      cat("\nAccepted figure output types are: png, jpeg, and tiff.\n")
      stop()
    }
  }
}





