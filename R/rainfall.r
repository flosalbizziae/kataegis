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
### interMutDist()

#To generate a rainfall plot of the conversion types of all the muations.
rainfall<-function(x, y, name="all", type="png", assembly="hg19", width=1000, height=500, col=6){
  #x: a matrix for plot, which is the format of the readVCF()/readMAF output.
  #y: a dataframe output from the interMutDist() function.
  #name: the sample name for the plot output, default is all.
  #type: the plot output format, default is png, other are jpeg, tiff, and eps
  #assembly: the reference assembly version for analysis, default is human hg19, supported assembly incl. hg38, hg16-19, and mm7-10.
  #width: the width of the output image, default is 1000.
  #height: the height of the output image, default is 500.
  #col: the color palatte number of the points.

  if(!assembly %in% c("hg38","hg19","hg18","hg17","hg16","mm7","mm8","mm9","mm10")){
    stop("assembly must be one of hg38, hg19, hg18, hg17 hg16 or mm10, mm9, mm8, mm7",call.=FALSE)
  }

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

  if(!is.numeric(col)&& col<6){
    cat("\n\nThe col should be a integaer!\n\n")
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
      }
      else{
        ######### change 2020Nov15 by Xue Lin #########
        #consider to display chromosomes or genomes automatically
        #calculate the length of each chromosome
        #chrlen<-getChrLen(assembly)
        #genomelen<-sum(chrlen)
        #plotlen<-(chrlen/genomelen)*width
        #indics<-suppressWarnings(sort(as.numeric(gsub("chr", '', names(chrlen))), na.last=T, index.return=T))
        #chr_sorted<-chrlen[indics$ix]
        #plotlen<-plotlen[indics$ix]
        #chr_cor_x<-c()
        #xlab_cor<-c()
        #for(i in 1:length(chr_sorted)){
        # chr_cor_x[i]=sum(plotlen[i-1:i])
        #  xlab_cor[i]=chr_cor_x[i]+plotlen[i]/2
        # }
        #names(chr_cor_x)<-names(chr_sorted)
        #pxcor<-c() #x coordinate indicates the genomics location of each mutation, and the y indicates the distance of the mutation to the previous mutation.
        #pcol<-c()
        #colpool<-rainbow(col+1)
        #x<-x[order(x[,1],as.numeric(x[,2])),] #sort according to the genomic locations
        #for(i in 1:length(y[,1])){
        #  if(x[i+1,1]==y[i,1]){
        #   ref=x[i+1,3]
        #    alt=x[i+1,4]
        #    #print(ref)
        #    #print(alt)
        #   #calculate the x and y coordinates
        #   chr<-paste("chr",x[i+1,1],sep='')
        #   pxcor[i]=chr_cor_x[chr]+(plotlen[chr]*y[i,2]/chrlen[chr])
        #  }
        ###############################################
        chrlen<-getChrLen(assembly)
        chrIDs<-paste("chr", x[,1],sep='')
        chrIDs<- chrIDs[!duplicated(chrIDs)]
        common_ids<-intersect(names(chrlen),chrIDs)
        chrlen<-chrlen[common_ids]
        indics<-suppressWarnings(sort(as.numeric(gsub("chr", '', common_ids)), na.last=TRUE, index.return=TRUE))
        chr_sorted<-chrlen[indics$ix]
        genomelen<-sum(chrlen)
        plotlen<-(chrlen/genomelen)*width
        plotlen<-plotlen[indics$ix]
        chr_cor_x<-c()
        xlab_cor<-c()
        for(i in seq_len(length(chr_sorted))){
          #chr_cor_x[i]=sum(plotlen[seq(i-1,i,by=1)])
          chr_cor_x[i]=sum(plotlen[i-1:i])
          xlab_cor[i]=chr_cor_x[i]+plotlen[i]/2
        }
        names(chr_cor_x)<-names(chr_sorted)
        #calculate the x and y coordinates
        pxcor<-c() #x coordinate indicates the genomics location of each mutation, and the y indicates the distance of the mutation to the previous mutation.
        pcol<-c()
        ntype<-c()
        colpool<-rainbow(col+1)
        x<-x[order(x[,1],as.numeric(x[,2])),] #sort according to the genomic locations
        j=1
        for(i in 2:length(x[,1])){
          if(as.character(x[i,1])!=as.character(x[i-1,1])){
            j=j+1
          }else{
            #cat(i)
            #cat("\t")
            #cat(y[i-j,1])
            #cat("\t")
            #cat(i-j)
            #cat("\t")
            #cat(x[i,1])
            #cat("\n")
            ref=x[i,3]
            alt=x[i,4]
            chr<-paste("chr",x[i,1],sep='')
            pxcor[i-j]=chr_cor_x[chr]+(plotlen[chr]*y[i-j,2]/chrlen[chr])
          }

          #assign the color
          if((ref=="C" && alt=="A")|(ref=="G" && alt=="T")){
            pcol[i]=colpool[1]
            ntype[1]="C>A"
          }else if((ref=="C" && alt=="G")|(ref=="G" && alt=="C")){
            pcol[i]=colpool[2]
            ntype[2]="C>G"
          }else if((ref=="C" && alt=="T")|(ref=="G" && alt=="A")){
            pcol[i]=colpool[3]
            ntype[3]="C>T"
          }else if((ref=="T" && alt=="A")|(ref=="A" && alt=="T")){
            pcol[i]=colpool[4]
            ntype[4]="T>A"
          }else if((ref=="T" && alt=="C")|(ref=="A" && alt=="G")){
            pcol[i]=colpool[5]
            ntype[5]="T>C"
          }else if((ref=="T" && alt=="G")|(ref=="A" && alt=="C")){
            pcol[i]=colpool[6]
            ntype[6]="T>G"
          }else{
            pcol[i]="grey" #grey color means you have indels other than point mutations or characters unrecoganizable
            ntype[7]="other"
          }
        }


        #pycor<-c()#y coordinates for points
        #for(i in 1:length(y[,3])){
          #pycor[i]=((y[i,3]-10**nchar(y[i,3]))/(10**(nchar(y[i,3])+1)-10**(nchar(y[i,3]))))*(height/ymaxexp)+(height/ymaxexp)*nchar(y[i,3])
          #pycor[i]=log10(y[i,3])
        #}
        pycor<-log10(y[,3])#y coordinates for points
        ytop=ceiling(max(pycor))
        yunit=height/ytop
        pycor<-pycor*yunit

        ylab_cor<-c()#y mark coordinates
        ylab<-c()
       # ymaxexp=nchar(max(y[,3]))+1
        #for(i in seq(0,ymaxexp)){
          #ylab[i+1]=10**i
          #ylab_cor[i+1]=(height/ymaxexp)*i
        #}
        ylab_cor<-seq(0,ytop)*yunit
        ylab<-10**seq(0,ytop)

        #create the image
        if(type=="png"){
          png(paste(name,"_rainfall.png",sep=''),width=width,height=height)
        }
        if(type=="jpeg"){
          jpeg(paste(name,"_rainfall.jpg",sep=''),width=width,height=height)
        }
        if(type=="tiff"){
          tiff(paste(name,"_rainfall.tiff",sep=''),width=width,height=height)
        }
        #if(type=="eps"){
         # eps(paste(name,"_rainfall.eps",sep=''),width=width,height=height)
        #}
        ###########changed by Xue Lin on 2020Nov15################
        #plot(x=chr_cor_x,seq(1,length(chr_cor_x)), ylim=c(0,height), type="n", xlab='', ylab="Inter-mutational Distance",  xaxt='n', yaxt='n')
        #plot(x=chr_cor_x,seq(1,24), ylim=c(0,height), type="n", xlab='', ylab="Inter-mutational Distance",  xaxt='n', yaxt='n')
        #plot(x=chr_cor_x,seq(1,24), ylim=c(summary(y[,3])[1],summary(y[,3])[3]), type="n", xlab='', ylab="Inter-mutational Distance",  xaxt='n', yaxt='n')
        #for(i in chr_cor_x){
        #  if(i>0){
        #    abline(v=i, col="grey", lty=2)
        #  }
        #}
        #axis(side=1, at=xlab_cor, labels=names(chr_sorted), tck=0, cex.axis=0.9)#x axis
        #axis(side=2, at=ylab_cor, labels=ylab, cex.axis=0.9)#y axis
        #points(x=pxcor,y=pycor, pch=20, col=pcol)#plot the points of the mutations
        ##########################################################
        #plot(x=chr_cor_x,seq(1,length(chr_cor_x)), ylim=c(0,height), type="n", xlab='', ylab="Inter-mutational Distance", xaxt='n', yaxt='n')
        #plot(x=seq(1,width), ylim=c(0,height), type="n", xlab='', ylab="Inter-mutational Distance", xaxt='n', yaxt='n')
        plot(x=seq(1,width), ylim=c(0,height), type="n", xlab='', ylab="Inter-mutational Distance", axes=FALSE)
        box()
       # if(length(chr_cor_x)>1){
          for(i in chr_cor_x){
            if(i>0){
             abline(v=i, col="grey", lty=2)
            }
          }
        #}else{
         # xlab_cor=(xlab_cor/width)*(length(chr_cor_x)+1)-1
        #  pxcor=(pxcor/width)*(length(chr_cor_x)+1)-1
        #}

        axis(side=1, at=xlab_cor, labels=names(chr_sorted), tck=0, cex.axis=0.9, lwd=0)
        axis(side=2, at=ylab_cor, labels=ylab, cex.axis=0.9)#y axis
        points(x=pxcor,y=pycor, pch=20, col=pcol)#plot the points of the mutations
        legend("top",ntype,fill=c(colpool[1:6],"grey"),horiz=TRUE,bty='n',xpd=TRUE,border=NA,inset=-0.1)
        dev.off()
      }
    }
  }
}
