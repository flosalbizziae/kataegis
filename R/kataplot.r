####################################################################
## Author: Xue Lin, Jian Li.
## Maintainer: Xue Lin <xue.lin@njmu.edu.cn>
## License: Artistic 2.0
## Part of the Kataegis package
## Reference: To be published
####################################################################

##Input:
###the output of the kata and readVCF/MAF function

##Required by:
### None

##Requires:
### rainfall()
### fcobarplot()
### conbarplot()


#The kataplot will generate a rainfall plot for the mutations, a barplot of the nucleotide conversion of the mutations, and as well as a barplot for the flanking region component of the core kataegis foci.
kataplot<-function(x, y, z, n=20, fbar=TRUE, cbar=TRUE, rain=TRUE, colr=6, fbw=800, fbh=400, rw=1000, rh=500, assembly="hg19", name="all", type="png"){
	#x: is the matrix/list of matrix of the read in data from readVCF or readMAF
	#y: is the data.frame/list of data.frame of the kataegis data
  #z: a dataframe output from the interMutDist() function.
	#n: is the bin-size, i.e. number of the division of the flanking region totally, default is 20, which means that the upstream is devided into 10 and the downstream is devided into 10
	#fbar: the flanking region nucleotide spectra barplot
	#cbar: the core region mutation types barplot
	#rain: the rain-fall plot of all the raw data, and the poisson will also be tested
	#colr: the color set of the barplot and rainfall plot, default is the color palettes from grDevices, rainbow(6),which is "#FF0000" "#FFFF00" "#00FF00" "#00FFFF" "#0000FF" "#FF00FF". You can also use any R colors here.
	#fbw: the width of the flanking region nucleotides types barplot
	#fbh: the height of the flanking region nucleotides types barplot
	#rw: the width of the rainfall plot
	#rh: the height of the rainfall plot
  #assembly: the reference assembly version for analysis, default is human hg19, supported assembly incl. hg38, hg16-19, and mm7-10.
  #name: the sample name for the plot output, default is all.
  #type: the plot output format, default is png, other are jpeg, tiff, and eps

	#check parameters
  if(!assembly %in% c("hg38","hg19","hg18","hg17","hg16","mm7","mm8","mm9","mm10")){
    stop("assembly must be one of hg38, hg19, hg18, hg17 hg16 or mm10, mm9, mm8, mm7",call.=FALSE)
  }

  if(!is.integer(colr) & colr<4){
    cat("\n\nThe colr should be an interger no less than 4 which will be used for generating colorpallete for plots!\n\n")
    stop()
  }

  if(!is.integer(fbw) & fbw>3000){
    cat("\n\nThe fbw should be an integer no larger than 3000!\n\n")
    stop()
  }

  if(!is.integer(fbh) & fbh>3000){
    cat("\n\nThe fbh should be an integer no larger than 3000\n\n")
    stop()
  }

  if(!is.integer(rw) & rw>3000){
    cat("\nThe rw should be an integer no larger than 3000!\n")
    stop()
  }

  if(!is.integer(rh) & rh>3000){
    cat("\nThe rh should be an integer no larger than 3000!\n")
    stop()
  }

	if(!is.integer(n) & n>20){
		cat("\n\nError! The n should be an interger which should be less than 20!\n\n")
		stop()
	}

	if(!is.logical(fbar)){
		cat("\n\nError! The fbar should be a logical parameter!\n\n")
		stop()
	}

  if(!is.logical(cbar)){
    cat("\n\nError! The cbar should be a logical parameter!\n\n")
    stop()
  }

  if(!is.logical(rain)){
    cat("\n\nError! The rain should be a logical parameter!\n\n")
    stop()
  }

  #plot comes here
  if(is.data.frame(y) && is.matrix(x) && is.data.frame(z))	{
    if(length(y[1,])==6 && all(colnames(y)==c("chr", "arm", "start.pos", "end.pos", "n.mut", "dist.mean")) && length(x[1,])==5 && all(colnames(x)==c("#CHROM", "POS", "REF", "ALT", "FILTER")) &&length(z[1,])==3 && all(colnames(z)==c("CHR","POS","DISTANCE"))){
      cat("\nThe data-sets x and y both seems OK!\n" )
      #start plotting
      if(cbar){
        cat("\n\nStart the kataegis foci content bar plot...\n\n")
        conbarplot(x=x, name=name, cols=colr, type=type)

      }
      if(fbar){
        cat("\n\nStart the kataegis foci flanking region spectra barplot...\n\n")
        fcobarplot(x=x, y=y, name=name, type=type, col=colr, n=n, width=fbw, height=fbh)
        }
      if(rain){
        cat("\n\nStart the kataegis foci rainfall plot...\n\n")
        rainfall(x=x, y=z, name=name, type=type, assembly=assembly, col=colr, width=rw, height=rh)
      }
    }
    else{
      cat("\nError! The data-sets you input seems NOT OK!\n")
      cat("\nPlease make sure to input the right data-ests(refer to the documentations).\n")
      stop()
    }
  }
  else if(is.list(x) && is.list(y) && is.list(z)){
    if(length(y)==0){
      cat("\nThere is no kataegis called, thus no plot will be output.\n")
      stop()
    }else{
      if(length(y)<length(x)){
        sample_count=length(x)
        kata_sample_count=length(y)
        number_warn<-paste("\nYou have ", sample_count, ", but only ", kata_sample_count, " have kataegis called, which will be plotted.\n",sep='')
        cat(number_warn)
      }
      #start plotting here
      key<-names(y)
      #for(i in 1:length(key)){ #changed by XL on 2020-11-19
      for(i in seq_len(length(key))){
        origin<-x[[key[i]]]
        kataout<-y[[key[i]]]
        distout<-z[[key[i]]]

        if(cbar){
          cat("\n\nStart the kataegis foci content bar plot...\n\n")
          conbarplot(x=origin,name=key[i],cols=colr,type=type)
        }
        if(fbar){
          cat("\n\nStart the kataegis foci flanking region spectra barplot...\n\n")
          fcobarplot(x=origin,y=kataout,name=key[i],type=type,col=colr,n=n,width=fbw,height=fbh)
        }
        if(rain){
          cat("\n\nStart the kataegis foci rainfall plot...\n\n")
          rainfall(x=origin,y=distout,name=key[i],type=type,assembly=assembly,col=colr,width=rw,height=rh)
        }
      }
    }
  }else{
    cat("\nError! x and y should be data frame and matrix or both list!\n")
    stop()
  }
	#rain-fall plot: rainfall.r

	#the flanking and core region barplot: fcobarplot.r

  #the summarizing of mutation types barplot: conbarplot.r
}
