
####################################################################
## Author: Xue Lin and Jian Li
## Maintainer: Xue Lin xue.lin@njmu.edu.cn
## License: Artistic 2.0
## Part of the Kataegis package
## Reference: Herited from the copynumber package pcf fucntion
## Reference: Nilsen and Liest?l et al. (2012), BMC Genomics
####################################################################


##Required by:
### callKataegis.r


##Requires:
### None


## Main function for pcf-analysis to be called by the user

pcf <- function(data,pos.unit="bp",arms=NULL,Y=NULL,kmin=5,gamma=40,normalize=TRUE,fast=TRUE,assembly="hg19",digits=4,return.est=FALSE,save.res=FALSE,file.names=NULL,verbose=TRUE){

  #Check pos.unit input:
  if(!pos.unit %in% c("bp","kbp","mbp")){
    stop("pos.unit must be one of bp, kbp and mbp",call.=FALSE)
  }

  #Check assembly input:
  #if(!assembly %in% c("hg19","hg18","hg17","hg16","mm7","mm8","mm9")){ #revised by XL
  if(!assembly %in% c("hg38","hg19","hg18","hg17","hg16","mm7","mm8","mm9","mm10")){
    stop("assembly must be one of hg38, hg19, hg18, hg17 hg16 or mm10, mm9, mm8, mm7",call.=FALSE)
  }

  #Is data a file:
  isfile.data <- is.character(data)#changed by XL 2020-11-19

  #Check data input:
  if(!isfile.data){
    #Input could come from winsorize and thus be a list; check and possibly retrieve data frame wins.data
    ################pullOutContent################
    #revised 2019-12-01
    pullOutContent <- function(res,what="segments"){

      #check if input is data frame or list
      if(!is.data.frame(res)){
        #res could either be a list containing the two segmentation elements segments and estimates, a list containing the two winsorize elements wins.data and wins.outliers, a list containing several segments as data frames, or a list containing several lists with segmentation results
        if("segments" %in% names(res)){
          #can assume that the list contains the output from one of the segmentation algorithms and has names segments and estimates
          #pick out the desired data frame depending on input in what:
          if(what=="estimates"){
            if("logR_estimates" %in% names(res)){
              #Segmentation results come from aspcf which returns a different name for estimates
              res <- res$logR_estimates
            }else{
              res <- res$estimates
            }
          }else if(what=="segments"){
            res <- res$segments
          }
        }else if("wins.data" %in% names(res)){
          #can assume that the list contains output from winsorize and has names wins.data and wins.outliers
          #pick out the desired data frame depending on input in what:
          if(what %in% c("wins.data","estimates")){
            res <- res$wins.data
          }else if(what=="wins.outliers"){
            res <- res$wins.outliers
          }
        }else{
          #assume that the list contains different segmentation results
          #if each element in the list is itself a list containing segments and estimates need to pull out the segments (functions that take estimates as input does not have the option of specifying list of estimates)
          nSeg <- length(res)
          for(l in seq_len(nSeg)){
            if(!is.data.frame(res[[l]])){
              #pull out the segments data frame
              res[[l]] <- res[[l]]$segments
            }
          }
        }
      }
      #If a data frame, res should be already be a data frame with the correct information, and does not need modification

      return(res)

    }
    ################pullOutContent################

    data <- pullOutContent(data,what="wins.data")
    stopifnot(ncol(data)>=3)  #something is missing in input data
    #Extract information from data:
    chrom <- data[,1]
    position <- data[,2]
    nSample <- ncol(data)-2
    sampleid <- colnames(data)[-c(seq_len(2))] #sampleid <- colnames(data)[-c(1:2)]

  }else{
    #data is a datafile which contains data
    f <- file(data,"r")  #open file connection
    head <- scan(f,nlines=1,what="character",quiet=TRUE,sep="\t") #Read header
    if(length(head)<3){
      stop("Data in file must have at least 3 columns",call.=FALSE)
    }
    sampleid <- head[-c(seq_len(2))]#sampleid <- head[-c(1:2)]
    nSample <- length(sampleid)

    #Read just the two first columns to get chrom and pos
    chrom.pos <- read.table(file=data,sep="\t",header=TRUE,colClasses=c(rep(NA,2),rep("NULL",nSample)),as.is=TRUE)  #chromosomes could be character or numeric
    chrom <- chrom.pos[,1]
    position <- chrom.pos[,2]
  }


  #Make sure chrom is not factor:
  if(is.factor(chrom)){
    #If chrom is factor; convert to character
    chrom <- as.character(chrom)
  }

  #Make sure chromosomes are numeric (replace X and Y by 23 and 24)
  #####################numericChrom##########################
  ##revised 2019-12-01
  numericChrom <- function(chrom){
   if(!is.numeric(chrom)){
      if(is.factor(chrom)){
        #If chrom is factor; need to convert to character first
        chrom <- as.character(chrom)
      }
      #Replace X by 23:
      chrx <- c(which(chrom=="x"),which(chrom=="X"))
      chrom[chrx] <- 23
      #Replace Y by 24
      chry <- c(which(chrom=="y"),which(chrom=="Y"))
      chrom[chry] <- 24

      chrom <- as.numeric(chrom)
    }
    return(chrom)
  }
  #####################numericChrom##########################

  num.chrom <- numericChrom(chrom)
  nProbe <- length(num.chrom)

  #Make sure position is numeric:
  if(!is.numeric(position)){
    stop("input in data column 2 (posistions) must be numeric",call.=FALSE)
  }
  #Get character arms:
  #######################getAmrms########################
  #revised 2019-12-01
  getArms <- function(chrom, pos, pos.unit="bp", cyto.data){

  #Make sure chromosomes are numeric:


  chrom <- numericChrom(chrom)

	nProbe <- length(chrom)
	chrom.list <- unique(chrom)
	nChrom <- length(chrom.list)

	#Get local stopping posistions for each p-arm and each chromosome from cytoband data
	######################getArmandChromStop##########################
	##revised 2019-12-01
	getArmandChromStop <- function(cyto.data, unit){

	#Sort cyto.data by chromosome number; let be represented by X=23 and Y=24:
	chrom <- cyto.data[,1]
	use.chrom <- gsub("chr","",chrom)  #Remove 'chr' from chrom-strings
	use.chrom[use.chrom=="X"] <- "23"	#Replace X by 23
	use.chrom[use.chrom=="Y"] <- "24"	#Replace Y by 24
	num.chrom <- as.numeric(use.chrom)	#Convert to numeric chromosomes

	#Order such that chromosomes are in increasing order from 1:24:
	ord.chrom <- order(num.chrom)
	cyto.data <- cyto.data[ord.chrom,,drop=FALSE]

	#Get chromosome stopping positions:
	chrom <- cyto.data[,1]
	chrom.stop <- which(chrom[seq_len(length(chrom))-1]!=chrom[seq(2,length(chrom),by=1)]) #chrom.stop <- which(chrom[1:length(chrom)-1]!=chrom[2:length(chrom)])
	chrom.stop <- c(chrom.stop,length(chrom))  #include last chromstop as well

	#Get p-arm stopping positions:
	arm.char <- substring(cyto.data[,4],1,1)   #Retrive first character in name which identifies p and q arms
	arm.stop <- which(arm.char[seq_len(length(arm.char))-1]!=arm.char[seq(2,length(arm.char),by=1)])#arm.stop <- which(arm.char[1:length(arm.char)-1]!=arm.char[2:length(arm.char)])
	p.stop <- arm.stop[-which(arm.stop%in%chrom.stop)]  #Remove qstops

	pos.chromstop <- cyto.data[chrom.stop,3]  #Local stopping position for each chromosome
	pos.pstop <- cyto.data[p.stop,3]		#Local stopping position for each p-arm

	#Factor used to convert positions into desired unit
	f <- switch(unit,
		bp = 1,
		kbp = 10^(-3),
		mbp = 10^(-6))

	return(list(pstop=pos.pstop*f,chromstop=pos.chromstop*f))

    }
	######################getArmandChromStop##########################

	l <- getArmandChromStop(cyto.data=cyto.data,unit=pos.unit)
	pStop <- l$pstop
	chromStop <- l$chromstop

	#Intitialize
	arms <- rep(NA,nProbe)

	for(i in seq_len(nChrom)){
		#Find corresponding arm numbers:
		c <- chrom.list[i]
		ind.c <- which(chrom==c)

		arms[ind.c] <- "q"
		p.arm <- ind.c[pos[ind.c]<=pStop[c]]
		arms[p.arm] <- "p"   #p-arm

	}

	return(arms)
    }
   #######################getAmrms########################

	if(is.null(arms)){
    arms <- getArms(num.chrom,position,pos.unit,get(assembly))
	}else{
    stopifnot(length(arms)==nProbe)
	}

	#Translate to numeric arms:
	##################numericArms#######################
	#revised 2019-12-01
	numericArms <- function(chrom,char.arms){
        p.arm <- which(char.arms=="p")
        q.arm <- which(char.arms=="q")
        arms <- rep(NA,length(char.arms))
        arms[p.arm] <- chrom[p.arm]*2-1
        arms[q.arm] <- chrom[q.arm]*2

        return(arms)
    }
	##################numericArms#######################
	num.arms <- numericArms(num.chrom,arms)
	#Unique arms:
	arm.list <- unique(num.arms)
	nArm <- length(arm.list)

  #Check Y input:
  if(!is.null(Y)){
    stopifnot(class(Y)%in%c("matrix","data.frame","character"))
    isfile.Y <- is.character(Y) #changed by XL on 2020-11-19
    if(!isfile.Y){
      ncol.Y <- ncol(Y)
      nrow.Y <- nrow(Y)
    }else{
      f.y <- file(Y,"r")
      ncol.Y <- length(scan(f.y,nlines=1,what="character",quiet=TRUE,sep="\t"))
      nrow.Y <- nrow(read.table(file=Y,sep="\t",header=TRUE,colClasses=c(NA,rep("NULL",ncol.Y-1)),as.is=TRUE))
    }
    if(nrow.Y!=nProbe || ncol.Y!=nSample+2){
      stop("Input Y does not represent the same number of probes and samples as found in input data",call.=FALSE)
    }
  }

  #save user's gamma
	gamma0 <- gamma

	sd <- rep(1,nSample) #sd is only used if normalize=TRUE, and then these values are replaced by MAD-sd
	#If number of probes in entire data set is less than 100K, the MAD sd-estimate is calculated using all obs for every sample
  #Only required if normalize=T
  if(nProbe<100000 && normalize){
    #calculate MAD-sd for each sample:
    for(j in seq_len(nSample)){
      if(!isfile.data){
        sample.data <- data[,j+2]
      }else{
        cc <- rep("NULL",nSample+2)
        cc[j+2] <- "numeric"
        #only read data for the j'th sample
        sample.data <- read.table(file=data,sep="\t",header=TRUE,colClasses=cc)[,1]
      }

      ###############getMad#################
      #revised 2019-12-01
      getMad <- function(x,k=25){

	  #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
	  x <- x[x!=0]

      #Calculate runMedian
      medianFilter <- function(x,k){
      n <- length(x)
      filtWidth <- 2*k + 1

      #Make sure filtWidth does not exceed n
      if(filtWidth > n){
        if(n==0){
          filtWidth <- 1
        }else if(n%%2 == 0){
          #runmed requires filtWidth to be odd, ensure this:
          filtWidth <- n - 1
        }else{
          filtWidth <- n
        }
      }
      ##
      runMedian <- runmed(x,k=filtWidth,endrule="median")

      return(runMedian)

     }
     runMedian <- medianFilter(x,k)

     dif <- x-runMedian
     SD <- mad(dif)

	 return(SD)
    }
    ###############getMad#################

      sd[j] <- getMad(sample.data[!is.na(sample.data)],k=25)   #Take out missing values before calculating mad
		}
  }#endif


	#Initialize
	pcf.names <- c("chrom","pos",sampleid)
	seg.names <- c("sampleID","chrom","arm","start.pos","end.pos","n.probes","mean")

  segments <- data.frame(matrix(nrow=0,ncol=7))
  colnames(segments) <- seg.names
  if(return.est){
    pcf.est <- matrix(nrow=0,ncol=nSample)
	}
	if(save.res){
    if(is.null(file.names)){
      #Create directory where results are to be saved
      dir.res <- "pcf_results"
      if(!dir.res %in% dir()){
        dir.create(dir.res)
      }
      file.names <- c(paste(dir.res,"/","estimates.txt",sep=""),paste(dir.res,"/","segments.txt",sep=""))
    }else{
      #Check that file.names is the correct length
      if(length(file.names)<2){
        stop("'file.names' must be of length 2", call.=FALSE)
      }
    }
  }

  #estimates must be returned from routines if return.est or save.res
  yest <- any(return.est,save.res)

  #run PCF separately on each chromosomearm:
  for(c in seq_len(nArm)){
    probe.c <- which(num.arms==arm.list[c])
    pos.c <- position[probe.c]
    this.arm <- unique(arms[probe.c])
    this.chrom <- unique(chrom[probe.c])
    #Empty result matrix and dataframe for this arm
    segments.c <- data.frame(matrix(nrow=0,ncol=7))
    colnames(segments.c) <- seg.names
    if(yest){
      pcf.est.c <-  matrix(nrow=length(probe.c),ncol=0)
    }

    #Get data for this arm
    if(!isfile.data){
        arm.data <- data[probe.c,-c(seq_len(2)),drop=FALSE]#arm.data <- data[probe.c,-c(1:2),drop=FALSE]
    }else{
      #Read data for this arm from file; since f is a opened connection, the reading will start on the next line which has not already been read
      #two first columns are skipped
      arm.data <- read.table(f,nrows=length(probe.c),sep="\t",colClasses=c(rep("NULL",2),rep("numeric",nSample)))
    }

    #Make sure data is numeric:
    if(any(!sapply(arm.data,is.numeric))){
      stop("input in data columns 3 and onwards (copy numbers) must be numeric",call.=FALSE)
    }
    #Get Y-values for this arm
    if(!is.null(Y)){
      if(!isfile.Y){
          arm.Y <- Y[probe.c,-c(seq_len(2)),drop=FALSE]#arm.Y <- Y[probe.c,-c(1:2),drop=FALSE]
      }else{
        arm.Y <- read.table(f.y,nrows=length(probe.c),sep="\t",colClasses=c(rep("NULL",2),rep("numeric",nSample)))
      }
      #Make sure Y is numeric:
      if(any(!sapply(arm.Y,is.numeric))){
        stop("input in Y columns 3 and onwards (copy numbers) must be numeric",call.=FALSE)
      }
    }

    #Run PCF separately for each sample:
    for(i in seq_len(nSample)){
      sample.data <- arm.data[,i]

			#Remove probes with missing obs; Only run pcf on non-missing values
			obs <- !is.na(sample.data)
			obs.data <- sample.data[obs]

		  if(yest){
			 #Initialize:
			 yhat <- rep(NA,length(probe.c))
			}

      if(length(obs.data)>0){  ##Make sure there are observations on this arm, this sample! If not, estimates are left NA as well

  			#If number of probes in entire data set is >= 100K, the MAD sd-estimate is calculated using obs in this arm for this sample.
        #Only required if normalize=T
        if(nProbe>=100000 && normalize){
          sd[i] <- getMad(obs.data,k=25)
        }

        #Scale gamma by variance if normalize is TRUE
        use.gamma <- gamma0
        if(normalize){
          use.gamma <- gamma0*(sd[i])^2
        }

    		#Must check that sd!=0 and sd!!=NA -> use.gamma=0/NA. If not, simply calculate mean of observations
  			if(use.gamma==0 || is.na(use.gamma)){
  			  if(yest){
            res <- list(Lengde=length(obs.data),sta=1,mean=mean(obs.data),nIntervals=1,yhat=rep(mean(obs.data)))
          }else{
            res <- list(Lengde=length(obs.data),sta=1,mean=mean(obs.data),nIntervals=1)
          }

  			}else{
    			# Compute piecewise constant fit
    			#run fast approximate PCF if fast=TRUE and number of probes>400, or exact PCF otherwise
    			if(!fast || length(obs.data)<400){
    				#Exact PCF:
    			 	res <- exactPcf(y=obs.data,kmin=kmin,gamma=use.gamma,yest=yest)
    			}else{
            #Run fast PCF:
            res <- selectFastPcf(x=obs.data,kmin=kmin,gamma=use.gamma,yest=yest)
    			}#endif


  			}#endif

  			#Retrieve segment info from results:
  			seg.start <- res$sta
			  seg.stop <- c(seg.start[-1]-1,length(obs.data))
			  seg.npos <- res$Lengde
  			seg.mean <- res$mean
  			nSeg <- res$nIntervals
        if(yest){
          yhat[obs] <- res$yhat
        }
        #Find genomic positions for start and stop of each segment:
			  pos.start <- pos.c[obs][seg.start]        #pos.c[seg.start]   NOTE: THIS ERROR WAS CORRECTED 11/3-2013
			  pos.stop <- pos.c[obs][seg.stop]          #pos.c[seg.stop]    NOTE: THIS ERROR WAS CORRECTED 11/3-2013

        #Handle missing values:
  			if(any(!obs)){
  			  #first find nearest non-missing neighbour for missing probes:
  		  ##################findNN################
  		  #revised 2019-12-01
  		  findNN <- function(pos,obs){

  		  	ind.obs <- which(obs)
	  		  pos.obs <- pos[obs]
	  		  pos.na <- pos[!obs]

  		  	#Find distances between position of missing obs and positions of non-missing obs
	  		  d <- sapply(pos.na,FUN="-",y=pos.obs)
	  		  if(!is.matrix(d)){
    		    d <- matrix(d,nrow=length(pos.obs),ncol=length(pos.na))
	  		  }
	  		  d <- abs(d)
	  		  nn <- apply(d,2,which.min)  #find which observed obs is the nearest
	  		  nn <- ind.obs[nn]  #get index for the nearest

	  		  return(nn)

  		  }

  		  ##################findNN################
          nn <- findNN(pos=pos.c,obs=obs)

          #Include probes with missing values in segments where their nearest neighbour probes are located
          #new.res <- handleMissing(nn=nn,pos=pos.c,obs=obs,pos.start=pos.start,pos.stop=pos.stop,seg.npos=seg.npos)
          #pos.start <- new.res$pos.start
          #pos.stop <- new.res$pos.stop
          #seg.npos <- new.res$seg.npos

          if(yest){
            yhat[!obs] <- yhat[nn]
          }
  			}

			}else{
			  warning(paste("pcf is not run for sample ",i," on chromosome arm ",this.chrom,this.arm," because all observations are missing. NA is returned.",sep=""),immediate.=TRUE,call.=FALSE)
			  seg.start <- 1
			  seg.stop <- length(pos.c)
		    pos.start <- pos.c[seg.start]
		    pos.stop <- pos.c[seg.stop]
		    nSeg <- 1
		    seg.mean <- NA
		    seg.npos <- length(pos.c)
		  }


      #Check if mean segment-value should be replaced by Y-values (possibly wins.data):
			if(!is.null(Y)){
				seg.mean <- rep(NA,nSeg)
        #Use observed data to calculate segment mean (recommended)
				sample.y <- arm.Y[,i]
				#Make sure data is numeric:
        for(s in seq_len(nSeg)){
					seg.mean[s] <- mean(sample.y[seg.start[s]:seg.stop[s]],na.rm=TRUE)
				}
			}
			#Round:
			seg.mean <- round(seg.mean,digits=digits)
			#Create table with relevant segment-information
			seg.arm <- rep(this.arm,nSeg)
      seg.chrom <- rep(this.chrom,nSeg)

			#Add results for this sample to results for other samples in data frame:
      seg <- data.frame(rep(sampleid[i],nSeg),seg.chrom,seg.arm,pos.start,pos.stop,seg.npos,seg.mean,stringsAsFactors=FALSE)
			colnames(seg) <- seg.names
			segments.c <- rbind(segments.c,seg)

      if(yest){
        #Rounding:
 			  yhat <- round(yhat,digits=digits)
        pcf.est.c <- cbind(pcf.est.c,yhat)
      }
    }#endfor


    #Should results be written to files or returned to user:
    if(save.res){
      if(c==1){
        #open connection for writing to file
        w1 <- file(file.names[1],"w")
        w2 <- file(file.names[2],"w")
      }
      #Write estimated PCF-values file for this arm:
      write.table(data.frame(chrom[probe.c],pos.c,pcf.est.c,stringsAsFactors=FALSE), file = w1,col.names=if(c==1) pcf.names else FALSE,row.names=FALSE,quote=FALSE,sep="\t")

      #Write segments to file for this arm
      write.table(segments.c,file=w2,col.names=if(c==1) seg.names else FALSE,row.names=FALSE,quote=FALSE,sep="\t")
    }


    #Append to results for other arms:
    segments <- rbind(segments,segments.c)
    if(return.est){
      pcf.est <- rbind(pcf.est,pcf.est.c)
    }

    if(verbose){
      cat(paste("pcf finished for chromosome arm ",this.chrom,this.arm,sep=""),"\n")
    }
 	}#endfor

	if(isfile.data){
    #Close connection
    close(f)
  }
  if(!is.null(Y)){
    if(isfile.Y){
      #Close connection
      close(f.y)
    }
  }


	if(save.res){
    close(w1)
    close(w2)
    cat(paste("pcf-estimates were saved in file",file.names[1]),sep="\n")
    cat(paste("segments were saved in file",file.names[2]),sep="\n")
	}
	if(return.est){
    pcf.est <- data.frame(chrom,position,pcf.est,stringsAsFactors=FALSE)
	  colnames(pcf.est) <- pcf.names
		return(list(estimates=pcf.est,segments=segments))
	}else{
	 return(segments)
	}

}#endfunction





### EXACT PCF-ALGORITHM
exactPcf <- function(y, kmin=5, gamma, yest) {
  ## Implementation of exact PCF by Potts-filtering
	## x: input array of (log2) copy numbers
	## kmin: Mininal length of plateaus
	## gamma: penalty for each discontinuity
	## yest: logical, should estimates for each pos be returned

  	N <- length(y)
  	yhat <- rep(0,N);
  	if (N < 2*kmin) {
	     if (yest) {
 		     return(list(Lengde = N, sta = 1, mean = mean(y), nIntervals=1, yhat=rep(mean(y),N)))
       } else {
 		     return(list(Lengde = N, sta = 1, mean = mean(y), nIntervals=1))
 	     }
  	}
  	initSum <- sum(y[seq_len(kmin)])
  	initAve <- initSum/kmin;
  	bestCost <- rep(0,N)
  	bestCost[kmin] <- (-initSum*initAve)
  	bestSplit <- rep(0,N)
  	bestAver <- rep(0,N)
  	bestAver[kmin] <- initAve
  	Sum <- rep(0,N)
  	Aver <- rep(0,N)
  	Cost <- rep(0,N)
  	kminP1=kmin+1
  	for (k in (kminP1):(2*kmin-1)) {
    		Sum[kminP1:k]<-Sum[kminP1:k]+y[k]
    		Aver[kminP1:k] <- Sum[kminP1:k]/((k-kmin):1)
    		bestAver[k] <- (initSum+Sum[kminP1])/k
    		bestCost[k] <- (-k*bestAver[k]^2)
  	}
  	for (n in (2*kmin):N) {
   		yn <- y[n]
   		yn2 <- yn^2
   		Sum[kminP1:n] <- Sum[kminP1:n]+yn
   		Aver[kminP1:n] <- Sum[kminP1:n]/((n-kmin):1)
   		nMkminP1=n-kmin+1
   		Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)]-Sum[kminP1:nMkminP1]*Aver[kminP1:nMkminP1]+gamma
   		Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
   		cost <- Cost[Pos]
   		aver <- Aver[Pos]
   		totAver <- (Sum[kminP1]+initSum)/n
   		totCost <-  (-n*totAver*totAver)
   		if (totCost < cost) {
       		Pos <- 1
          cost <- totCost
          aver <- totAver
   		}
   		bestCost[n] <- cost
   		bestAver[n] <- aver
   		bestSplit[n] <- Pos-1
 	}
 	n <- N
	antInt <- 0
	if(yest){
 		while (n > 0) {
   			yhat[(bestSplit[n]+1):n] <- bestAver[n]
   			n <- bestSplit[n]
			antInt <- antInt+1
 		}
 	} else {
 		while (n > 0) {
	   		n <- bestSplit[n]
			antInt <- antInt+1
 		}
 	}
	n <- N
	lengde <- rep(0,antInt)
	start <- rep(0,antInt)
	verdi <- rep(0,antInt)
	oldSplit  <- n
	antall <- antInt
	while (n > 0) {
    start[antall] <- bestSplit[n]+1
		lengde[antall] <- oldSplit-bestSplit[n]
		verdi[antall] <- bestAver[n]
		n <- bestSplit[n]
		oldSplit <- n
		antall <- antall-1
 	}
	if (yest) {
 		return(list(Lengde = lengde, sta = start, mean = verdi, nIntervals=antInt, yhat=yhat))
 	} else {
 		return(list(Lengde = lengde, sta = start, mean = verdi, nIntervals=antInt))
 	}
}


# Choose fast version
selectFastPcf <- function(x,kmin,gamma,yest){
	xLength <- length(x)
	if (xLength< 1000) {
		result<-runFastPcf(x,kmin,gamma,0.15,0.15,yest)
	} else {
	if (xLength < 15000){
		result<-runFastPcf(x,kmin,gamma,0.12,0.05,yest)
	} else  {
		result<-runPcfSubset(x,kmin,gamma,0.12,0.05,yest)
		}
	}
	return(result)
}

# Start fast version 1, moderately long sequences
runFastPcf <- function(x,kmin,gamma,frac1,frac2,yest){
	antGen <- length(x)
	mark<-filterMarkS4(x,kmin,8,1,frac1,frac2,0.02,0.9)
	mark[antGen]=TRUE
	dense <- compact(x,mark)
	result<-PottsCompact(kmin,gamma,dense$Nr,dense$Sum,yest)
  return(result)
}

# Start fast version 2, very long sequences
runPcfSubset <- function(x,kmin,gamma,frac1,frac2,yest){
	SUBSIZE <- 5000
	antGen <- length(x)
	mark<-filterMarkS4(x,kmin,8,1,frac1,frac2,0.02,0.9)
	markInit<-c(mark[seq_len(SUBSIZE-1)],TRUE)#markInit<-c(mark[1:(SUBSIZE-1)],TRUE)
	compX<-compact(x[seq_len(SUBSIZE)],markInit)#compX<-compact(x[1:SUBSIZE],markInit)
	mark2 <- rep(FALSE,antGen)
	mark2[seq_len(SUBSIZE)] <- markWithPotts(kmin,gamma,compX$Nr,compX$Sum,SUBSIZE)#mark2[1:SUBSIZE] <- markWithPotts(kmin,gamma,compX$Nr,compX$Sum,SUBSIZE)
  mark2[4*SUBSIZE/5]<-TRUE
	start <- 4*SUBSIZE/5+1
	while(start + SUBSIZE < antGen){
		slutt<-start+SUBSIZE-1
		markSub<-c(mark2[seq_len(start-1)],mark[start:slutt])#markSub<-c(mark2[1:(start-1)],mark[start:slutt])
		markSub[slutt] <- TRUE
		compX<-compact(x[seq_len(slutt)],markSub)#compX<-compact(x[1:slutt],markSub)
		mark2[seq_len(slutt)] <- markWithPotts(kmin,gamma,compX$Nr,compX$Sum,slutt)#mark2[1:slutt] <- markWithPotts(kmin,gamma,compX$Nr,compX$Sum,slutt)
    start <- start+4*SUBSIZE/5
		mark2[start-1]<-TRUE
	}
	markSub<-c(mark2[seq_len(start-1)],mark[seq(start,antGen,by=1)])#markSub<-c(mark2[1:(start-1)],mark[start:antGen])
	compX<-compact(x,markSub)
	result <- PottsCompact(kmin,gamma,compX$Nr,compX$Sum,yest)

  return(result)
}


# Calculations for fast version 1
PottsCompact <- function(kmin, gamma, nr, res, yest) {

  ## Potts filtering on compact array;
	## kmin: minimal length of plateau
	## gamma: penalty for discontinuity
	## nr: number of values between breakpoints
	## res: sum of values between breakpoints
	## sq: sum of squares of values between breakpoints

	N <- length(nr)
	Ant <- rep(0,N)
	Sum <- rep(0,N)
	Cost <- rep(0,N)
	if (sum(nr) < 2*kmin){
    estim <- sum(res)/sum(nr)
    return(estim)
	}
	initAnt <- nr[1]
	initSum <- res[1]
	initAve <- initSum/initAnt
	bestCost <- rep(0,N)
	bestCost[1] <- (-initSum*initAve)
  bestSplit <- rep(0,N)
	k <- 2
	while(sum(nr[seq_len(k)]) < 2*kmin) {#while(sum(nr[1:k]) < 2*kmin) {
	Ant[seq(2,k,by=1)] <- Ant[seq(2,k,by=1)]+nr[k]#Ant[2:k] <- Ant[2:k]+nr[k]
	Sum[seq(2,k,by=1)]<-Sum[seq(2,k,by=1)]+res[k]#Sum[2:k]<-Sum[2:k]+res[k]
    bestCost[k] <- (-(initSum+Sum[2])^2/(initAnt+Ant[2]))
    k <- k+1
	}
	for (n in k:N) {
    Ant[2:n] <- Ant[2:n]+nr[n]
 		Sum[2:n] <- Sum[2:n]+res[n]
    limit <- n
	  while(limit > 2 & Ant[limit] < kmin) {limit <- limit-1}
    Cost[seq(2,limit,by=1)] <- bestCost[seq_len(limit)-1]-Sum[seq(2,limit,by=1)]^2/Ant[seq(2,limit,by=1)]#Cost[2:limit] <- bestCost[1:limit-1]-Sum[2:limit]^2/Ant[2:limit]
    Pos <- which.min(Cost[2:limit])+ 1
    cost <- Cost[Pos]+gamma
    totCost <- (-(Sum[2]+initSum)^2/(Ant[2]+initAnt))
    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
    }
    bestCost[n] <- cost
    bestSplit[n] <- Pos-1
  }
  if (yest) {
    yhat<-rep(0,N)
    res<-findEst(bestSplit,N,nr,res,TRUE)
  } else {
    res<-findEst(bestSplit,N,nr,res,FALSE)
  }
  return(res)
}


# Accumulate information between potential breakpoints
compact <- function(y,mark){
  ## accumulates numbers of observations, sums and
  ## sums of squares between potential breakpoints
  ## y:  array to be compacted
  ## mark:  logical array of potential breakpoints
    tell <- seq(seq_len(length(y)))	#tell <- seq(1:length(y))
	cCTell <- tell[mark]
	Ncomp <- length(cCTell)
	lowTell <- c(0,cCTell[seq_len(Ncomp-1)])#lowTell <- c(0,cCTell[1:(Ncomp-1)])
	ant <- cCTell-lowTell
	cy <- cumsum(y)
	cCcy <- cy[mark]
	lowcy <- c(0,cCcy[seq_len(Ncomp-1)])#lowcy <- c(0,cCcy[1:(Ncomp-1)])
	sum <- cCcy-lowcy
	return(list(Nr=ant,Sum=sum))
}

#Retrieve segment information
findEst <- function(bestSplit,N,Nr,Sum,yest){
	n <- N
	lengde <- rep(0,N)
	antInt <- 0
	while (n>0){
		antInt <- antInt+1
		lengde[antInt] <- n-bestSplit[n]
		n <- bestSplit[n]
	}
	lengde <- lengde[antInt:1]
	lengdeOrig <- rep(0,antInt)
	startOrig <- rep(1,antInt+1)
	verdi <- rep(0,antInt)
	start <- rep(1,antInt+1)
	for(i in seq_len(antInt)){
		start[i+1] <- start[i]+lengde[i]
		lengdeOrig[i] <- sum(Nr[start[i]:(start[i+1]-1)])
		startOrig[i+1] <- startOrig[i]+lengdeOrig[i]
		verdi[i] <- sum(Sum[start[i]:(start[i+1]-1)])/lengdeOrig[i]
	}

	if(yest){
		yhat <- rep(0,startOrig[antInt+1]-1)
		for (i in seq_len(antInt)){
			yhat[startOrig[i]:(startOrig[i+1]-1)]<-verdi[i]
		}
		startOrig <- startOrig[seq_len(antInt)]#startOrig <- startOrig[1:antInt]
		return(list(Lengde=lengdeOrig,sta=startOrig,mean=verdi,nIntervals=antInt,yhat=yhat))
	} else {
	    startOrig <- startOrig[seq_len(antInt)]#startOrig <- startOrig[1:antInt]
		return(list(Lengde=lengdeOrig,sta=startOrig,mean=verdi,nIntervals=antInt))
	}

}


#Help function for very long sequence version
markWithPotts <- function(kmin, gamma, nr, res, subsize) {

  ## Potts filtering on compact array;
	## kmin: minimal length of plateau
	## gamma: penalty for discontinuity
	## nr: number of values between breakpoints
	## res: sum of values between breakpoints
	## sq: sum of squares of values between breakpoints

  N <- length(nr)
 	Ant <- rep(0,N)
 	Sum <- rep(0,N)
 	Cost <- rep(0,N)
	markSub <- rep(FALSE,N)
	initAnt <- nr[1]
	initSum <- res[1]
	initAve <- initSum/initAnt
	bestCost <- rep(0,N)
	bestCost[1] <- (-initSum*initAve)
  bestSplit <- rep(0,N)
	k <- 2
	while(sum(nr[seq_len(k)]) < 2*kmin) {#while(sum(nr[1:k]) < 2*kmin) {
    Ant[2:k] <- Ant[2:k]+nr[k]
    Sum[2:k]<-Sum[2:k]+res[k]
    bestCost[k] <- (-(initSum+Sum[2])^2/(initAnt+Ant[2]))
		k <- k+1
 	}
 	for (n in k:N) {
    Ant[2:n] <- Ant[2:n]+nr[n]
		Sum[2:n] <- Sum[2:n]+res[n]
		limit <- n
		while(limit > 2 & Ant[limit] < kmin) {limit <- limit-1}
	Cost[seq(2,limit,by=1)] <- bestCost[seq_len(limit)-1]-Sum[seq(2,limit,by=1)]^2/Ant[seq(2,limit,by=1)]#Cost[2:limit] <- bestCost[1:limit-1]-Sum[2:limit]^2/Ant[2:limit]
    Pos <- which.min(Cost[2:limit])+ 1
	  cost <- Cost[Pos]+gamma
    totCost <- (-(Sum[2]+initSum)^2/(Ant[2]+initAnt))
    if (totCost < cost) {
	    Pos <- 1
 			cost <- totCost
 		}
 		bestCost[n] <- cost
 		bestSplit[n] <- Pos-1
	  markSub[Pos-1] <- TRUE
  }
  help<-findMarks(markSub,nr,subsize)
  return(help=help)
}


#Find breakpoints on original scale
findMarks <- function(markSub,Nr,subsize){
	## markSub: marks in compressed scale
	## NR: number of observations between potenstial breakpoints
	mark <- rep(FALSE,subsize)  ## marks in original scale
	if(sum(markSub)<1) {
    return(mark)
  } else {
		N <- length(markSub)
		ant <- seq_len(N)#ant <- seq(1:N)
		help <- ant[markSub]
		lengdeHelp <- length(help)
		help0 <- c(0,help[seq_len(lengdeHelp-1)])#help0 <- c(0,help[1:(lengdeHelp-1)])
		lengde <- help-help0
		start <- 1
		oldStart <- 1
		startOrig <- 1
		for(i in seq_len(lengdeHelp)){
			start <- start+lengde[i]
			lengdeOrig <- sum(Nr[oldStart:(start-1)])
			startOrig <- startOrig+lengdeOrig
			mark[startOrig-1]<-TRUE
			oldStart<-start
		}
		return(mark)
	}
}



## marks potential breakpoints, partially by a two 6*L and 6*L2 highpass
## filters (L>L2), then by a filter seaching for potential kmin long segments
filterMarkS4 <- function(x,kmin,L,L2,frac1,frac2,frac3,thres){
  lengdeArr <- length(x)
  xc <- cumsum(x)
	xc <- c(0,xc)
	ind11 <- seq_len(lengdeArr-6*L+1)
	ind12 <- ind11+L
	ind13 <- ind11+3*L
	ind14 <- ind11+5*L
	ind15 <- ind11+6*L
	cost1 <- abs(4*xc[ind13]-xc[ind11]-xc[ind12]-xc[ind14]-xc[ind15])
	cost1 <- c(rep(0,3*L-1),cost1,rep(0,3*L))
	##mark shortening in here
	in1 <- seq_len(lengdeArr-6)
	in2 <- in1+1
	in3 <- in1+2
	in4 <- in1+3
	in5 <- in1+4
	in6 <- in1+5
	in7 <- in1+6
	test <- pmax(cost1[in1],cost1[in2],cost1[in3],cost1[in4],cost1[in5],cost1[in6],cost1[in7])
	test <- c(rep(0,3),test,rep(0,3))
	cost1B <- cost1[cost1>=thres*test]
	frac1B <- min(0.8,frac1*length(cost1)/length(cost1B))
	limit <- quantile(cost1B,(1-frac1B),names=FALSE)
	mark <- (cost1>limit)&(cost1>0.9*test)


	ind21 <- seq_len(lengdeArr-6*L2+1)
	ind22 <- ind21+L2
	ind23 <- ind21+3*L2
	ind24 <- ind21+5*L2
	ind25 <- ind21+6*L2
	cost2 <- abs(4*xc[ind23]-xc[ind21]-xc[ind22]-xc[ind24]-xc[ind25])
	limit2 <- quantile(cost2,(1-frac2),names=FALSE)
	mark2 <- (cost2>limit2)
	mark2 <- c(rep(0,3*L2-1),mark2,rep(0,3*L2))
	if(3*L>kmin){
		mark[kmin:(3*L-1)] <- TRUE
		mark[(lengdeArr-3*L+1):(lengdeArr-kmin)] <- TRUE
	}else{
		mark[kmin] <- TRUE
		mark[lengdeArr-kmin] <- TRUE
	}

	if(kmin>1){
		ind1 <- seq_len(lengdeArr-3*kmin+1)
		ind2 <- ind1+3*kmin
		ind3 <- ind1+kmin
		ind4 <- ind1+2*kmin
		shortAb <- abs(3*(xc[ind4]-xc[ind3])-(xc[ind2]-xc[ind1]))
		in1 <- seq_len(length(shortAb)-6)
		in2 <- in1+1
		in3 <- in1+2
		in4 <- in1+3
		in5 <- in1+4
		in6 <- in1+5
		in7 <- in1+6
		test <- pmax(shortAb[in1],shortAb[in2],shortAb[in3],shortAb[in4],shortAb[in5],shortAb[in6],shortAb[in7])
		test <- c(rep(0,3),test,rep(0,3))
		cost1C <- shortAb[shortAb>=thres*test]
		frac1C <- min(0.8,frac3*length(shortAb)/length(cost1C))
		limit3 <- quantile(cost1C,(1-frac1C),names=FALSE)
		markH1 <- (shortAb>limit3)&(shortAb>thres*test)
		markH2 <- c(rep(FALSE,(kmin-1)),markH1,rep(FALSE,2*kmin))
		markH3 <- c(rep(FALSE,(2*kmin-1)),markH1,rep(FALSE,kmin))
		mark <- mark|mark2|markH2|markH3
	} else {
		mark <- mark|mark2
	}

	if(3*L>kmin){
		mark[seq_len(kmin-1)] <- FALSE
		mark[kmin:(3*L-1)] <- TRUE
		mark[(lengdeArr-3*L+1):(lengdeArr-kmin)] <- TRUE
		mark[(lengdeArr-kmin+1):(lengdeArr-1)] <- FALSE
		mark[lengdeArr] <- TRUE
	}else{
		mark[seq_len(kmin-1)] <- FALSE
		mark[(lengdeArr-kmin+1):(lengdeArr-1)] <- FALSE
		mark[lengdeArr] <- TRUE
		mark[kmin] <- TRUE
		mark[lengdeArr-kmin] <- TRUE
	}

	return(mark)
}
