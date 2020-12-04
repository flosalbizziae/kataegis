####################################################################
## Author: Xue Lin, Jian Li.
## Maintainer: Xue Lin <xue.lin@njmu.edu.cn>
## License: Artistic 2.0
## Part of the Kataegis package
## Reference: To be published
####################################################################

#Function to read in the VCF files to control the format:

##Input:
### file: the name of the vcf file

##Required by:
### None

##Requires:
### None

#Read in VCF file:
readVCF<-function(file){

  if(!file.exists(file)){
    cat(paste("\n\nThe file ", file, " does not exists in your directory, please make sure before your input again.\n\n",sep=''))
    stop()
  }

	if(grepl(pattern="*\\.vcf$",x=file,ignore.case=TRUE)){
		if(file.size(file) < 2147483648){ #file size should not exceed 2GB
			tmp<-as.character(readLines(file))
			title<-tmp[which(grepl("^\\#",tmp))]
			title=title[length(title)]
			tmp<-tmp[which(!(grepl("^\\#",tmp)))]
			title<-strsplit(title,"\\s+")
			if(length(title[[1]])==1){
				title<-strsplit(title,"\t")
			}
			keep<-c()
			j=0
			newheader=c()
			for(i in seq_len(length(title[[1]]))){
				if(title[[1]][i]=="#CHROM"|title[[1]][i]=="POS"|title[[1]][i]=="REF"|title[[1]][i]=="ALT"|title[[1]][i]=="FILTER"){
					j=j+1
					keep[j]=i
					newheader[j]=title[[1]][i]
					if(title[[1]][i]=="FILTER"){#only the variations passing all the filters will be kept.
						filter=i
					}
				}
			}
			if(j==5){
				x=matrix(nrow=length(tmp),ncol=length(keep))
				n=0
				warn1=0
				for(i in seq_len(length(tmp))){
					line=strsplit(tmp[i],"\\s+")[[1]]
					if(length(line)==1){
						line=strsplit(tmp[i],"\t")[[1]]
					}
					if(line[filter]!="PASS"){
						warn1=1
					}
					n=n+1
					for(m in seq_len(length(keep))){
						x[n,m]=line[keep[m]]
					}
				}
				colnames(x)<-newheader
				#x<-x[-which(apply(x,1,function(y) all(is.na(y)))),]
				if(warn1){
					cat("Warning! The filter information of your file is not in detail, please make sure all the variants has passed all the filters.\n")
				}
				return(x) #return a matrix
			}
			else{
				cat("\nError! The header of the data in your VCF file is not parsable, please make sure you have the following columns in your data:\n\n")
				cat("\n#CHROM\tThe chromosome\n\n")
				cat("\nPOS\tThe start position of the variant\n\n")
				cat("\nREF\tThe reference sequence\n\n")
				cat("\nALT\tThe variant sequence\n\n")
				cat("\nFILTER\tThe filtered for the variants, which would be PASS for example\n\n")
			}
		}
		else{
			cat("\nWarning! Your file size exceed 2GB, Please filter out low quality variants or delete the annotations first.\n\n")
		}
	}
	else{
		cat("\nError! Please make sure to input files in VCF format.\n\n")
	}
}

