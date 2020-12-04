####################################################################
## Author: Xue Lin, Jian Li.
## Maintainer: Xue Lin <xue.lin@njmu.edu.cn>
## License: Artistic 2.0
## Part of the Kataegis package
## Reference: To be published
####################################################################

#Function to read in the MAF files to control the format:

##Input:
### file: the name of the maf file

##Required by:
### None

##Requires:
### None

#To read in the GDC maf format file
readMAF<-function(file, split=FALSE){ #the split option is used to control split samples or not
  if(!file.exists(file)){
    cat(paste("\n\nThe file ", file, " does not exists in your directory, please make sure before your input again.\n\n",sep=''))
    stop()
  }

	if(grepl(pattern="*\\.maf$",x=file,ignore.case=TRUE)){
		if(file.size(file)<2147483648){ #file size should not exceed 2GB

			tmp<-as.character(readLines(file))
			samples<-tmp[which(grepl("^\\#",tmp))]
			sample<-strsplit(samples,"\\s+")
			sample<-sample[[1]][2:length(sample[[1]])]

			title<-tmp[which(grepl("Start_Position",tmp))]
			title<-strsplit(title,"\\s+")

			sample_index=match("Tumor_Sample_Barcode",title[[1]])

			keep<-c()
			j=0
			newheader=c()
			for(i in seq_len(length(title[[1]]))){
				if(title[[1]][i]=="Chromosome"|title[[1]][i]=="Start_Position"|title[[1]][i]=="End_Position"|title[[1]][i]=="Reference_Allele"|title[[1]][i]=="Tumor_Seq_Allele2"){
					j=j+1
					keep[j]=i
					newheader[j]=title[[1]][i]
				}
			}
			if(j==5){
				x=matrix(nrow=length(tmp)-2,ncol=j+1)
				n=0
				for(i in 3:length(tmp)){
					line=strsplit(tmp[i],"\t")[[1]]
					n=n+1
					for(m in seq_len(length(keep))){
						if(title[[1]][keep[m]]=="Chromosome"){
							x[n,1]=line[keep[m]]
						}
						if(title[[1]][keep[m]]=="Start_Position"){
							start=as.numeric(line[keep[m]])
							end=as.numeric(line[keep[m]+1])
							pos=start+(end-start)/2 #get the central of the mutation location
							x[n,2]=pos
						}
						if(title[[1]][keep[m]]=="Reference_Allele"){
							x[n,3]=line[keep[m]]
							x[n,4]=line[keep[m]+2]
							x[n,5]="."
						}
					}
					x[n,j+1]=line[sample_index]
				}
				colnames(x)<-c("#CHROM","POS","REF","ALT","FILTER","SAMPLE_ID")

				if(split){ #if split is on then the file will be splitted to each sample
					result<-list()
					for(i in seq_len(length(sample))){
						result[[(sample[i])]]=x[which(x[,j+1]==sample[i]),seq_len(j)]
					}
					return(result)
				}
				else{ #if split is off then the file will be treated as one whole sample
					return(x[,seq_len(j)])
				}
			}
			else{
				cat("\nError! The header of the data in your MAF file is not parsable, please make sure you have the following columns in your data:\n\n")
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
		cat("\nError! Please make sure to input files in MAF format.\n\n")
	}
}
