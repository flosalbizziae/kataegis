\name{readMAF}
\alias{readMAF}
\title{
Read In MAF File
}
\description{
This function can be used to read in the variation file in MAF format in a merge mode or a seperated mode.
}
\usage{
readMAF(file, split = FALSE)
}
\arguments{
  \item{file}{
The argument file, of which the columns of 'Chromosome', 'Start_Position',
'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', are obliged. Detailed
information of the MAF format, please refer to the Mutation Annotation Format
(MAF) GDC Version 1.0.0 Specification.

As the users usually used this function in PC, considering the memory capacity
of PCs and the usual size of the VCF files, the function only accept the files
of no larger than 2GB.
}
  \item{split}{
The boolean value to define whether the MAF file is splitted according to the
samples ID. If the argument is TRUE, the MAF will be splitted into a list of
several matrixes of the variants. Default is FALSE. See 'Details'.
}
}
\details{
The argument split is TRUE, the MAF file will be splitted to each sample. In
other word, after the readMAF function with the splitted argument set TRUE, each
sample will has its own matrix returned, with the columns exactly like the
readVCF does, and all the matrixes will be gathered in a list, in which each
matrix can be accessed by its corresponding sample ID.
}
\value{
If the split argument is not set, a matrix will be returned, and the structures
is exactly like the matrix returned by readVCF. Otherwise a list of matrix will
be returned, and the structure will be as followings:
List
|_sample_ID_1
| |_matrix_1
|_sample_ID_2
| |_matrix_2
|_sample_ID_3
  |_matrix_3
}
\references{
https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
}
\author{
Xue Lin, Jian Li
}
\note{
It is usually adviable to filter according to the filters,and convert your file
into the GDC MAF v1.0.0 format before you read in your variants for further analysis.
And the .maf file suffix is strongly recommendated.
}
\seealso{
readMAF
}
\examples{
#Read in MAF file and merge all the samples in one file:
fpath <- system.file("extdata", "test_mut.maf", package="kataegis")
dat<-readMAF(fpath)

#Read in MAF file and split the maf according the samples' ID:
list_dat<-readMAF(fpath, split=TRUE)

#Or you can also import maf file by typing
#dat<-readMAF("path/to/your/file.maf")
#list_dat<-readMAF("path/to/your/file.maf", split=TRUE)
}
\keyword{ IO }
