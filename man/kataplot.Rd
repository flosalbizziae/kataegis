\name{kataplot}
\alias{kataplot}
\title{
Kataegis Region Visualization
}
\description{
This function will use the read-in data, the inter-mutational distance data, and
the kataegis segments to visualize the nucleotide spectra of the kataegis foci,
the necleotide conversion type, and the relation of the distribution of the
variants on the genome and the inter-mutational distances.
}
\usage{
kataplot(x, y, z, n = 20, fbar = TRUE, cbar = TRUE, rain = TRUE, colr = 6, fbw = 800,
fbh = 400, rw = 1000, rh = 500, assembly = "hg19", name = "all", type = "png")
}
\arguments{
  \item{x}{
the matrixm or a list of matrixes of the read in data from readVCF or readMAF
}
  \item{y}{
the dataframe or a list of dataframes of the kataegis data
}
  \item{z}{
the dataframe or a list of dataframes output from the interMutDist() function
}
  \item{n}{
the number of the division of the flanking region totally, default is 20, which
means that the upstream is devided into 10 and the downstream is devided into 10
}
  \item{fbar}{
the flanking region nucleotide spectra barplot
}
  \item{cbar}{
the core region mutation types barplot
}
  \item{rain}{
the rain-fall plot of all the raw data, and the poisson will also be tested
}
  \item{colr}{
the color set of the barplot and rainfall plot, default is the color palettes
from grDevices, rainbow(6),which is "#FF0000" "#FFFF00" "#00FF00" "#00FFFF"
"#0000FF" "#FF00FF". You can also use any R colors here. Theis value should be
equal or larger than 6.
}
  \item{fbw}{
the width of the flanking region nucleotides types barplot
}
  \item{fbh}{
the height of the flanking region nucleotides types barplot
}
  \item{rw}{
the width of the rainfall plot
}
  \item{rh}{
the height of the rainfall plot
}
  \item{assembly}{
the reference assembly version for analysis, default is human hg19, supported
assembly incl. hg38, hg16-19, and mm7-10.
}
  \item{name}{
the sample name for the plot output, default is all.
}
  \item{type}{
the plot output format, default is png, other are jpeg, and tiff
}
}
\value{
The details of the read in data please refere the manual page of readVCF(), readMAF(),
interMutDist() and kata().
}
\references{
Lundmil B. Alexandrov, Serena Nik-Zaina, David C. Wedge, et.al(2013). Signatures
of mutational processes in human cancer. Nature 500, 415–421 (2013);
doi:10.1038/nature12477.
}
\author{
Xue Lin, Jian Li
}
\note{
The output format and width and height of graphs could be set according to the
graph type, and the size of your data. The rainfall plot will use more memory,
in this case you could consider for good filtering at the very beginning of the
analysis.
}
\seealso{
kata, readVCF, readMAF
}
\examples{
#library(kataegis)

##Read in MAF file and merge all the samples in one file:
fpath <- system.file("extdata", "test_mut.maf", package="kataegis")
dat<-readMAF(fpath)

#Read in MAF file and split the maf according the samples' ID:
list_dat<-readMAF(fpath, split=TRUE)

#Or you can also import maf file typing
#dat<-readMAF("path/to/your/file.maf")
#list_dat<-readMAF("path/to/your/file.maf", split=TRUE)

##Read in VCF file:
fpath <- system.file("extdata", "test_mut.vcf", package="kataegis")
dat<-readVCF(fpath)

#Or you can also import vcf file typing
#dat<-readVCF("path/to/your/file.vcf")

##Calculate the inter-mutational distances
dist<-interMutDist(dat)
#list_dist<-interMutDist(list_dat)

#Call kataegis region with default parameters
#seg<-kata(dist)
#or
#list_seg<-kata(list_dist)

##Plot the result of the analysis
#kataplot(dat, seg, dist)
#or
#kataplot(lsit_dat, list_seg, list_dist)
}
\keyword{ graph }% use one of  RShowDoc("KEYWORDS")
