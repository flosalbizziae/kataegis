\name{kata}
\alias{kata}
\title{
Call The Kataegis Regions Of Your Data
}
\description{
The function will segment the inter-mutational distances using the pieceswise
constant fitting (PCF), then the kataegis will be called according to the
freuency of the variants' occurrance in a certain length of genomics region.
}
\usage{
kata(data, kmin = 2, gamma = 25, assembly = "hg19", len = 1000, nmut = 6, verbose = TRUE)
}
\arguments{
  \item{data}{
a dataframe or a list of dataframe of the inter-mutational distances calculated
and output by the interMutDist.
}
  \item{kmin}{
minimal number of variants in each segment, default is 2.
}
  \item{gamma}{
penalty for each discontinuity in the curve, default is 25.
}
  \item{assembly}{
the reference assembly version for analysis, default is human hg19, supported
assembly incl. hg38, hg16-19, and mm7-10.
}
  \item{len}{
the average inter-mutational distance to define the kataegis region.
}
  \item{nmut}{
the number of the variants contained in segments to define the kataegis region.
}
  \item{verbose}{
a verbose mode while segmenting the variants and calling kataegis.
}
}
\details{
If default parameters are used, the kataegis region (hypermutational region) will
be defined as a putative region with more than 6 variants withn a average inter-
muational distance of 1000bp.

The kmin and gamma are parameters used for PCF, and the default 2 and 25 were
trained on the set of kataegis foci that had been manually identified, curated
and validated using orthogonal sequencing platforms. See 'References'.
}
\value{
A dataframe or a list of dataframe of the kataegis region defined with the
following columns:
  \item{chr }{chromosome of the kataegis region locate}
  \item{arm }{the arm of the chromosome in which the kataegis region locate}
  \item{start.pos }{the start position of the kataegis of the chromosome}
  \item{end.pos }{the end position of the kataegis of the chromosome}
  \item{nmut }{the number of variants of the segment}
  \item{dsit.mean }{the average inter-mutational distance of the segment}
}
\references{
Lundmil B. Alexandrov, Serena Nik-Zaina, David C. Wedge, et.al(2013). Signatures
of mutational processes in human cancer. Nature 500, 415–421 (2013);
doi:10.1038/nature12477.

Gro Nilsen, Knut Liestol and Ole Christian Lingjaerde (2013). copynumber:
Segmentation of single- and multi-track copy number data by penalized least
squares regression. R package version 1.24.0.

Nilsen, G. & Liestol, K. et al. (2012) Copynumber: Efficient algorithms for
single- and multi-track copy number segmentation. BMC Genomics 13(1):591.
}
\author{
Xue Lin, Jian Li
}
\note{
Please be notified that the default parameters for identifying a kataegis region
is more than 6 variants within 1Kb genomic region, which is putative and it is
usually advisable to use the default parameters, but you can also set the
parameters according to your specific genomes.
}
\seealso{
interMutDistance
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
#lsit_dat<-readMAF("path/to/your/file.maf", split=TRUE)

##Read in VCF file:
fpath <- system.file("extdata", "test_mut.vcf", package="kataegis")
dat<-readVCF(fpath)

#Or you can also import vcf file typing
#dat<-readVCF("path/to/your/file.vcf")

##calculate the inter-mutational distances
dist<-interMutDist(dat)
#list_dist<-interMutDist(list_dat)

##call kataegis region with default parameters
#seg<-kata(dist)
#seg<-kata(list_dist)
}
\keyword{ math }% use one of  RShowDoc("KEYWORDS")
