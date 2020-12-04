\name{readVCF}
\alias{readVCF}
\title{
Read In VCF File
}
\description{
This function can be used to read in the variation file in VCF format.
}
\usage{
readVCF(file)
}
\arguments{
  \item{file}{
a text file of VCF format.
}
}
\details{
The argument file, of which the columns of '#CHROM', 'POS', 'ID', 'REF', 'ALT',
'QUAL', and 'FILTER' are obliged. Detailed information of the VCF format please
refer to the Variant Call Format (VCF) Version 4.2 Specification.

The function will filter the variants according to the 'FILTER' column. Only all
the variants pass the filters that were set. If the 'FILTER' field has string
other than 'PASS', the function will rise a warning, and ask the user to decide
whether to stop and filter before use or continue with the analysis.

As the users usually used this function in PC, considering the memory capacity
of PCs and the usual size of the VCF files, the function only accept the files
of no larger than 2GB.
}
\value{
The function return a matrix of the variants passed all the filters with the
following columns:
  \item{#CHROM }{chromosome: An identifier from the reference genome or an
  angle-bracketed ID String (“<ID>”) pointing to a contig in the assembly file
  (cf. the ##assembly line in the header). All entries for a specific CHROM
  should form a contiguous block within the VCF file. (String, no white-space
  permitted, Required).}
  \item{POS }{position: The reference position, with the 1st base having
  position 1. Positions are sorted numerically, in increasing order, within
  each reference sequence CHROM. It is permitted to have multiple records with
  the same POS. Telomeres are indicated by using positions 0 or N+1, where N is
  the length of the corresponding chromosome or contig. (Integer, Required)}
  \item{REF }{reference base(s): Each base must be one of A,C,G,T,N (case
  insensitive). Multiple bases are permitted. The value in the POS field refers
  to the position of the first base in the String. For simple insertions and
  deletions in which either the REF or one of the ALT alleles would otherwise
  be null/empty, the REF and ALT Strings must include the base before the event
  (which must be reflected in the POS field), unless the event occurs at
  position 1 on the contig in which case it must include the base after the
  event; this padding base is not required (although it is permitted) for
  e.g. complex substitutions or other events where all alleles have at least one
  base represented in their Strings. If any of the ALT alleles is a symbolic
  allele (an angle-bracketed ID String “<ID>”) then the padding base is required
  and POS denotes the coordinate of the base preceding the polymorphism. Tools
  processing VCF files are not required to preserve case in the allele Strings.
  (String, Required).}
  \item{ALT }{alternate base(s): Comma separated list of alternate non-reference
  alleles. These alleles do not have to be called in any of the samples. Options
  are base Strings made up of the bases A,C,G,T,N,*, (case insensitive) or an
  angle-bracketed ID String (“<ID>”) or a breakend replacement string as
  described in the section on breakends. The ‘*’ allele is reserved to indicate
  that the allele is missing due to a upstream deletion. If there are no
  alternative alleles, then the missing value should be used. Tools processing
  VCF files are not required to preserve case in the allele String, except for
  IDs, which are case sensitive. (String; no whitespace, commas, or angle-brackets
  are permitted in the ID String itself)}
  \item{FILTER }{filter status: PASS if this position has passed all filters,
  i.e., a call is made at this position. Otherwise, if the site has not passed
  all filters, a semicolon-separated list of codes for filters that fail.
  e.g. “q10;s50” might indicate that at this site the quality is below 10 and
  the number of samples with data is below 50 percent  of the total number of
  samples.  ‘0’ is reserved and should not be used as a filter String. If filters
  have not been applied, then this field should be set to the missing value.
  (String, no white-space or semi-colons permitted) If the filter information
  is not in detail here, the program will take in all the variants and print out
  a warning message.}
}
\references{
Info, Q. F. (2020). The Variant Call Format ( VCF ) Version 4 . 2 Specification The VCF specification, 1–28.
}
\author{
Xue Lin, Jian Li
}
\note{
It is usually adviable to filter according to the filters,and convert your file
into the VCF 4.0 format before you read in your variants for further analysis.
And the .vcf file suffix is strongly recommendated.
}
\seealso{
readMAF
}
\examples{
#Read in VCF file:
fpath <- system.file("extdata", "test_mut.vcf", package="kataegis")
dat<-readVCF(fpath)

#Or you can also import vcf file typing
#dat<-readVCF("path/to/your/file.vcf")
}
\keyword{ IO }% use one of  RShowDoc("KEYWORDS")
