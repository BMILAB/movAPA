% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_funclib_GFF.r
\name{parseGenomeAnnotation}
\alias{parseGenomeAnnotation}
\title{Parse genome annotation}
\usage{
parseGenomeAnnotation(anAnno)
}
\arguments{
\item{anAnno}{can be a list of anno.rna/anno.need, or .rda/.rdata/.gff3/.gtf file name, or TxDB object.}
}
\value{
a parsed genome annotation object, which is a list of three elements (anno.rna, anno.need, anno.frame) and can be used for annotatePAC().
}
\description{
parseGenomeAnnotation parses different types of genome annotations.
}
\details{
Due to the complex GFF3/GTF/TxDB structure of different genome annotation files from different species,
this function may not be always applicable for any given file. You may need to check mannually.
}
\examples{
## Way1: Based on an annotation file in gff3 format,
## You can dowonload annotation from Ensemble Plants
#Prepare the annotation
#wget -c ftp://ftp.ensemblgenomes.org/pub/plants/release-44/gff3/
#        arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.44.gff3.gz
#gunzip  Arabidopsis_thaliana.TAIR10.44.gff3.gz
gff.path <- "/path/Arabidopsis_thaliana.TAIR10.44.gff3"
anno <- parseGenomeAnnotation(anAnno=gff.path)

##way2: load from a .rda file (already processed file)
anno <- parseGenomeAnnotation("anno.rda")

##Way3: Based on a TxDb object generated from BioMart.
# Parse Arabidopsis Txdb
library(TxDb.Athaliana.BioMart.plantsmart28)
anno <- parseGenomeAnnotation(TxDb.Athaliana.BioMart.plantsmart28)
# Parse mm10 Txdb
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
anno <- parseGenomeAnnotation(TxDb.Mmusculus.UCSC.mm10.ensGene)
}
\seealso{
[annotatePAC()] to annotate a PACdataset.

Other genome annotation functions: 
\code{\link{parseGff}()}
}
\concept{genome annotation functions}
