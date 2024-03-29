% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_funclib_movAPA.r
\name{ext3UTRPACds}
\alias{ext3UTRPACds}
\title{Extend 3'UTR}
\usage{
ext3UTRPACds(pacds, ext3UTRlen, extFtr = "3UTR")
}
\arguments{
\item{pacds}{a PACdataset object, which should already be annotated by annotatePAC().}

\item{ext3UTRlen}{an integer denoting the extended length, like 500bp.}

\item{extFtr}{set the ftr column of PACs in extended 3UTR as extFtr.}
}
\value{
A PACdataset with @anno$ftr column changed and $toStop column added.
PACs that are located in intergenic region but within ext3UTRlen bp of the upstream gene end are considered as 3'UTR PACs.
For extended 3'UTR PACs, the columns of ftr_start, ftr_end, gene_start, gene_end, gene_stop_codon are set as NA.
Note: 3'UTR extension is defined based on the three_extend column, so even for non-coding genes, there may also be PACs in the extended 3'UTR.
}
\description{
ext3UTRPACds extends annotated 3'UTR by a give length for a PACdataset object.
}
\details{
PACs in extended 3UTR region (@anno$ftr=intergenic) are set as 3UTR, then the 3UTR length (toStop) is calculated for all 3'UTR PACs.
}
\examples{
data(PACds)
ext3UTRPACds(PACds, 500)
}
\concept{PACdataset functions
\code{\link{annotatePAC}} to annotate a PACdataset.}
