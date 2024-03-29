% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_funclib_movAPA.r
\name{mergePACds_v0}
\alias{mergePACds_v0}
\title{Merge multiple PACdatasets for RNA-seq}
\usage{
mergePACds_v0(PACdsList, d = 24)
}
\arguments{
\item{PACdsList}{is a PACdataset, or a list of multiple PACdataset objects. The PACds@anno should have columns chr/strand/coord.
If there is no colData in PACds, then the sample label will be set as groupN.
If PACdsList is a PACdataset, then will treat it as PA and group nearby PAs into PACs.}

\item{d}{distance to group nearby PACds, default is 24 nt.}
}
\value{
A merged PACdataset. The counts slot stores counts of merged samples.
If sample names from different PACdataset objects are any duplicated, then the sample name will be added a suffix .i for each sample in PACds[[i]].
The colData slot stores the merge sample annotation from the first column of each @colData.
The anno slot contains these columns: chr, strand, coord, tottag, UPA_start, UPA_end, nPA, maxtag.
}
\description{
mergePACds_v0 groups nearby PACs from single/multiple PACdataset objects.
This function is the mergePACds in previous version of movAPA, which is not recommanded to use now because it is not very suitable for (large) single-cell data.
This function is particularly useful for grouping nearby cleavage sites into PACs.
It is also useful When you have multiple PA or PAC files, each file is from one sample.
Then you need to merge these PACds into one PACds for DE or other analyses.
But after grouping and/or merging, you may need call annotatePAC to annotate the merged PACs by a GFF annotation.
}
\examples{
## Group PA into PACs
data(PACds)
PACds@counts=rbind(PACds@counts, PACds@counts)
PACds@anno=rbind(PACds@anno, PACds@anno)
ds=mergePACds_v0(PACds, d=24)
## merge two PACds
ds1=makeExamplePACds()
ds2=makeExamplePACds()
ds=mergePACds_v0(list(ds1, ds2), d=24)
}
\concept{PACdataset functions
\code{\link{annotatePAC}} to annotate a PACdataset; [rbind()] to combine multiple PACdatasets of the same format.}
