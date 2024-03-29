% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_funclib_movAPA.r
\name{movAPAindexDiff}
\alias{movAPAindexDiff}
\title{Get APA index difference between conditions}
\usage{
movAPAindexDiff(
  PACds,
  group,
  method,
  doReps = "avg",
  choose2PA = NULL,
  RUD.includeNon3UTR = FALSE,
  clearPAT = 0,
  nPermute = 0
)
}
\arguments{
\item{PACds}{a PACdataset.}

\item{group}{a sample group name to get all conditions under the group.}

\item{method}{can be WUL, RUD, SLR.
\itemize{
\item{"WUL"}: Weighted UTR length (Jia 2017; Ultisky) (3UTRAPA, choose2PA=NULL/PD/MOST)
\item{"RUD"}: Relative Usage of Distal Poly(A) Sites (Ji 2009) (3UTRAPA with/without non-3UTR PACs, choose2PA=NULL/PD/MOST)
\item{"SLR"}: short to long ratio (Begik 2017) (3UTRAPA, choose2PA=PD/most). Similar to RUD (relative expression difference) (W Li 2015; Arefeen 2018). RED is equivalent to divide the SLR index but just choose PD/most PA.
\item{"RS"}: Patrick et al. 2016, sum(abs_ratio/2) (3UTRAPA if choose2PA!=null; allAPA if =NULL)
}}

\item{doReps}{Can be `avg` or `pool` to merge replicates.}

\item{choose2PA}{can be NULL, PD (proximal, distal, farest), most (most abundant two PAs).}

\item{RUD.includeNon3UTR}{only used when method=RUD, to use all PA or just 3UTR PA.}

\item{clearPAT}{set PAC counts to 0 when it is < clearPAT.
First call movAPAindex(clearPAT), if method=RS then subsetPACds(clearPAT).
If doReps!=NULL, then first doReps, and then filter by clearPAT.}

\item{nPermute}{If >0, then will get pvalue by bootstraping.
nPermute>0 must be with choose2PA=PD or MOST, because pvalue can be obtained for only two PACs.}
}
\value{
A data frame with columns being condition pairs like [cond1.cond2], and values are log2(index1/index2).
}
\description{
movAPAindexDiff calculates APA index difference between conditions by different metrics.
}
\details{
This function is different from APAswitching/UTR_trend. It does not filter PACds by cond1/2 but for all conditions.
So there might be NaN/Inf for some condition pairs when some PACs are not expressed in two specific condtions but expressed in other conds.
}
\examples{
data(PACds)
## RED + only for PD
diff=movAPAindexDiff(PACds, group='len', method="SLR", choose2PA="PD", RUD.includeNon3UTR=FALSE)
## RS+only for PD
diff=movAPAindexDiff(PACds, group='len', method="RS", choose2PA="PD")
## Get pvalue
geneSLRsig=movAPAindexDiff(PACds, group='len', method='SLR',
                          choose2PA="PD", doReps='avg', clearPAT=0, nPermute=10)
}
\seealso{
Other comparison functions: 
\code{\link{movAPAindex}()},
\code{\link{movAPAswitch}()},
\code{\link{movDEGene}()},
\code{\link{movDEPAC}()},
\code{\link{movPAindex}()},
\code{\link{movUTRtrend}()},
\code{\link{plotCummPAindex}()}
}
\concept{comparison functions}
