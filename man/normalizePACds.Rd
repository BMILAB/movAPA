% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_funclib_movAPA.r
\name{normalizePACds}
\alias{normalizePACds}
\title{Normalize a PACdataset}
\usage{
normalizePACds(pacds, method = "TPM")
}
\arguments{
\item{pacds}{a PACdataset.}

\item{method}{normalization method, can be TPM/CPM/DEseq/DESeq2/EdgeR.
CPM=TPM (after normalization, libsize of all columns are 1e6), DESeq=DESeq2.}
}
\value{
A normalized PACdataset with pacds@counts normalized.
}
\description{
normalizePACds returns a normalized PACds with @counts normalized.
}
\details{
This function is used for normalization with different methods.
}
\examples{
data(PACds)
PACdsNorm=normalizePACds(PACds, method='TPM')
}
\seealso{
Other PACdataset functions: 
\code{\link{PACdataset-class}},
\code{\link{PACds}},
\code{\link{annotateByPAS}()},
\code{\link{annotatePAC}()},
\code{\link{createPACdataset}()},
\code{\link{get3UTRAPAds}()},
\code{\link{get3UTRAPApd}()},
\code{\link{length}()},
\code{\link{makeExamplePACds}()},
\code{\link{mergePACds}()},
\code{\link{plotPACdsStat}()},
\code{\link{rbind}()},
\code{\link{readPACds}()},
\code{\link{removePACdsIP}()},
\code{\link{scPACds}},
\code{\link{subscript_operator}},
\code{\link{summary}()},
\code{\link{writePACds}()}
}
\concept{PACdataset functions}
