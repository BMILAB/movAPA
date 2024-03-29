% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_funclib_movAPA.r
\name{plotSeqLogo}
\alias{plotSeqLogo}
\title{Plot seqlogo for poly(A) signals.}
\usage{
plotSeqLogo(pas)
}
\arguments{
\item{pas}{a vector of k-grams, normally from annotateByPAS().}
}
\value{
A plot of seqlogo.
}
\description{
plotSeqLogo plots seqlogo for poly(A) signals by utilizing the motifStack package.
}
\examples{
data(PACds)
library("BSgenome.Oryza.ENSEMBL.IRGSP1")
bsgenome <- BSgenome.Oryza.ENSEMBL.IRGSP1
PACdsPAS=annotateByPAS(PACds, bsgenome, grams='V1',
                       from=-50, to=-1, label=NULL)
pas=PACdsPAS@anno$V1_gram[!is.na(PACdsPAS@anno$V1_gram)]
plotSeqLogo(pas)
}
\seealso{
Other APA signal functions: 
\code{\link{annotateByPAS}()},
\code{\link{faFromPACds}()},
\code{\link{getVarGrams}()},
\code{\link{kcount}()},
\code{\link{plotATCGforFAfile}()}
}
\concept{APA signal functions}
