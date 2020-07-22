\name{BSgenome.Oryza.ENSEMBL.IRGSP1}
\docType{package}

\alias{BSgenome.Oryza.ENSEMBL.IRGSP1-package}
\alias{BSgenome.Oryza.ENSEMBL.IRGSP1}
\alias{Oryza}

\title{Full genome sequences for Oryza sataiva japonica (ENSEMBL version IRGSP1.0)}

\description{
  Full genome sequences for Oryza sataiva japonica as provided by Ensembl (IRGSP1.0, November, 2018).
}

\details{
  
}

\note{
  This BSgenome data package was made from the following source data files:
  \preformatted{
Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz from ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/oryza_sativa/dna/
  }

  See \code{?\link[BSgenome]{BSgenomeForge}} and the BSgenomeForge
  vignette (\code{vignette("BSgenomeForge")}) in the \pkg{BSgenome}
  software package for how to make a BSgenome data package.
}

\author{The Bioconductor Dev Team}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects and the
          \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.
    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.
    \item The BSgenomeForge vignette (\code{vignette("BSgenomeForge")})
          in the \pkg{BSgenome} software package for how to make a BSgenome
          data package.
  }
}

\examples{
BSgenome.Oryza.ENSEMBL.IRGSP1
genome <- BSgenome.Oryza.ENSEMBL.IRGSP1
head(seqlengths(genome))


## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------
## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
