% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/movVizSc.R
\name{movVizSC}
\alias{movVizSC}
\title{Visualize PACs of a gene in single cells}
\usage{
movVizSC(
  scPACds,
  gene,
  cellGroupName,
  cellGroupColors = NULL,
  txdb = NULL,
  PAwidth = 50,
  showRatio = FALSE
)
}
\arguments{
\item{scPACds}{a PACdataset storing single cell PACs. The PACds@anno should have columns chr/strand/coord.
If there is no colData in PACds, then the sample label will be set as groupN.}

\item{gene}{a gene to plot.}

\item{cellGroupName}{group name of cell type in column names of scPACds@colData.}

\item{cellGroupColors}{a named vector denoting colors for cell groups, e.g., c('SC="yellow",ES="red",RS="blue").
its names must be the same as the labels of cell groups.}

\item{txdb}{txdb annotation, a character of file name for .gff3 or .rda, or a list. See parseGenomeAnnotation().}

\item{PAwidth}{expand the width of a PAC by PAwidth for the plot.}

\item{showRatio}{default is FALSE. If TRUE then show PAC ratio instead of PAC count.}
}
\value{
a plot with the top panel showing expression levels of PACs in single cells, the middle panel showing expression levels or ratios of PACs among cell types, the bottome panel showing the gene model and PAC.
}
\description{
movVizSC visualizes PACs of a gene in single cells.
}
\details{
This function is useful for visualizing distributions of expression levels of PACs among single cells and cell types. (Note: the code of this function was adapted from the Millefy R package [Ozaki, H., Hayashi, T., Umeda, M., et al. Millefy: visualizing cell-to-cell heterogeneity in read coverage of single-cell RNA sequencing datasets. BMC Genomics 2020;21(1):177. ])
}
\examples{
\dontrun{
data("scPACds")
load('txdbmm10.rda')
#gene <- unique(scPACds@anno$gene)[2]
gene='ENSMUSG00000019969'
movVizSC(scPACds, gene, cellGroupName="celltype", cellGroupColors=NULL, txdb=txdbmm10)
movVizSC(scPACds, gene, cellGroupName="celltype",
        cellGroupColors = c(SC="yellow",ES="red",RS="blue"), txdb=txdbmm10, showRatio = T)
movVizSC(scPACds, gene, cellGroupName="celltype", PAwidth=100, txdb=NULL, showRatio = T)
}
}
\seealso{
[movViz()] to visualize a PACdataset.
}
