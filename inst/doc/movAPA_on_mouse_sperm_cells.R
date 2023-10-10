## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 6,
  fig.height = 5.5,
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(movAPA, warn.conflicts = FALSE, quietly=TRUE)
data(scPACds)
summary(scPACds)
head(scPACds@counts[1:2,1:5]) 
head(scPACds@anno, n=2)
head(scPACds@colData, n=2) 
levels(scPACds@colData$celltype)

## ----parseTxdb, eval=FALSE----------------------------------------------------
#  library(TxDb.Mmusculus.UCSC.mm10.ensGene)
#  txdbmm10 <- TxDb.Mmusculus.UCSC.mm10.ensGene
#  
#  scPACds=createPACdataset(scPACds@counts, scPACds@anno, scPACds@colData, forceSparse = TRUE)
#  scPACds=annotatePAC(scPACds, txdbmm10)

## ----testExt3UTR--------------------------------------------------------------
testExt3UTR(scPACds, seq(1000, 10000, by=1000))

## ----ext3UTR------------------------------------------------------------------
table(scPACds@anno$ftr)
scPACds=ext3UTRPACds(scPACds, ext3UTRlen=2000)
table(scPACds@anno$ftr)

## ----normalize----------------------------------------------------------------
scPACdsNorm2=normalizePACds(scPACds, method='TPM') 
head(Matrix::colSums(scPACdsNorm2@counts))

## ----subsetPACds1-------------------------------------------------------------
scPACdsFlt=subsetPACds(scPACds, totPACtag=20, choosePA=NULL,
                       noIntergenic=TRUE, verbose=TRUE) 
summary(scPACdsFlt)

## ----subsetPACds2-------------------------------------------------------------
scPACdsFlt=get3UTRAPAds(scPACdsFlt, sortPA=TRUE, choose2PA=NULL)
summary(scPACdsFlt)

## ----pool---------------------------------------------------------------------
scPACdsCt=subsetPACds(scPACds, group='celltype', pool=TRUE)

## ----movStat------------------------------------------------------------------
scPACdsCtStat=movStat(scPACdsCt, minPAT=c(20, 40), ofilePrefix=NULL)

## -----------------------------------------------------------------------------
scPACdsCtStat$pat20

## ----plotPACdsStat, fig.dim=c(6,4)--------------------------------------------
plotPACdsStat(scPACdsCtStat, pdfFile=NULL)

## ----movStat_sc---------------------------------------------------------------
scPACdsStat=movStat(scPACds, minPAT=c(1,5), ofilePrefix='scPACds.stat')

## -----------------------------------------------------------------------------
scPACdsStat$pat1['total',]

## -----------------------------------------------------------------------------
summary(scPACdsStat$pat1$nPAC[1:(nrow(scPACdsStat$pat1)-1)])

## -----------------------------------------------------------------------------
summary(scPACdsStat$pat1$nPAT[1:(nrow(scPACdsStat$pat1)-1)])

## ----message=FALSE------------------------------------------------------------
d=scPACdsStat$pat1[, c('nPAC','nPAT','nGene','nAPAgene','APAextent')]
d$cell=rownames(d)
d=d[1:(nrow(d)-1), ]
d=reshape2::melt(d)

## ----plot_sc, warning=FALSE, message=FALSE, fig.dim=c(7,7)--------------------
library(ggplot2, quietly = TRUE)
sp <- ggplot(data=d, aes(x=cell, y=value, color=variable)) +
  geom_bar(stat='identity', width=1.1)
sp = sp + theme_bw() +  guides(color=FALSE) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())
sp + facet_wrap( ~ variable, ncol=2, scales="free")

## ----movDEPAC-----------------------------------------------------------------
DEqPAC=movDEPAC(scPACdsCt, method='chisq', group='celltype',
                minSumPAT=50, chisqPadjust=TRUE)

## ----movStat_DE, message=FALSE------------------------------------------------
stat=movStat(object=DEqPAC, padjThd=0.05, valueThd=0.95)
stat$nsig
head(stat$ovp)
head(stat$siglist[[1]])

## ----movSelect_DE-------------------------------------------------------------
sel=movSelect(aMovRes=DEqPAC, condpair='SC.RS', 
              padjThd=0.05, valueThd=0.95, 
              out='full', PACds=scPACdsCt)
head(sel, n=2)

## ----eval=-3------------------------------------------------------------------
head(stat$tf01)
## Output stat results into files: "DEqPAC.plots.pdf" and 'DEqPAC.stat'.
outputHeatStat(heatStats=stat, ostatfile='DEqPAC.stat', plotPre='DEqPAC')

## -----------------------------------------------------------------------------
gene='ENSMUSG00000019969'
gp=scPACds[scPACds@anno$gene==gene, ]
cbind(gp@anno$ftr, Matrix::rowSums(gp@counts))

## ----gff----------------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdbmm10 <- TxDb.Mmusculus.UCSC.mm10.ensGene
gff=parseGenomeAnnotation(txdbmm10)

## ----movViz, fig.dim=c(7, 5), message=FALSE-----------------------------------
movViz(object=scPACdsCt, gene=gene, txdb=gff, group='celltype')

## ----movViz_nogff, fig.dim=c(7, 5), message=FALSE-----------------------------
movViz(object=DEqPAC, gene=gene, txdb=NULL,
       PACds=scPACdsCt,
       collapseConds=FALSE, padjThd=0.01,
       showRatio=FALSE, showAllPA=FALSE)

## ----fig.dim=c(7, 7), message=FALSE-------------------------------------------
movViz(object=DEqPAC, gene=gene, txdb=NULL, PACds=scPACdsCt, collapseConds=T, 
       padjThd=0.01, showPV=TRUE, showRatio=TRUE, showAllPA=T)

## ----movAPAswitch, message=FALSE----------------------------------------------
sw=movAPAswitch(PACds=scPACds, group='celltype',
                avgPACtag=0, avgGeneTag=0,
                only3UTR=TRUE, mergeReps='pool',
                aMovDEPACRes=NULL, DEPAC.padjThd=NULL, nDEPAC=0,
                mindist=0, fisherThd=0.05, logFCThd=0, cross=FALSE,
                selectOne='fisherPV')

head(sw@fullList$SC.RS, n=2)
## Filter results with padj<0.05.
swstat=movStat(object=sw, padjThd=0.05, valueThd=0)
## Number of significant switching genes between cell types.
swstat$nsig 
head(swstat$siglist$SC.RS)

## ----movAPAswitch2------------------------------------------------------------
swDEq=movAPAswitch(PACds=scPACds, group='celltype',
                   avgPACtag=0, avgGeneTag=0,
                   only3UTR=TRUE, mergeReps='pool',
                   aMovDEPACRes=DEqPAC, DEPAC.padjThd=0.05, nDEPAC=1, 
                   mindist=50, fisherThd=0.05, logFCThd=1, 
                   cross=FALSE, selectOne='fisherPV')

## ----message=FALSE------------------------------------------------------------
swDEqStat=movStat(object=swDEq, padjThd=0.01, valueThd=1, upThd=NULL, dnThd=NULL)
swDEqStat$nsig

## ----movAPAswitch_compare, message=FALSE, fig.dim=c(6,4)----------------------
nsig=as.data.frame(cbind(swstat$nsig, swDEqStat$nsig))
colnames(nsig)=c('Fisher.test', 'Fisher.test+DE')
nsig$cellType=rownames(nsig)
nsig=reshape2::melt(nsig, variable.name='Method', value.name="SwitchingEvents")
ggplot(data=nsig, aes(x=cellType, y=SwitchingEvents, fill=Method)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("Switching Gene#") + theme_bw()

## ----movRes2heatmapResults----------------------------------------------------
heat=movRes2heatmapResults(swDEq)

## -----------------------------------------------------------------------------
heat=subsetHeatmap(heat, padjThd=0.001, valueThd=2)
nrow(heat@value)

## -----------------------------------------------------------------------------
heat@value=heat@value[rowSums(is.na(heat@value))==0, ]
heat@value=heat@value[order(rowMeans(abs(heat@value)), decreasing =T ), ]
nrow(heat@value)

## ----fig.dim=c(6,4), message=FALSE, warning=FALSE-----------------------------
plotHeatmap(heat@value[1:20, ],  show_rownames=TRUE, plotPre=NULL)
heat@value['ENSMUSG00000021281',]

## -----------------------------------------------------------------------------
swDEq@fullList$SC.RS[swDEq@fullList$SC.RS$gene=='ENSMUSG00000021281',]

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
gene='ENSMUSG00000021281'
movViz(object=swDEq, gene=gene, txdb=gff, PACds=scPACdsCt)

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
movViz(object=swDEq, gene=gene, txdb=NULL, PACds=scPACdsCt)

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
movViz(object=swDEq, gene=gene, txdb=NULL, PACds=scPACdsCt, collapseConds=TRUE,
       conds=NULL, highlightConds=NULL, showRatio=TRUE, linkPAs=TRUE, 
       padjThd=0.01, showAllPA=FALSE, showPV=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  ds=get3UTRAPAds(scPACds, sortPA=TRUE, choose2PA='PD')
#  gpi=movAPAindex(ds, method="GPI")
#  head(gpi[1:5, 1:5])
#  gpi=gpi[rowSums(is.na(gpi))==0, ]

## ----movAPAindex_GPI----------------------------------------------------------
ds2=subsetPACds(scPACdsCt, group='celltype', pool=TRUE)
ds2=get3UTRAPAds(ds2, sortPA=TRUE, choose2PA='PD')
gpi2=movAPAindex(ds2, method="GPI")
summary(gpi2)

## ----fig.dim=c(4,3), warning=FALSE, message=FALSE-----------------------------
plotCummPAindex(PAindex=gpi2, groupName='cell type', xlab='GPI') 

## ----movAPAindex_shan, fig.dim=c(4,3), warning=FALSE, message=FALSE-----------
shan=movPAindex(scPACdsCt, method="shan")
plotCummPAindex(PAindex=shan[, -c(1:3)], 
                groupName='cell type', xlab='Shannon index') 

## ----movVizSC, fig.dim=c(7, 5), message=FALSE---------------------------------
gene='ENSMUSG00000019969'

movVizSC(scPACds, gene, cellGroupName='celltype', txdb=NULL,
         cellGroupColors=NULL, showRatio = F)

## Plot the expression levels of cell types as ratios, and specify colors for cell types.
movVizSC(scPACds, gene, cellGroupName='celltype', txdb=NULL,
         cellGroupColors=c(SC="yellow", ES="red",RS="blue"), showRatio = T)

## -----------------------------------------------------------------------------
sessionInfo()

