## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 6,
  fig.height = 5.5,
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## ----include=FALSE, eval=TRUE-------------------------------------------------
library(movAPA)
## source('/media/bmi/My Passport/scPACext_HC_288cells/movAPA-latest/movAPA/R/R_funclib_movAPA.r')

## ---- eval=-1-----------------------------------------------------------------
library(movAPA, warn.conflicts = FALSE, quietly=TRUE)
data(scPACds)
head(scPACds@counts[1:2,1:5]) 
head(scPACds@anno, n=2)
head(scPACds@colData, n=2) 
levels(scPACds@colData$celltype)

## ----eval=FALSE---------------------------------------------------------------
#  library(TxDb.Mmusculus.UCSC.mm10.ensGene)
#  txdbmm10 <- parseGenomeAnnotation(TxDb.Mmusculus.UCSC.mm10.ensGene)
#  table(txdbmm10$anno.need$type)
#  gff=txdbmm10
#  save(gff, file='txdbmm10.rda')

## ----include=FALSE, eval=TRUE-------------------------------------------------
load("/media/bmi/My Passport/scPACext_HC_288cells/movAPA/movAPA/txdbmm10.rda")
gff=txdbmm10

## -----------------------------------------------------------------------------
table(scPACds@anno$ftr)
scPACds=ext3UTRPACds(scPACds, ext3UTRlen=2000)
table(scPACds@anno$ftr)

## -----------------------------------------------------------------------------
scPACdsNorm2=normalizePACds(scPACds, method='TPM') 
head(colSums(scPACdsNorm2@counts))

## -----------------------------------------------------------------------------
scPACdsFlt=subsetPACds(scPACds, totPACtag=20, choosePA=NULL,
                       noIntergenic=TRUE, verbose=TRUE) 

## -----------------------------------------------------------------------------
scPACdsFlt=get3UTRAPAds(scPACdsFlt, sortPA=TRUE, choose2PA=NULL)
length(scPACdsFlt)

## -----------------------------------------------------------------------------
scPACdsCt=subsetPACds(scPACds, group='celltype', pool=TRUE)

## -----------------------------------------------------------------------------
scPACdsCtStat=movStat(scPACdsCt, minPAT=c(20, 40), ofilePrefix=NULL)

## -----------------------------------------------------------------------------
scPACdsCtStat$pat20

## ----fig.dim=c(6,4), eval=FALSE-----------------------------------------------
#  plotPACdsStat(scPACdsCtStat, pdfFile=NULL)

## -----------------------------------------------------------------------------
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
d=melt(d)

## ----warning=FALSE, message=FALSE, fig.dim=c(7,7)-----------------------------
sp <- ggplot(data=d, aes(x=cell, y=value, color=variable)) +
  geom_bar(stat='identity', width=1.1)
sp = sp + theme_bw() +  guides(color=FALSE) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())
sp + facet_wrap( ~ variable, ncol=2, scales="free")

## -----------------------------------------------------------------------------
DEqPAC=movDEPAC(scPACdsCt, method='chisq', group='celltype',
                minSumPAT=50, chisqPadjust=TRUE)

## ----message=FALSE------------------------------------------------------------
stat=movStat(object=DEqPAC, padjThd=0.05, valueThd=0.95)
stat$nsig
head(stat$ovp)
head(stat$siglist[[1]])

## -----------------------------------------------------------------------------
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
cbind(gp@anno$ftr, rowSums(gp@counts))

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
movViz(object=scPACdsCt, gene=gene, txdb=gff, group='celltype')

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
movViz(object=DEqPAC, gene=gene, txdb=NULL,
       PACds=scPACdsCt,
       collapseConds=FALSE, padjThd=0.01,
       showRatio=FALSE, showAllPA=FALSE)

## ----fig.dim=c(7, 7),message=FALSE--------------------------------------------
movViz(object=DEqPAC, gene=gene, txdb=NULL, PACds=scPACdsCt, collapseConds=T, 
       padjThd=0.01, showPV=TRUE, showRatio=TRUE, showAllPA=T)

## ----message=FALSE------------------------------------------------------------
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

## -----------------------------------------------------------------------------
swDEq=movAPAswitch(PACds=scPACds, group='celltype',
                   avgPACtag=0, avgGeneTag=0,
                   only3UTR=TRUE, mergeReps='pool',
                   aMovDEPACRes=DEqPAC, DEPAC.padjThd=0.05, nDEPAC=1, 
                   mindist=50, fisherThd=0.05, logFCThd=1, 
                   cross=FALSE, selectOne='fisherPV')

## ----message=FALSE------------------------------------------------------------
swDEqStat=movStat(object=swDEq, padjThd=0.01, valueThd=1, upThd=NULL, dnThd=NULL)
swDEqStat$nsig

## ----message=FALSE, fig.dim=c(6,4)--------------------------------------------
nsig=as.data.frame(cbind(swstat$nsig, swDEqStat$nsig))
colnames(nsig)=c('Fisher.test', 'Fisher.test+DE')
nsig$cellType=rownames(nsig)
nsig=melt(nsig, variable.name='Method', value.name="SwitchingEvents")
ggplot(data=nsig, aes(x=cellType, y=SwitchingEvents, fill=Method)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("Switching Gene#") + theme_bw()

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
ds2=subsetPACds(scPACdsCt, group='celltype', pool=TRUE)
ds2=get3UTRAPAds(ds2, sortPA=TRUE, choose2PA='PD')
gpi2=movAPAindex(ds2, method="GPI")
summary(gpi2)

## ----fig.dim=c(4,3), warning=FALSE, message=FALSE-----------------------------
plotCummPAindex(PAindex=gpi2, groupName='cell type', xlab='GPI') 

## ----fig.dim=c(4,3), warning=FALSE, message=FALSE-----------------------------
shan=movPAindex(scPACdsCt, method="shan")
plotCummPAindex(PAindex=shan[, -c(1:3)], 
                groupName='cell type', xlab='Shannon index') 

## -----------------------------------------------------------------------------
cellCounts=colSums(scPACds@counts)
cellCounts=as.data.frame(cbind(cell=colnames(scPACds@counts),
                               celltype=as.character(scPACds@colData$celltype),
                               count=cellCounts))
cellCounts$count=as.integer(cellCounts$count)
head(cellCounts)
shownCells=cellCounts %>% dplyr::group_by(celltype) %>%
  dplyr::filter(rank(dplyr::desc(count)) <=10)
shownCells=shownCells[order(shownCells$celltype), ]
unique(shownCells$celltype)
shownCells=shownCells$cell

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
gene='ENSMUSG00000019969'
movViz(scPACds[, shownCells], gene=gene, txdb=NULL,
       group='group', collapseConds=FALSE,
       simple=TRUE, showYaxis=FALSE, geneHeight=0.1, 
       trackOrder=shownCells, ylimits=c(NA,NA),
       showBox=TRUE, boxWidth=10, showRatio=F, title=NULL)

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
gene='ENSMUSG00000019969'

movVizSC(scPACds, gene, cellGroupName='celltype', txdb=NULL,
         cellGroupColors=NULL, showRatio = F)

## Plot the expression levels of cell types as ratios, and specify colors for cell types.
movVizSC(scPACds, gene, cellGroupName='celltype', txdb=NULL,
         cellGroupColors=c(SC="yellow", ES="red",RS="blue"), showRatio = T)

## -----------------------------------------------------------------------------
sessionInfo()

