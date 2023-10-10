## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 6,
  fig.height = 5.5,
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## ----load_PACds---------------------------------------------------------------
library(movAPA, warn.conflicts = FALSE, quietly=TRUE)
data("PACds")
PACds
summary(PACds)
# Transform the older version of PACdataset to newer version; the counts slot was converted from data.frame to matrix or dgCMatrix.
# PACds@counts=asAnyMatrix(PACds@counts)

## ----message=FALSE, eval=FALSE------------------------------------------------
#  devtools::load_all("/media/bmi/My Passport/scPACext_HC_288cells/movAPA/movAPA/BSgenome.Oryza.ENSEMBL.IRGSP1")

## ----bsgenome, message=FALSE--------------------------------------------------
library("BSgenome.Oryza.ENSEMBL.IRGSP1", quietly = TRUE)
bsgenome <- BSgenome.Oryza.ENSEMBL.IRGSP1

## ----parseGff, eval=FALSE-----------------------------------------------------
#  gffFile="Oryza_sativa.IRGSP-1.0.42.gff3"
#  gff=parseGff(gffFile)
#  save(gff, file='Oryza_sativa.IRGSP-1.0.42.gff.rda')

## -----------------------------------------------------------------------------
load('Oryza_sativa.IRGSP-1.0.42.gff.rda')

## ----removePACdsIP------------------------------------------------------------
PACdsIP=removePACdsIP(PACds, bsgenome, returnBoth=TRUE, 
                      up=-10, dn=10, conA=6, sepA=7)
length(PACdsIP$real)
length(PACdsIP$ip)

## ----results='hide'-----------------------------------------------------------
PACdsClust=mergePACds(PACds, d=100)

## -----------------------------------------------------------------------------
summary(PACds)
summary(PACdsClust)

## -----------------------------------------------------------------------------
## Constuct another demo PACdataset for merging.
PACds2=PACds
PACds2@anno$coord = PACds2@anno$coord + sample(-50:50, 1)

## You may also change the sample names and group names.
# rownames(PACds2@colData)=paste0(rownames(PACds2@colData),'v2')
# PACds2@colData$group=paste0(PACds2@colData$group,'v2')
# colnames(PACds2@counts)=paste0(colnames(PACds2@counts),'v2')
## Construct a list of PACds to be merged.
PACdsList=list(pac1=PACds, pac2=PACds2)

## ----mergePACds_multi---------------------------------------------------------
## Merge two PACdatasets, nearby PACs within 24bp of each other 
## will be merged into one PAC.
pp=mergePACds(PACdsList, d=24)
summary(pp)

## ----normalizePACds-----------------------------------------------------------
## Here normalization method TMM (or EdgeR) is used, 
## while you may also choose TPM or DESeq.
PACds=normalizePACds(PACds, method='TMM')

## Library sizes after normalization.
colSums(PACds@counts)

## ----eval=FALSE---------------------------------------------------------------
#  load('Oryza_sativa.IRGSP-1.0.42.gff.rda')

## ----annotatePAC--------------------------------------------------------------
PACds1=PACds
PACds1@anno[,c('gene','ftr','gene_type','ftr_start','ftr_end')]=NULL
PACds1=annotatePAC(PACds1, gff)

## ----eval=FALSE---------------------------------------------------------------
#  writePACds(PACds1, file='rice_pac_data.txt',
#             colDataFile = 'rice_pac_data.coldata.txt')

## ----testExt3UTR--------------------------------------------------------------
testExt3UTR(PACds1, seq(1000, 10000, by=1000))

## ----ext3UTR------------------------------------------------------------------
table(PACds1@anno$ftr)
PACds1=ext3UTRPACds(PACds1, ext3UTRlen=2000)
table(PACds1@anno$ftr)

## -----------------------------------------------------------------------------
PACds1=subsetPACds(PACds, group='group', pool=TRUE)
head(PACds1@counts)

## ----movStat------------------------------------------------------------------
pstats=movStat(PACds1, minPAT=c(1, 5, 10, 20, 50, 60), ofilePrefix=NULL)
names(pstats)
pstats$pat10

## ----eval=FALSE---------------------------------------------------------------
#  plotPACdsStat(pstats, pdfFile='PACds_stat.pdf', minPAT=c(5,10))

## ----eval=FALSE---------------------------------------------------------------
#  plotPACdsStat(pstats, pdfFile='PACds_stat_anther_embryo.pdf',
#                minPAT=c(5,10), conds=c('anther1','embryo1'))

## ----eval=FALSE---------------------------------------------------------------
#  plotPACdsStat(pstats, pdfFile='PACds_stat_total.pdf',
#                minPAT=c(5,10), conds=c('total'))

## ----eval=FALSE---------------------------------------------------------------
#  plotPACdsStat(pstats, pdfFile='PACds_stat_total_PAT10.pdf',
#                minPAT=c(10), conds=c('total'))

## ----fig.dim = c(6, 4), fig.align='left'--------------------------------------
plotPACdsStat(pstats, pdfFile=NULL, minPAT=c(5,10))

## -----------------------------------------------------------------------------
PACdsPAS=annotateByPAS(PACds, bsgenome, grams='AATAAA', 
                       from=-50, to=-1, label=NULL)
summary(PACdsPAS@anno$AATAAA_dist)

## -----------------------------------------------------------------------------
PACdsPAS=annotateByPAS(PACds, bsgenome, grams='V1', 
                       from=-50, to=-1, label=NULL)
table(PACdsPAS@anno$V1_gram)

## -----------------------------------------------------------------------------
PACdsPAS=annotateByPAS(PACds, bsgenome, 
                       grams=c('AATAAA','ATTAAA','GATAAA','AAAA'),
                       from=-50, to=-1, label='GRAM')
table(PACdsPAS@anno$GRAM_gram)

## -----------------------------------------------------------------------------
PACdsPAS=annotateByPAS(PACds, bsgenome, 
                       grams=c('AATAAA','ATTAAA','GATAAA','AAAA'),
                       priority=c(1,2,3,3), 
                       from=-50, to=-1, label='GRAM')
table(PACdsPAS@anno$GRAM_gram)

## ----message=FALSE, warning=FALSE---------------------------------------------
pas=PACdsPAS@anno$GRAM_gram[!is.na(PACdsPAS@anno$GRAM_gram)]
plotSeqLogo(pas)

## -----------------------------------------------------------------------------
v=getVarGrams('mm')
priority=c(1,2,rep(3, length(v)-2))

## -----------------------------------------------------------------------------
PACdsMM=annotateByPAS(PACds, bsgenome, grams=v, 
                      priority=priority, 
                      from=-50, to=-1, label='mm')

## ----results='hide'-----------------------------------------------------------
library(magrittr)
library(dplyr)
pas=as.data.frame(cbind(region=PACdsMM@anno$ftr, PAS=PACdsMM@anno$mm_gram))
pas$PAS[is.na(pas$PAS)]='NOPAS'
pas$PAS[pas$PAS %in% v[-c(1:2)]]='Variants'
n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
n=merge(n, n2)
n$PAC=n$nPAC/n$nTot
n=n[n$PAS!='NOPAS', ]
n$PAS=factor(n$PAS, levels=rev(c('AATAAA', 'ATTAAA','Variants', 'NOPAS')))
n$region=factor(n$region, 
                levels=c('3UTR','Ext_3UTR', 'intergenic','intron','CDS','5UTR'))

## ----fig.dim = c(6, 4), fig.align='left'--------------------------------------
library(ggplot2)
ggplot(data=n, aes(x=region, y=PAC, fill=PAS)) + 
  geom_bar(stat="identity") + 
  ylab("PAC Fraction") + theme_bw()

## ----eval=FALSE---------------------------------------------------------------
#  ## Extract the sequence of PACs, from UPA_start to UPA_end.
#  faFromPACds(PACds, bsgenome, what='pac', fapre='pac')
#  
#  ## Extract upstream 300 bp ~ downstream 100 bp around PACs,
#  ## where the position of PAC is 301.
#  faFromPACds(PACds, bsgenome, what='updn', fapre='updn',
#              up=-300, dn=100)
#  
#  ## Divide PACs into groups of genomic regions and then extract sequences for each group.
#  faFromPACds(PACds, bsgenome, what='updn', fapre='updn',
#              up=-100, dn=100, byGrp='ftr')
#  
#  ## Extract sequences for only 3UTR PACs.
#  faFromPACds(PACds, bsgenome, what='updn', fapre='updn',
#              up=-300, dn=100, byGrp=list(ftr='3UTR'))
#  
#  ## Extract sequences for only 3UTR PACs and separate sequences by strand.
#  faFromPACds(PACds, bsgenome, what='updn', fapre='updn',
#              up=-300, dn=100,
#              byGrp=list(ftr='3UTR', strand=c('+','-')))
#  
#  ## Extract sequences of genomic regions where PACs are located.
#  faFromPACds(PACds, bsgenome, what='region', fapre='region', byGrp='ftr')

## ----eval=FALSE---------------------------------------------------------------
#  ## The suggested signal regions when species is 'chlamydomonas_reinhardtii'.
#  files=faFromPACds(PACds, bsgenome, what='updn', fapre='Chlamy.NUE',
#                    up=-25, dn=-5, byGrp='ftr')
#  files=faFromPACds(PACds, bsgenome, what='updn', fapre='Chlamy.FUE',
#                    up=-150, dn=-25, byGrp='ftr')
#  files=faFromPACds(PACds, bsgenome, what='updn', fapre='Chlamy.CE',
#                    up=-5, dn=5, byGrp='ftr')
#  files=faFromPACds(PACds, bsgenome, what='updn', fapre='Chlamy.DE',
#                    up=-5, dn=30, byGrp='ftr')
#  
#  ## The suggested signal regions when species is plant.
#  ## In Arabidopsis or rice, signal regions are: FUE -200~-35, NUE -35~-10, CE -10~15.
#  files=faFromPACds(PACds, bsgenome, what='updn', fapre='plants.NUE',
#                    up=-35, dn=-10, byGrp='ftr')
#  files=faFromPACds(PACds, bsgenome, what='updn', fapre='plants.FUE',
#                    up=-200, dn=-35, byGrp='ftr')
#  files=faFromPACds(PACds, bsgenome, what='updn', fapre='plants.CE',
#                    up=-10, dn=15, byGrp='ftr')

## ----faFromPACds--------------------------------------------------------------
faFiles=faFromPACds(PACds, bsgenome, what='updn', fapre='updn', 
                    up=-300, dn=100, byGrp='ftr')

## ----plotATCGforFAfile, fig.dim=c(6,4)----------------------------------------
faFiles=c("updn.3UTR.fa", "updn.Ext_3UTR.fa", "updn.intergenic.fa", "updn.intron.fa")
## Plot single nucleotide profiles using the extracted sequences and merge all plots into one.
plotATCGforFAfile(faFiles, ofreq=FALSE, opdf=FALSE, 
                  refPos=301, mergePlots = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  plotATCGforFAfile (faFiles='updn.intron.fa',
#                     ofreq=FALSE, opdf=FALSE, refPos=301)
#  plotATCGforFAfile (faFiles='updn.intron.fa',
#                     ofreq=FALSE, opdf=FALSE, refPos=NULL,
#                     filepre ='NUE', start=250, end=350)

## ----eval=FALSE---------------------------------------------------------------
#  plotATCGforFAfile (faFiles, ofreq=TRUE, opdf=TRUE, refPos=301,
#                     filepre='singleBasePlot', mergePlots = TRUE)

## -----------------------------------------------------------------------------
## Count top 10 hexamers (k=6) in the NUE region 
## (normally from 265~295 if the PAC is at 301).
fafile='updn.3UTR.fa'
kcount(fafile=fafile, k=6, from=265, to=295, topn=10)

## Count given hexamers.
kcount(fafile=fafile, grams=c('AATAAA','ATTAAA'), 
       from=265, to=295, sort=FALSE)

## Count AATAAA and its 1nt variants in a given region.
kcount(fafile=fafile, grams='v1', from=265, to=295, sort=FALSE)

## ----message=FALSE------------------------------------------------------------
p=subsetPACds(PACds, group='group', pool=TRUE, totPACtag=20)

## -----------------------------------------------------------------------------
paShan=movPAindex(p, method='shan') 
## Show some rows with low H value (which means high overall tissue-specificity).
head(paShan[paShan$H<0.2742785, ], n=2)

## -----------------------------------------------------------------------------
paShan2=movPAindex(p, method='shan', shan.ratio = TRUE) 
head(paShan2, n=2)

## -----------------------------------------------------------------------------
paGeo=movPAindex(p, method='geo')
head(paGeo, n=2)

## -----------------------------------------------------------------------------
paRatio=movPAindex(p, method='ratio')
head(paRatio)

## ----fig.dim=c(6,6), message=FALSE, warning=FALSE-----------------------------
paShanHm=paShan[, -(1:3)]
paShanHm=paShanHm[rowSums(is.na(paShanHm))==0, ]
library(ComplexHeatmap, quietly = TRUE)
Heatmap(paShanHm, show_row_names=FALSE,  cluster_columns = FALSE, 
        heatmap_legend_param = list(title = 'Tissue\nspecificity'))

## -----------------------------------------------------------------------------
paShan=movPAindex(PACds, method='shan')

## ----fig.dim=c(6,6)-----------------------------------------------------------
## Plot heamap to show the consistency among replicates.
paShanHm=paShan[, -(1:3)]
paShanHm=paShanHm[rowSums(is.na(paShanHm))==0, ]
Heatmap(paShanHm, show_row_names=FALSE,  cluster_columns = TRUE, 
        heatmap_legend_param = list(title = 'Tissue\nspecificity'))

## ----eval=FALSE---------------------------------------------------------------
#  pd=get3UTRAPApd(pacds=p, minDist=50, maxDist=1000, minRatio=0.05, fixDistal=FALSE, addCols='pd')
#  rud=movAPAindex(pd, method="smartRUD", sRUD.oweight=TRUE)
#  head(rud$rud)
#  head(rud$weight)
#  geneRUD=rud$rud
#  geneRUD=geneRUD[rowSums(is.na(geneRUD))==0, ]
#  head(geneRUD, n=2)
#  Heatmap(geneRUD, show_row_names=FALSE,  cluster_columns = F,
#          heatmap_legend_param = list(title = 'RUD'))

## -----------------------------------------------------------------------------
geneWUL=movAPAindex(p, method="WUL", choose2PA=NULL)
head(geneWUL, n=2)

## ----eval=FALSE---------------------------------------------------------------
#  ## Remove NA rows before plotting heatmap.
#  geneWUL=geneWUL[rowSums(is.na(geneWUL))==0, ]
#  Heatmap(geneWUL, show_row_names=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  geneRUD=movAPAindex(p, method="RUD",
#                      choose2PA=NULL, RUD.includeNon3UTR=TRUE)
#  geneRUD=geneRUD[rowSums(is.na(geneRUD))==0, ]
#  head(geneRUD, n=2)
#  Heatmap(geneRUD, show_row_names=FALSE,  cluster_columns = F,
#          heatmap_legend_param = list(title = 'RUD'))

## ----eval=FALSE---------------------------------------------------------------
#  geneSLR=movAPAindex(p, method="SLR", choose2PA='PD')
#  head(geneSLR, n=2)
#  geneSLR=geneSLR[rowSums(is.na(geneSLR))==0, ]
#  Heatmap(geneSLR, show_row_names=FALSE)

## ----fig.dim=c(6,6)-----------------------------------------------------------
geneGPI=movAPAindex(PACds, method="GPI", choose2PA='PD')
head(geneGPI)
geneGPI=geneGPI[rowSums(is.na(geneGPI))==0, ]
Heatmap(geneGPI, show_row_names=FALSE,  cluster_columns = TRUE, 
        heatmap_legend_param = list(title = 'GPI'))

## ----results='hide', message=FALSE--------------------------------------------
library(DESeq2)
## Subset two conditions first.
pacds=subsetPACds(PACds, group='group', cond1='anther', cond2='embryo')
## Detect DE genes using DESeq2 method, 
## only genes with total read counts in all samples >=50 are used.
DEgene=movDEGene(PACds=pacds, method='DESeq2', group='group', minSumPAT=50)

## -----------------------------------------------------------------------------
stat=movStat(object=DEgene, padjThd=0.05, valueThd=0.5)
stat$nsig
head(stat$siglist$anther.embryo)

## ----results='hide', message=FALSE--------------------------------------------
DEgene=movDEGene(PACds=PACds, method='DESeq2', group='group', minSumPAT=50)
stat=movStat(object=DEgene, padjThd=0.05, valueThd=1)

## -----------------------------------------------------------------------------
## Number of DE genes in each pair of conditions.
stat$nsig
## Overlap between condition pairs.
stat$ovp

## ----eval=FALSE---------------------------------------------------------------
#  outputHeatStat(heatStats=stat, ostatfile='DEgene.stat', plotPre='DEgene')

## ----message=FALSE------------------------------------------------------------
selFull=movSelect(DEgene, condpair='embryo.anther', padjThd=0.05, valueThd=1, 
                  out='full', PACds=PACds)
head(selFull)

## -----------------------------------------------------------------------------
sel=movSelect(DEgene, condpair='anther.embryo', 
              padjThd=0.05, valueThd=1, out='pv')
head(sel)

## -----------------------------------------------------------------------------
sel=movSelect(DEgene, condpair='embryo.anther', 
              padjThd=0.05, upThd=0.5, out='gene')
head(sel)

## ----results='hide', message=FALSE--------------------------------------------
DEPAC=movDEPAC(PACds, method='DESeq2', group='group', minSumPAT=20)
DEXPAC=movDEPAC(PACds, method='DEXseq', group='group', minSumPAT=20)
DEqPAC=movDEPAC(PACds, method='chisq', group='group', minSumPAT=20)

## ----results='hide'-----------------------------------------------------------
library(ggplot2)
## Get significant DE results.
stat1=movStat(object=DEPAC, padjThd=0.05, valueThd=1)
stat2=movStat(object=DEXPAC, padjThd=0.05, valueThd=1)
stat3=movStat(object=DEqPAC, padjThd=0.05, valueThd=0.95)

## ----fig.dim = c(6, 4)--------------------------------------------------------
## Count the number of DE PACs by different methods.
nsig=as.data.frame(cbind(stat1$nsig, stat2$nsig, stat3$nsig))
colnames(nsig)=c('DESeq2','DEXseq','Chisq.test')
nsig$tissueA.tissueB=rownames(nsig)
nsig
## Plot a barplot.
nsig=reshape2::melt(nsig, variable.name='Method')
ggplot(data=nsig, aes(x=tissueA.tissueB, y=value, fill=Method)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("DE PAC#") + theme_bw()

## ----eval=FALSE---------------------------------------------------------------
#  ## First subset PACs in two conditions.
#  PACds1=subsetPACds(PACds, group='group',
#                     cond1='anther', cond2='embryo', choosePA='apa')
#  ## Detect DE PACs.
#  DEPAC1=movDEPAC(PACds1, method='DESeq2', group='group', minSumPAT=10)
#  DEXPAC1=movDEPAC(PACds1, method='DEXseq', group='group', minSumPAT=10)
#  DEqPAC1=movDEPAC(PACds1, method='chisq', group='group', minSumPAT=10)

## ----results='hide'-----------------------------------------------------------
stat=movStat(object=DEPAC, padjThd=0.05, valueThd=1)

## -----------------------------------------------------------------------------
## Number of DE PACs between conditions.
stat$nsig
## Overlap of DE PACs between different pairs of conditions.
head(stat$ovp)
## DE PAC list
head(stat$siglist[[1]])

## ----eval=FALSE---------------------------------------------------------------
#  library(VennDiagram, quietly = TRUE)
#  x=venn.diagram(stat$siglist, fill=brewer.pal(3, "Set1"), cex=2,
#                 cat.fontface=4,  filename='DEPAC.venn')

## ----eval=FALSE---------------------------------------------------------------
#  stat=movStat(object=DEqPAC, padjThd=0.05, valueThd=0.95)

## ----message=FALSE, warning=FALSE---------------------------------------------
## Here method is DEXseq, so the valueThd (log2FC) threshold is automatelly determined. 
sel=movSelect(aMovRes=DEXPAC, condpair='embryo.anther', 
              padjThd=0.1, out='full', PACds=PACds)
head(sel, n=2)

## You can also mannually set a log2FC threshold.
sel=movSelect(aMovRes=DEXPAC, condpair='embryo.anther', 
              padjThd=0.1, valueThd=2, out='pa');
head(sel)

## Filter only up-regulated PACs in embryo 
## (value=log2(embryo_this_others/anther_this_others)).
sel=movSelect(aMovRes=DEXPAC, condpair='embryo.anther', 
              padjThd=0.1, upThd=2, out='full', PACds=PACds)
head(sel, 2) 

## ----eval=FALSE---------------------------------------------------------------
#  ## Filter less results and plot the heatmap clearly.
#  stat=movStat(object=DEPAC, padjThd=0.001, valueThd=8)
#  outputHeatStat(heatStats=stat, ostatfile='DEPAC.stat', plotPre='DEPAC',
#                 show_rownames = TRUE)

## -----------------------------------------------------------------------------
gene='Os05g0541900'
gp=PACds[PACds@anno$gene==gene, ]
cbind(gp@anno$ftr, rowSums(gp@counts))

## ----fig.dim=c(7, 5), message=FALSE, warning=FALSE----------------------------
movViz(object=DEPAC, gene=gene, txdb=gff, PACds=PACds, collapseConds=FALSE, 
       padjThd=0.01, showRatio=FALSE, showAllPA=TRUE) 

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
movViz(object=DEPAC, gene=gene, txdb=gff, PACds=PACds, collapseConds=TRUE, 
       padjThd=0.01, showPV=TRUE, showAllPA=FALSE, showRatio=F,
       conds=DEPAC@conds[c(1,3), ], highlightConds=DEPAC@conds[c(3), ])

## ----message=FALSE------------------------------------------------------------
utr=movUTRtrend(PACds, group='group', method='linearTrend', 
                avgPACtag=10, avgGeneTag=20)
## Number of genes for analyzing, including those not significant.
lapply(utr@fullList, nrow)
head(utr@fullList[["anther.embryo"]], n=2)

## -----------------------------------------------------------------------------
stat=movStat(object=utr, padjThd=0.1, valueThd=0)
stat$nsig

## ----eval=FALSE---------------------------------------------------------------
#  ## Only output gene ids.
#  out=movSelect(aMovRes=utr, condpair='anther.embryo',
#                padjThd=0.1, valueThd=0, out='gene')
#  ## Output PAC ids.
#  out=movSelect(aMovRes=utr, condpair='anther.maturePollen',
#                padjThd=0.1, valueThd=0, out='pa')
#  ## Output gene ids with padj and value.
#  out=movSelect(aMovRes=utr, condpair='anther.embryo',
#                padjThd=0.1, valueThd=0, out='pv')
#  ## Output full information with expression levels, 3UTR length,
#  ## read counts of each PA in each sample, etc.
#  out=movSelect(aMovRes=utr, condpair='anther.embryo',
#                padjThd=0.1, valueThd=0, out='full')
#  ## Output full information for 3UTR lengthening genes from anther to embryo (change=1).
#  out=movSelect(aMovRes=utr, condpair='anther.embryo',
#                padjThd=0.1, upThd=0, out='full')

## -----------------------------------------------------------------------------
## Output full information for 3UTR shortening genes from anther to embryo (change=-1).
out=movSelect(aMovRes=utr, condpair='anther.embryo', 
              padjThd=0.1, dnThd=0, out='full')
head(out, n=2)

## ----results='hide', message=FALSE--------------------------------------------
DEXPAC=movDEPAC(PACds, method='DEXseq', group='group', minSumPAT=10)
swDEX=movUTRtrend(PACds, group='group', method='DEX',
                   avgPACtag=10, avgGeneTag=20,
                   aMovDEPACRes=DEXPAC, DEPAC.padjThd=0.01,
                   mindist=50, fisherThd=0.01, logFCThd=1, selectOne='farest')

## ----message=FALSE------------------------------------------------------------
stat=movStat(object=swDEX, padjThd=0.01, valueThd=1)
stat$nsig

out=movSelect(aMovRes=swDEX, condpair='anther.embryo', 
              padjThd=0.01, valueThd=1, out='full')
head(out, n=2)

## ----results='hide'-----------------------------------------------------------
swLinear=movUTRtrend(PACds, group='group',method='linearTrend',
                     avgPACtag=10, avgGeneTag=20)
swDEX=movUTRtrend(PACds, group='group', method='DEX',
                 avgPACtag=10, avgGeneTag=20,
                 aMovDEPACRes=DEXPAC, DEPAC.padjThd=0.01,
                 mindist=50, fisherThd=0.01, logFCThd=1, selectOne='fisherPV')

swDE=movUTRtrend(PACds, group='group', method='DE',
                  avgPACtag=10, avgGeneTag=20,
                  aMovDEPACRes=DEPAC, DEPAC.padjThd=0.01,
                  mindist=50, fisherThd=0.01, logFCThd=1, selectOne='fisherPV')

## ----results='hide'-----------------------------------------------------------
stat1=movStat(object=swLinear, padjThd=0.1, valueThd=0)
stat2=movStat(object=swDEX, padjThd=0.01, valueThd=1)
stat3=movStat(object=swDE, padjThd=0.01, valueThd=1)

## ----fig.dim=c(6, 4)----------------------------------------------------------
nsig=as.data.frame(cbind(stat1$nsig, stat2$nsig, stat3$nsig))
colnames(nsig)=c('LinearTrend','DE-DEXseq','DE-DESeq')
nsig$tissueA.tissueB=rownames(nsig)
nsig
nsig=reshape2::melt(nsig, variable.name='Method')
ggplot(data=nsig, aes(x=tissueA.tissueB, y=value, fill=Method)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("3\'UTR switching events #") + theme_bw()


## -----------------------------------------------------------------------------
gene='Os02g0759700'
gp=PACds[PACds@anno$gene==gene, ]
cbind(gp@anno$ftr, rowSums(gp@counts))

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
movViz(object=swDE, gene=gene, txdb=NULL, PACds=PACds, showRatio=TRUE, 
       padjThd=0.01, showAllPA=TRUE)

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
movViz(object=swDE, gene=gene, txdb=gff, PACds=PACds, collapseConds=TRUE,
       conds=swDE@conds, highlightConds=swDE@conds[c(1,3), ], showRatio=TRUE,
       linkPAs=TRUE, padjThd=0.01, showAllPA=FALSE, showPV=TRUE)

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
movViz(object=swDE, gene=gene, txdb=NULL, PACds=PACds, collapseConds=TRUE,
       conds=swDE@conds[1, ], highlightConds=NULL, showRatio=TRUE, linkPAs=TRUE,
       padjThd=0.01, showAllPA=FALSE, showPV=FALSE)

## ----message=FALSE------------------------------------------------------------
stat=movStat(object=swDE, padjThd=0.01, valueThd=1)
stat$nsig

## ----eval=FALSE---------------------------------------------------------------
#  outputHeatStat(heatStats=stat, ostatfile='3UTR_switching_DE.stat',
#                 plotPre='3UTR_switching_DE', show_rownames = TRUE)

## -----------------------------------------------------------------------------
heat=movRes2heatmapResults(swDE)
heatUp=subsetHeatmap(heat, padjThd=0.05, valueThd=1)

## ----fig.dim=c(6, 6)----------------------------------------------------------
plotHeatmap(heatUp@value, show_rownames=TRUE, plotPre=NULL, cluster_rows=TRUE) 

## -----------------------------------------------------------------------------
fl=swDE@fullList$anther.embryo
fl[fl$gene=='Os06g0682633',]

## ----eval=FALSE---------------------------------------------------------------
#  DEXPAC=movDEPAC(PACds, method='DEXseq', group='group', minSumPAT=10)

## ----message=FALSE, results='hide'--------------------------------------------
swDEX=movAPAswitch(PACds, group='group',aMovDEPACRes=DEXPAC,
                   avgPACtag=5, avgGeneTag=10,
                   only3UTR=TRUE,
                   DEPAC.padjThd=0.1, nDEPAC=1,
                   mindist=50, fisherThd=0.1, logFCThd=0.5, 
                   cross=FALSE, selectOne=NULL)

## ----message=FALSE------------------------------------------------------------
stat=movStat(object=swDEX, padjThd=0.1, valueThd=1)
stat$nsig

## ----message=FALSE------------------------------------------------------------
sel=movSelect(aMovRes=swDEX, condpair='anther.embryo', 
              padjThd=0.1, valueThd=1, out='full')
head(sel, n=2)

## ----message=FALSE, results='hide'--------------------------------------------
swDE=movAPAswitch(PACds, group='group', aMovDEPACRes=DEXPAC,
                   avgPACtag=10, avgGeneTag=20,
                   only3UTR=FALSE,
                   DEPAC.padjThd=0.1, nDEPAC=1,
                   mindist=50, fisherThd=0.1, logFCThd=0.5, 
                  cross=FALSE, selectOne=NULL)

## ----message=FALSE------------------------------------------------------------
stat=movStat(object=swDE, padjThd=0.1, valueThd=1)
stat$nsig

## ----message=FALSE------------------------------------------------------------
sw=movSelect(aMovRes=swDE, condpair='anther.embryo', 
             padjThd=0.01, valueThd=1, out='full')
head(sw[order(sw$fisherPV), ], n=2)

## ----message=FALSE------------------------------------------------------------
genes=movSelect(aMovRes=swDE, condpair='anther.embryo', 
                padjThd=0.01, valueThd=1, out='gene')
swPAC=subsetPACds(PACds, genes=genes, verbose=TRUE)
table(swPAC@anno$ftr)

PAs=movSelect(aMovRes=swDE, condpair='anther.embryo', padjThd=0.01,
              valueThd=1, out='pa')
swPAC=subsetPACds(PACds, PAs=PAs, verbose=TRUE)
table(swPAC@anno$ftr)

## -----------------------------------------------------------------------------
gene='Os05g0451900'
gp=PACds[PACds@anno$gene==gene, ]
cbind(gp@anno$ftr, rowSums(gp@counts))

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
movViz(object=swDE, gene=gene, txdb=gff, PACds=PACds, 
       showRatio=TRUE, padjThd=0.01, showAllPA=TRUE)

## ----fig.dim=c(7, 5), message=FALSE-------------------------------------------
movViz(object=swDE, gene=gene, txdb=gff, PACds=PACds, collapseConds=TRUE,
       conds=swDE@conds, highlightConds=swDE@conds[c(2,3), ], showRatio=TRUE,
       linkPAs=TRUE, padjThd=0.01, showAllPA=FALSE)

## -----------------------------------------------------------------------------
sessionInfo()

