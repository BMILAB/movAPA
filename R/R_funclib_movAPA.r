#' @importFrom GenomicRanges GRanges reduce resize mcols "mcols<-" findOverlaps
#' follow precede width start end strand seqnames gaps
#' @importFrom IRanges IRanges gaps
#' @importFrom ggplot2 element_rect element_line aes theme_bw theme element_text ggplot geom_line xlab ylab
#' element_blank ggtitle facet_grid facet_wrap aes_string stat_ecdf labs annotate geom_point geom_vline
#' scale_colour_manual labs margin scale_y_continuous theme_classic scale_fill_gradientn guides .data
#' geom_text position_dodge
#' @importFrom magrittr %>%
#' @importFrom methods setClassUnion
#' @importFrom Matrix as.matrix colMeans colSums rowMeans rowSums sparseMatrix t
#' @importFrom grDevices boxplot.stats colorRampPalette dev.off pdf
#' @importFrom methods as new slot
#' @importFrom stats quantile aggregate as.formula chisq.test coef fisher.test na.omit p.adjust
#' pnorm prcomp prop.trend.test relevel rnorm runif sd
#' @importFrom utils combn read.table write.csv write.table
#' @importClassesFrom Matrix dgCMatrix dgeMatrix
#' @importFrom utils head
NULL

options(stringsAsFactors=F)

# ------------------- Results desc-------------------------------
# DESeq2:
# change of absolute expression
# 1200.0140.posup: under 1200/0140@pos=up, log2(1200/0140)

# DEXseq:
# change of relative expression
# 1200.0140: under 1200/0140, log2(1200_this_others/0140_this_others)
# 1200: under 1200, relative exon usage (the same trend as 1200_this_others/0140_this_others, but values are different)
# all condition pairs have just one padj, not each pair one padj!
# get FC between a condition pair by estimatefoldchange, but use getDEXlist() to get DElist

# UTR-linearTrend:
# 1200.0140, cor is from the FU(), if cor<0 <==> change=-1 <==> shorter in 0140 (avgUTRlen1 > avgUTRlen2)

# UTR-switching:
# 1200.0140 is logFC (as in TAPAS), the log2FC of PAij and sample12; logFC>0 <==> change=1 longer in 0140 (avgUTRlen1 < avgUTRlen2)

# For example:
# cond1: PH02Gene07028.PA4=22;PH02Gene07028.PA5=487;
# cond2: PH02Gene07028.PA4=13;PH02Gene07028.PA5=552;
# In switching, logFC=log2((22/487)/(13/552))=0.9397384 (change=1; cond1 to cond2 longer)

# Betas:
# 0350.0140.posdn means under the control (dn), the diff between 0350 and 0140;
# 0350.0140.posmd_vs_dn means, compared to dn, the *change of difference of 350-140* under md (>0 change to higher level, <0 change to less expression)

# ------------------- data -------------------------------
## package?movAPA
#' @details
#' The movAPA package is a comprehensive package for MOdeling and Visualization of dynamics of Alternative PolyAdenylation across biological samples in different species.
#' Functions of movAPA include pre-processing of poly(A) site clusters (PACs), poly(A) signal analysis, differential analyses, APA site switching analyses, etc.
#' The package is not species-specific, so it can be used for any species (providing qualified GFF annotation and genome references).
#' @references
#' Ye, W., et al. (2021) movAPA: Modeling and visualization of dynamics of alternative polyadenylation across biological samples, Bioinformatics, 37, 2470-2472.
#'
#' Zhu, S., et al. (2020) PlantAPAdb: A Comprehensive Database for Alternative Polyadenylation Sites in Plants, Plant Physiol., 182, 228-242.
#'
#' Wu, X., et al. (2011) Genome-wide landscape of polyadenylation in Arabidopsis provides evidence for extensive alternative polyadenylation, Proc. Natl. Acad. Sci. USA, 108, 12533-12538.
#' @docType package
#' @name movAPA-package
#' @rdname movAPA-package
#'
"_PACKAGE"


#' A demo PACdataset of rice tissues
#'
#' A dataset of the PACdataset class from rice japonica, containing 1233 PACs in 455 genes from three tissues
#' Biological samples are anther, embryo, and maturePollen, each with three replicates.
#' PAC counts have been normalized.
#'
#' @format A PACdataset with 1233 PACs.
#' @source <Zhu et al, 2019. Plant Physiol>
#' @family PACdataset functions
"PACds"

#' A demo PACdataset of mouse sperm cells
#'
#' A dataset of the PACdataset class from mouse sperm cells, containing 771 PACs from 396 genes located in chromosome 12.
#' There are total 2042 cells from three cell types (SC, Spermatocytes; RS, Round spermatids; ES, Elongating spermatids).
#' This dataset contains the gene Psen1 (ENSMUSG00000019969) presented in Shulman et al, 2019.
#' The original PAC dataset was obtained from Shulman et al, 2019, but re-annotated by the MM10 genome annotation (TxDb.Mmusculus.UCSC.mm10.ensGene).
#'
#' @format A PACdataset with 771 PACs from single cells.
#' @source <Shulman et al, 2019. Nucleic Acids Res>
#' @family PACdataset functions
"scPACds"

#' A demo peaks list of human
#'
#' demo_peaks stores a list containing 5w peaks in 2538 cells, including peaks.meta and peaks.demo. This file was obtained by scAPAtrap.
#'
#' @format A list including peaks.meta and peaks.demo.
#' @source One result from scAPAtrap.
"demo_peaks"


# ------------------- common funcs -------------------------------
.saveMtx2File<-function(mtx, file, append=T, row.names=F, col.names=T, sep="\t", quote=F) {
  write.table(mtx, file=file, append=append, col.names=col.names, row.names=row.names, quote=quote,sep=sep)
}

# Write log
#
# `LOG` writes log (text or table) to a file and print txt on the screen.
#
# This function Will automatelly add new line to a text.
# If text=NULL, then output ************************** and new line.
# if star=TRUE, then output star***+txt+star***
# @usage LOG(txt=NULL, file=NULL, append=T, star=FALSE)
# @param txt A character.
# @param file Output to a file, if NULL then output to the screen.
# @param append Append to a file.
# @param star Add **** before and after the txt
# @return NULL
# @examples
# ## write a header line to a.txt and print it on the screen
# LOG('write a header line', file='a.txt', append=FALSE, star=TRUE)
# ## write a line of star (*) to a.txt
# LOG(txt=NULL, file='a.txt')
# ## print a line on the screen
# LOG(txt='hello', file=NULL)
# @name LOG
LOG <- function(txt=NULL, file=NULL, append=T, star=FALSE){
  if (is.null(txt)) {
    write(paste0(rep('*',50),collapse = ''),file=file,append=append)
    return()
  }
  if (is.data.frame(txt) | is.matrix(txt) ) {
    if (!is.null(file)) {
      .saveMtx2File(txt,file,append, row.names=!is.null(rownames(txt)))
    }
    cat('save matrix [',colnames(txt),']\n')
  } else {
    if (!is.null(file)) {
      if (star) {
        if (append) {
          write("\n",file=file,append=append)
        }
        LOG(txt=NULL,file=file,append=append)
      }
      write(txt,file=file,append=append | star)
      if (star) {
        LOG(txt=NULL,file=file,append=TRUE)
      }
    }
    cat(txt,'\n')
  }
}

# Output all coef matrices from DEseq results
.outputResultsDDS <- function (dds, odir, sufix="") {
  dir.create(odir, showWarnings = FALSE)
  cat('>>> dir',odir,'\n')
  for (rs in DESeq2::resultsNames(dds)) {
    cat('>>>',rs,'\n')
    res=DESeq2::results(dds, name=rs, cooksCutoff=FALSE, independentFiltering=FALSE)
    write.csv(res, file=paste0(odir,'/',rs,sufix,'.csv'))
  }
}

# get specified columns of all data frames in a list, each column generates a single data frame
# @param resList A list, each element is one dataframe
# @param cols A column in each of all dfs in resList
# @Return A list[log2FC=aDf; padj=aDf] each item in cols generates one item in list
.getStatsFromResList<-function(resList, cols=c('log2FC'='log2FoldChange','padj'='padj','pvalue'='pvalue')) {
  ls=list()
  if (is.null(names(cols))) {
    names(cols)=cols
  }
  for (cl in cols) {
    mtx=lapply(resList,'[',cl)
    mtx=.asDf(mtx)
    colnames(mtx)=names(resList)
    ls[[names(cols)[which(cols==cl)]]]=mtx
  }
  return(ls)
}

# Whether A in B
# AinB(c('start','end'),c('start','x','y'),all=T)
# @param A a vector
# @param B a vector
# @param all whether need all A in B
# @return TRUE/FALSE
# @examples
# AinB(c('start','end'),c('start','x','y'),all=T)
AinB<-function(A, B, all=T){
  if (is.factor(A))   A=as.character(A)
  if (is.factor(B))   B=as.character(B)
  if (!(is.vector(A) & is.vector(B))) {
    return(F)
  }
  x=sum(A %in% B)
  if ((all & x==length(A)) | (!all & x>0)) {
    return(T)
  } else {
    return(F)
  }
}


## .autoDetectCommonString(c('xxACyy','xBy','xxBCyy')) --> c(x,y)
## .autoDetectCommonString(c('xxACyy','xxByy','xxBCC')) --> c(xx,yy)
## .autoDetectCommonString(c('xxACyy','xxByy','xxBCC'), sbj='xxByyyy') --> Byy
## .autoDetectCommonString(c('xxACyy','xBy','xxBC')) --> c(x,y)
## .autoDetectCommonString(c('x.xACyy','x.By','x.xBC')) --> c(x.,y)
## .autoDetectCommonString(c('x.xACyy','x.By','x.xBC'), beAll=FALSE) --> c(x.,y)
## .autoDetectCommonString(c('x.xACyy','x.By','x.xBC'), beAll=TRUE) --> c(x.,'')
## .autoDetectCommonString(c('AC','DB','EF')) --> c('','')
## beAll=TRUE, then cmchars must be in all strs (not just two strs)
.autoDetectCommonString<-function(strs, sbj=NULL, beAll=FALSE) {
  cms=''
  cme=''
  for (i in 1:(length(strs)-1)) {
    for (j in (i+1):length(strs)) {
      s1=strs[i]; s2=strs[j];
      #str1=s1; str2=s2;
      s1=unlist(strsplit(s1,split=''))
      s2=unlist(strsplit(s2,split=''))
      len=min(length(s1),length(s2))
      same=s1[1:len]==s2[1:len]

      newcms=''
      cmstart=which(!same)
      if (length(cmstart)==0) {
        newcms=paste(s1[1:len],collapse='')
      }
      if (length(cmstart)>0 & cmstart[1]>1) {
        cmstart=cmstart[1]-1
        if (cmstart>=1) {
          newcms=paste(s1[1:cmstart],collapse='')
        }
      }

      if (cms!='' & newcms!='' & newcms!=cms) {
        cms=.autoDetectCommonString(c(cms,newcms))[1]
      } else if (cms=='') {
        cms=newcms
      }

      same=s1[(length(s1)-len+1):length(s1)]==s2[(length(s2)-len+1):length(s2)]
      s1=s1[(length(s1)-len+1):length(s1)]
      cmstart=which(!same)
      newcms=''
      if (length(cmstart)==0) {
        newcms=paste(s1,collapse='')
      }
      if (length(cmstart)>0 ) {
        cmstart=cmstart[length(cmstart)]+1
        if (cmstart<=length(s1)) {
          newcms=paste(s1[cmstart:length(s1)],collapse='')
        }
      }

      if (cme!='' & newcms!='' &  newcms!=cme) {
        cme=.autoDetectCommonString(c(cme,newcms))[2]
      } else if (cme=='') {
        cme=newcms
      }
    }
  }

  if (beAll) {
    n1=sum(grepl(paste0('^',cms), strs))
    if (n1!=length(strs)) cms=''
    n1=sum(grepl(paste0(cme,'$'), strs))
    if (n1!=length(strs)) cme=''
  }

  if (is.null(sbj)) {
    return(c(cms,cme))
  } else {
    if (cms!='') {
      sbj=gsub(paste('^',cms,sep=''),'',sbj)
    }
    if (cme!='') {
      sbj=gsub(paste(cme,'$',sep=''),'',sbj)
    }
    return(sbj)
  }
}

#cms=.autoDetectCommonString(c('xxACyy','xxByy','xxBCC')) --> xx yy
#.removeCommonStr(c('xxACyy','xxByy','xxBCC'),cms) --> AC B BCC
.removeCommonStr<-function(sbj,cms,fixed=TRUE){
  if (fixed)  {
    if (length(cms)!=2) {
      stop(".removeCommonStr: cms not length 2")
    }
    if (cms[1]!='') {
      sbj=gsub(paste('^',cms[1],sep=''),'',sbj)
    }
    if (cms[2]!='') {
      sbj=gsub(paste(cms[2],'$',sep=''),'',sbj)
    }
    return(sbj)
  }
  for (cm in cms) {
    if (cm!='') {
      sbj=gsub(cm,'',sbj)
    }
  }
  return(sbj)
}



# -------- class PACdataset -------------
.isAnyMatrix<-function(data) {
  if (is.matrix(data)) return(TRUE)
  if (inherits(x = data, what = 'dgCMatrix')) return(TRUE)
  return(FALSE)
}


## if not matrix or 0%>0.5 or colN<=50, not sparse, return F
.isSparse<-function(data, perc=0.5, colN=50) {
  if (inherits(x = data, what = 'dgCMatrix')) return(TRUE)
  if (!.isAnyMatrix(data) & !is.data.frame(data)) return(FALSE)
  if (nrow(data)==0 | ncol(data)==0) return(FALSE)
  if (ncol(data)<=colN) return(FALSE)
  n0=sum(is.na(data) | data==0)
  N=nrow(data)*ncol(data)
  if (n0/N>=perc) return(TRUE)
  return(FALSE)
}


#' Convert to matrix or dgCMatrix
#'
#' asAnyMatrix convert data.frame/matrix/dgCMatrix to matrix or dgCMatrix.
#' For PACdataset@counts in previous movAPA, the data type is data.frame. Please call asAnyMatrix(PACds@counts) to make it compatible to new movAPA.
#' If dat is data-sparse or type-sparse, then return dgCMatrix
#' if dat is data.frame, then return matrix or dgCMatrix (when forceSparse or data-sparse)
#'
#' @param dat a data.frame/matrix/dgCMatrix or any object that can be converted to matrix using as.matrix.
#' @param forceSparse if TRUE, then force the dat to sparse matrix (dgCMatrix)
#' @return A matrix or dgCMatrix which can be used as the count slot of a PACdataset.
#' @name asAnyMatrix
#' @export
asAnyMatrix<-function(dat, forceSparse=FALSE) {
  if (inherits(x = dat, what = 'dgCMatrix')) return(dat)
  if (is.matrix(dat) & forceSparse) {
    if (mode(dat)=='character') { #convert char matrix
      dat <- matrix( as.numeric(dat), ncol = ncol(dat))
    }
    return(as(dat, 'dgCMatrix'))
  }
  if (is.data.frame(dat)) {
    dat=as.matrix(dat)
  }
  if (mode(dat)=='character') { #convert char matrix
    dat <- matrix( as.numeric(dat), ncol = ncol(dat), dimnames = dimnames(dat))
  }
  if (forceSparse | .isSparse(dat)) {
    return(as(as.matrix(dat), 'dgCMatrix'))
  }
  return(as.matrix(dat))
  #stop("asAnyMatrix error, dat not dgCMatrix/matirx/data.frame")
}

.asDf <-function(dat) {
  if (inherits(x = dat, what = 'dgCMatrix')) return(as.data.frame(as.matrix(dat), stringsAsFactors=FALSE))
  if (inherits(x = dat, what = 'dgeMatrix')) return(as.data.frame(as.matrix(dat), stringsAsFactors=FALSE))
  return( as.data.frame(dat, stringsAsFactors=FALSE) ) #tibble->df
}

## counts type
## @export AnyMatrix
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))

#' PACdataset object
#'
#' PACdataset is an S4 class to represent a poly(A) site list.
#' PACs are stored as PACdataset for filtering, statistics, and other operations.
#'
#' @slot counts a data frame denoting the counts of PAs. Each row is one PA, each column is one experiment.
#' @slot colData sample defination. Each row is one experiment, and the column is the sample group.
#' The number of rows is the same as the column number of the counts slot.
#' @slot anno detailed annotation information for PAs, with the same row number as the counts slot.
#' @slot supp a list storing additional data, like stopCodon.
#'
#' @name PACdataset-class
#' @family PACdataset functions
#' @exportClass PACdataset
PACdataset <- setClass(Class = "PACdataset",
                       slots=list(counts="AnyMatrix", anno="data.frame", colData="data.frame", supp="list"))


setMethod("show",
          "PACdataset",
          function(object) {
            cat('PAC#',nrow(object@counts),'\n')

            if ('ftr' %in% colnames(object@anno)) {
              cat('gene#',length(unique(object@anno$gene[getNonItgFtrId(object@anno$ftr)])),'\n')
              pn=table(object@anno$ftr)
              pn=cbind(nPAC=pn)
              print(pn)

              if ('3UTR' %in% unique(object@anno$ftr)) {
                avgutr=NA
                if ('toStop' %in% colnames(object@anno)) {
                  avgutr=mean(object@anno$toStop[object@anno$ftr=='3UTR'])
                } else if ('three_UTR_length' %in% colnames(object@anno)) {
                  avgutr=mean(object@anno$three_UTR_length[object@anno$ftr=='3UTR'])
                }
                if (!is.na(avgutr)) cat('Mean 3UTR length of PACs (bp):',floor(avgutr),'\n')
              }
            }

            cat('sample#',ncol(object@counts),'\n')
            cat(colnames(object@counts)[1:(min(ncol(object@counts),5))],'...\n')
            cat('groups:\n')
            cat(sprintf('%s...[%d x %d]%s','@colData',nrow(object@colData), ncol(object@colData),'\n'))
            if (nrow(object@colData)>0)
            print(head(object@colData[, 1:min(10, ncol(object@colData)), drop=F], 2))
            cat(sprintf('%s...[%d x %d]%s','@counts',nrow(object@counts), ncol(object@counts),'\n'))
            if (nrow(object@counts)>0)
            print(head(object@counts[, 1:min(10, ncol(object@counts)), drop=F],2))
            cat(sprintf('%s...[%d x %d]%s','@colData',nrow(object@colData), ncol(object@colData),'\n'))
            print(head(object@colData,2))
            cat(sprintf('%s...[%d x %d]%s','@anno',nrow(object@anno), ncol(object@anno),'\n'))
            print(head(object@anno,2))

            if (length(object@supp)>0) {
              cat(sprintf('%s...[%d]%s','@supp',length(object@supp),'\n'))
              cat(names(object@supp),'\n')
            }

          }
)

#' Summary of a PACdataset
#'
#' @param object a PACdataset
#' @return Print the summary information of a PACdataset.
#' @name summary
#' @family PACdataset functions
#' @export
setMethod("summary",
          "PACdataset",
          function(object) {
            cat('PAC#',nrow(object@counts),'\n')
            cat('sample#',ncol(object@counts),'\n')
            cat('summary of expression level of each PA\n')
            print(summary(rowSums(object@counts)))
            cat('summary of expressed sample# of each PA\n')
            print(summary(rowSums(object@counts>0)))

            if (all( c('gene', 'ftr') %in% colnames(object@anno))) {
              cat('gene#',length(unique(object@anno$gene[getNonItgFtrId(object@anno$ftr)])),'\n')
            } else if ('gene' %in% colnames(object@anno)) {
              cat('gene#',length(unique(object@anno$gene)),'\n')
            }

            if ('ftr' %in% colnames(object@anno)) {
              pn=table(object@anno$ftr)
              pn=cbind(nPAC=pn)
              print(pn)
            }
          }
)


#' Length of a PACdataset
#'
#' @param x a PACdataset
#' @return The length, i.e., the number of pAs, i.e., the number of rows in slot counts or anno.
#' @name length
#' @family PACdataset functions
#' @export
setMethod("length",
          "PACdataset",
          function(x) {
            return(nrow(x@counts))
          }
)

#' Subset rows and/or columns of a PACdataset
#'
#' Subset rows of a PACdataset, by overloading of subscript operator.
#'
#' @param x a PACdataset
#' @return A new PACdataset.
#' @family PACdataset functions
#' @examples
#' data(PACds)
#' PACds[1:100, 2:5]
#' PACds[1:100, ]
#' PACds[, 1:5]
#' PACds[, c(T,T,F)]
#' PACds[, c('anther1','anther2')]
#' @name subscript_operator
#' @export
setMethod("[", signature(x = "PACdataset"), def=function(x, i, j, ..., drop=FALSE) {
  if (!missing(i)) {
    x@counts=x@counts[i, , drop=drop]
    x@anno=x@anno[i,, drop=drop ]
  }
  if (!missing(j)) {
    x@counts=x@counts[, j, drop=drop ]
    cnames=colnames(x@counts)
    x@colData=x@colData[cnames, , drop=F]
  }
  for (i in 1:ncol(x@colData)) {
    x@colData[,i]=factor(x@colData[,i])
  }
  return(x)
}
)



## rbind multiple PACds
rbind.PACds <- function(...) {
  pars=list(...)
  if (length(pars)==0) return(pars)
  if (length(pars)==1) return(pars[[1]])
  if (sum(lapply(pars, class)=='PACdataset')!=length(pars)) {
    stop("rbind: not all vars are PACdataset!\n")
  }
  p=pars[[1]]
  for (i in 2:length(pars)) {
    p2=pars[[i]]
    if (ncol(p@counts)!=ncol(p2@counts)) stop("rbind: not same column number in @counts!")
    if (ncol(p@anno)!=ncol(p2@anno)) stop("rbind: not same column number in @anno!")
    if (!identical(p@colData,p2@colData)) stop("rbind: not the same @colData!")
    p@counts=rbind(p@counts, p2@counts)
    p@anno=rbind(p@anno, p2@anno)
    if (!identical(p@supp, p2@supp)) p@supp=list(p@supp, p2@supp)
  }
  return(p)
}


#' Combine multiple PACdatasets by rows
#'
#' @param ... one or more PACdataset object.
#' @return If params are PACdatasets with the same column definations, then will combine rows of PACds. See ?base::cbind for the value returned by the default methods.
#' @name rbind
#' @family PACdataset functions
#' @export
rbind <- function (...) {
  if (length(attr(list(...)[[1]], "class"))>0) {
    if (attr(list(...)[[1]], "class") == "PACdataset") return(rbind.PACds(...))
  }
  return(base::rbind(...))
}

# -------- class movRes -------------

#' movRes object
#'
#' movRes is an S4 class to represent comparison results among conditions by mov-like funcions (movDEPAC etc.).
#'
#' All other movRes classes, such as movDEGeneRes, movDEPACRes, movUTRTrendRes,
#' movAPASwitchRes, movIndexDiffRes are generic from movRes. Use movSelect, movStat,
#' and movViz to process the movRes class.
#' @slot type type of the movRes object, such as DEGENE, DEPAC.
#' @slot method define the method used for obtaining the movRes object, like DE, DEX, lineartrend, etc.
#' @slot group group label of the samples.
#' @slot conds a data frame with two columns specifying the two conditions (cond1 and cond2) in the sample group for the comparison analysis.
#' Rownames of the conds data frame are like "cond1.cond2".
#' @slot pairwise a list with two elements, padj and value, which are normally obtained from getHeatmapRes().
#' @family movRes objects
#' @seealso \code{\link{movSelect}} for select movRes, \code{\link{movStat}} for making statistics, \code{\link{movViz}} to visualize a gene.
setClass("movRes", slots=list(type="character", method="character",
                             group="character", conds="data.frame", pairwise="list"))


#' movDEGeneRes object
#'
#' An S4 class to represent DE gene results among conditions.
#'
#' The movDEGeneRes object is obtained by calling movDEGene function.
#' @family movRes objects
#' @seealso \code{\link{movDEGene}} for DE gene analysis on a PACdataset.
setClass("movDEGeneRes", slots=list(), contains="movRes", prototype = list(type="DEGENE", method='DESEQ2'))

#' movDEPACRes object
#'
#' movUTRTrendRes is an S4 class to represent DE PAC (poly(A) site cluster) results among conditions.
#'
#' The movDEPACRes object is obtained by calling movDEPAC function.
#' @slot PACusage When method=DEXseq, PACusage slot stores the exon (PAC) usage from DEXseq.
#' @family movRes objects
#' \code{\link{movDEPAC}} for DE PAC analysis on a PACdataset.
setClass("movDEPACRes", slots=list(PACusage="data.frame"), contains="movRes", prototype = list(type="DEPAC"))

#' movUTRTrendRes object
#'
#' movUTRTrendRes is an S4 class to represent 3' UTR switching results among conditions.
#'
#' The movUTRTrendRes object is obtained by calling movUTRtrend function.
#'
#' @slot fullList records the full information of 3' UTR switching results,
#' including 3' UTR length, poly(A) sites expression (PAs1 and PAs2), etc.
#' @family movRes objects
#' @seealso \code{\link{movUTRtrend}} for 3' UTR switching analysis on a PACdataset.
setClass("movUTRTrendRes", slots=list(fullList="list"), contains="movRes", prototype = list(type="UTRtrend"))

#' movAPASwitchRes object
#'
#' movAPASwitchRes is an S4 class to represent APA site switching results among conditions.
#'
#' The movAPASwitchRes object is obtained by calling movAPAswitch function.
#' movAPASwitchRes is for two PACs; movUTRTrendRes is for two or more PACs (dependent on the method used); movDEPACRes is only for one PAC.
#' @slot fullList records the full information of APA site switching results,
#' including 3' UTR length, poly(A) sites expression (PA1 and PA2), etc.
#' @family movRes objects
#' @seealso \code{\link{movAPAswitch}} for APA site switching analysis on a PACdataset.
setClass("movAPASwitchRes", slots=list(fullList="list"), contains="movRes", prototype = list(type="APAswitch"))


#' movIndexDiffRes object
#'
#' movIndexDiffRes is an S4 class to represent results of APA index difference among conditions.
#'
#' The movIndexDiffRes object is obtained by calling movAPAindexDiff function.
#'
#' @slot fullList records the full information of APA index difference results.
#' @family movRes objects
#' @seealso \code{\link{movAPAindexDiff}} for APA index difference analysis on a PACdataset.
setClass("movIndexDiffRes", slots=list(fullList="list"), contains="movRes", prototype = list(type="IndexDiff"))


setMethod("show",
          "movRes",
          function(object) {
            cat('class',class(object),'\n')
            cat('@type',object@type,'\n')
            cat('@method',object@method,'\n')
            cat('@group',object@group,'\n')
            cat('condpair#',nrow(object@conds),'\n')
            cat('nrow# in pairwise',nrow(object@pairwise$padj),'\n')
            cat('@conds\n')
            print(object@conds)
          }
)

# ---- movSelect -----


#' Select movRes
#'
#' movSelect is a generic function to subset results from a pair of conditions of a movRes object.
#'
#' movSelect(aMovRes) generic function is used for subsetting results from a movRes object.
#' This generic function is only internally used by other movSelect functions.
#' If the order of condpair is not the same as the order of condpair in aMovRes, then the value of aMovRes is multiplied by -1.
#' @param aMovRes a movRes object.
#' @param condpair a pair of conditions, can be c('cond1','cond2') or 'cond1.cond2' or 'cond2.cond1'.
#' @param padjThd cutoff to filter rows with aMovRes@pairwise$padj<padjThd.
#' @param valueThd cutoff to filter rows with |aMovRes@pairwise$value|>=valueThd.
#' @param upThd cutoff to filter up-regulated rows with aMovRes@pairwise$value>=upThd
#' @param dnThd cutoff to filter down-regulated rows with aMovRes@pairwise$value<=dnThd.
#' Only one parameter from valueThd, upThd, dnThd can be specified.
#' @param ... dependent on specific movSelect funtions.
#' @return A list.
#' Output of other movSelect functions is dependent on the specific movRes class and the out parameter.
#' @family mov-like functions
#' @export
setGeneric("movSelect", function(aMovRes, condpair, padjThd=NULL, valueThd=NULL, upThd=NULL, dnThd=NULL, ...) standardGeneric("movSelect"))

# Filter row index by amovRes$pairwise@padj and value
# If movDEPACRes is a DEXSEQ, the additional processing.
# return: list(rid=rid, heat=heat(already flipped), condname=condname)
.filterMovResRowID<-function(aMovRes, condpair, padjThd=NULL,
                             valueThd=NULL, upThd=NULL, dnThd=NULL) {

  if ((!is.null(valueThd) | !is.null(upThd)) & (!is.null(upThd) | !is.null(dnThd))  & (!is.null(valueThd) | !is.null(dnThd)) ) stop("valueThd, upThd, dnThd, only one shoud be not NULL")

  condname=getMovResPairwiseCondName(aMovRes, condpair,  check=TRUE, markFlip=TRUE)
  flip=condname$flip
  condname=condname$condname

  #if (is.null(valueThd) & is.null(padjThd) ) stop("valueThd, padjThd shoud not both be NULL")

  #if condname is flipped, then value is *-1
  flip=ifelse(flip, -1, 1)
  heat=data.frame(padj=aMovRes@pairwise$padj[, condname], value=aMovRes@pairwise$value[,condname]*flip)
  if (flip==-1) cat("Warning: condpair is flip of movRes@conds, so movRes@pairwise$value*(-1)\n")

  #If is DEXPAC, then record the row id of the max value in current condpair
  #if DEXPAC, then filter results also by maxLogFC per row
  #idx=which(fc==maxfc | fc>=logFC)
  #if (is.null(imax)) {
  imax=cbind(rmax=rep(F, nrow(aMovRes@pairwise$value)), rmin=rep(F, nrow(aMovRes@pairwise$value)), amax=rep(F, nrow(aMovRes@pairwise$value)))
  #}

  if (inherits(aMovRes, 'movDEPACRes') & toupper(aMovRes@method)=='DEXSEQ') {
    aMovRes@pairwise$value=aMovRes@pairwise$value*flip
    amax=apply(abs(aMovRes@pairwise$value), 1, max)
    maxfc=amax
    amax=abs(aMovRes@pairwise$value[,condname])==amax #all
    rmax=amax & aMovRes@pairwise$value[,condname]>0  #UP
    rmin=amax & aMovRes@pairwise$value[,condname]<0 #DOWN
    imax=cbind(rmax,rmin,amax)
    #If no valueThd is given for DEXPAC, the automatelly set one FC threshold
    #otherwise, the filtered DEPAC may be actually DE from other condpairs but not the given pair
    if (is.null(valueThd) & is.null(upThd) & is.null(dnThd)) {
      valueThd=min(maxfc)
      cat("Warning: movRes is DEXPAC, but valueThd/upThd/dnThd are all NULL, mannually set valueThd=min(maxfc)=", valueThd,'\n')
    }
  }

  imax=.asDf(imax)

  heat$value[is.na(heat$value)]=0

  if (length(heat$padj)==0) {
    return(list(rid=c()))
  }

  tf=rep(TRUE, nrow(heat))
  if (!is.null(padjThd)) {
    tf=!is.na(heat$padj)
    tf=tf & heat$padj<padjThd
    imax$rmax=tf & imax$rmax  #record row ids meet both imax and padj, these imax rows would surely be DEXPAC
    imax$rmin=tf & imax$rmin
    imax$amax=tf & imax$amax
  }
  if (!is.null(upThd)) {
    tf=tf & !is.na(heat$value) & heat$value>=upThd
    tf=tf | imax$rmax
    if (sum(imax$rmax)>0) cat("Warning: movRes is DEXPAC, also filter by rowMax(movRes@pairwise$value)\n")
  }
  if (!is.null(dnThd)) {
    tf=tf & !is.na(heat$value) & heat$value<=dnThd
    tf=tf | imax$rmin
    if (sum(imax$rmin)>0) cat("Warning: movRes is DEXPAC, also filter by rowMax(movRes@pairwise$value)\n")
  }
  if (!is.null(valueThd)) {
    tf=tf & !is.na(heat$value) & abs(heat$value)>=valueThd
    tf=tf | imax$amax
    if (sum(imax$amax)>0) cat("Warning: movRes is DEXPAC, also filter by rowMax(movRes@pairwise$value)\n")
  }

  rid=which(tf)
  return(list(rid=rid, heat=heat, condname=condname))
}


#' Select movRes
#'
#' movSelect(movDEGeneRes) subset results from a movDEGeneRes object.
#'
#' @param out should be [gene/pv/full] for movDEGeneRes, or [pa/pv/full] for movDEPACRes,
#' or [gene/pa/pv/full] for movUTRTrendRes and movAPASwitchRes, ignoring case.
#' out=gene, to output a gene list.
#' out=pv, to output a data frame with two columns, padj and value, with each row being a gene.
#' out=full, to output a data frame with full information. The content depends on the movRes class:
#' a data frame with detailed [counts + padj + value] for movDEGeneRes;
#'  a data frame of movUTRTrendRes@fullList for movUTRTrendRes;  a data frame of [PACds@anno + counts + value + padj] for movDEPACRes.
#' Note: If the given condpair is not in the same order as the movRes, for example, condpair=cond1.cond2, but in movRes it is cond2.cond1,
#' then the output value is movRes$value*(-1).
#' @param PACds a PACdataset object which is used for calculating movRes. If out=full, then need to specify PACds.
#' @examples
#' ## Examples for movSelect(movDEGeneRes) ##
#' data(PACds)
#' ## Detect DE genes in all pairwise conditions.
#' DEgene=movDEGene(PACds=PACds, method='DESeq2', group='group', minSumPAT=50)
#'
#' ## Select DE gene results with full information including the read counts in each sample.
#' selFull=movSelect(DEgene, condpair='embryo.anther',
#'               padjThd=0.05, valueThd=1, out='full', PACds=PACds)
#' head(selFull)
#'
#' ##  Select DE gene results with only padj and value. Here value is log2(embryo/anther).
#' sel1=movSelect(DEgene, condpair='embryo.anther', padjThd=0.05, valueThd=1, out='pv')
#' head(sel1)
#'
#' ## Here the condpair is Y.X, so the value is log2(anther/embryo)
#' sel2=movSelect(DEgene, condpair='anther.embryo', padjThd=0.05, valueThd=1, out='pv')
#' head(sel2)
#'
#' ## Output gene names of DE genes.
#' sel=movSelect(DEgene, condpair='embryo.anther', padjThd=0.05, upThd=0.5, out='gene')
#' head(sel)

#' @describeIn movSelect for a movDEGeneRes object.
#' @export
setMethod("movSelect", signature(aMovRes = "movDEGeneRes"), def=function(aMovRes, condpair, padjThd=NULL,
                                                                         valueThd=NULL, upThd=NULL, dnThd=NULL, out, PACds=NULL) {
  out=toupper(out)
  if (is.null(PACds) & out=='FULL') stop("movSelect: out=full but PACds=NULL\n")
  if (!(out %in% c('GENE','PV','FULL'))) stop("movSelect: error aMovRes class and out\n")

  filtered=.filterMovResRowID(aMovRes, condpair=condpair, padjThd=padjThd, valueThd=valueThd, upThd=upThd, dnThd=dnThd)
  rid=filtered$rid
  condname=filtered$condname
  heat=filtered$heat

  if (out=='FULL') { #Output [DEGENE=gene]+[DEPAC=anno]+[DEPAC/DEGENE=counts]+value, padj
    dds=PACds2geneDs(PACds)
    #allres=matrix[sample_counts, .log2FC, .padj]
    allres=heat
    allres$gene=rownames(aMovRes@pairwise$padj)
    allres=allres[rid, , drop=F]

    counts=SummarizedExperiment::assays(dds)$counts

    counts=cbind(gene=rownames(counts), counts)
    allres=.asDf(merge(counts,allres,by.x='gene',by.y='gene'))
    #allres=allres[order(allres$padj), ]
    return(allres)

  } #~full

  if (out=='PV') {
    pv=.asDf(heat[rid, , drop=F])
    rownames(pv)=rownames(aMovRes@pairwise$padj)[rid]
    #pv=pv[order(pv$padj), ]
    return(pv)
  }

  if (out=='GENE') {
    return(rownames(aMovRes@pairwise$padj)[rid])
  }
}
)


#' Select movRes
#'
#' movSelect(movDEPACRes) subset results from a movDEPACRes object.
#'
#' @examples
#' ## Examples for movSelect(movDEPACRes) ##
#' data(PACds)
#'
#' ## Detect DE PACs in all pairwise conditions
#' ## using DESeq2 method, only PACs with total read counts in all samples >=20 are used.
#' DEPAC=movDEPAC(PACds, method='DESeq2', group='group', minSumPAT=20)

#' ## Get PA list
#' movSelect(aMovRes=DEPAC, condpair='embryo.anther', padjThd=0.1, valueThd=0.5, out='PA')

#' ## Get data frame of [padj, value]
#' movSelect(aMovRes=DEPAC, condpair='embryo.anther', padjThd=0.1, valueThd=0.5, out='pv')

#' ## Detect DE PACs in all pairwise conditions using DEXseq method.
#' DEXPAC=movDEPAC(PACds, method='DEXseq', group='group', minSumPAT=20)
#' ## Output full PA list
#' full=movSelect(aMovRes=DEXPAC, condpair='anther.embryo', padjThd=0.1, out='full', PACds=PACds)

#' ## For method=DEXPAC, note that even valueThid=2,
#' ## the final list may have items of value<2, because these items meet the padj cutoff.
#' pv=movSelect(aMovRes=DEXPAC, condpair='embryo.anther',
#'              padjThd=0.1, valueThd=2, out='pv'); min(abs(pv$value))

#' ## When condpair is not the same order as in movRes, the output value is multiplied by -1
#' pv=movSelect(aMovRes=DEXPAC, condpair='anther.embryo',
#'              padjThd=0.1, valueThd=2, out='pv'); min(abs(pv$value))

#' ##Get list for only upregulation
#' movSelect(aMovRes=DEXPAC, condpair='anther.embryo', padjThd=0.1, upThd=2, out='pv')

#' ## If method is DEXPAC , when not set value thresholds,
#' ##will automatelly use min(maxfc) to filter DE PACs.
#' movSelect(aMovRes=DEXPAC, condpair='anther.embryo', padjThd=0.1, out='pv')
#' @describeIn movSelect for a movDEPACRes object.
#' @export
setMethod("movSelect", signature(aMovRes = "movDEPACRes"), def=function(aMovRes, condpair, padjThd=NULL, valueThd=NULL,
                                                                        upThd=NULL, dnThd=NULL, out, PACds=NULL) {

  if (!is.null(PACds)) PACds@counts=asAnyMatrix(PACds@counts)

  if (toupper(aMovRes@method)=='CHISQ') {
    if (!is.null(upThd) | !is.null(dnThd)) {
      stop("DE PACs are obtained from the Chisq method, the value matrix of aMovRes is the chisq.test's pvalue of the whole gene, cannot set upThd or dnThd!\n")
    }
  }

  out=toupper(out)
  if (is.null(PACds) & out=='FULL') stop("movSelect: out=full but PACds=NULL\n")

  if (!(out %in% c('PA','PV','FULL'))) stop("movSelect: error aMovRes class and out\n")

  filtered=.filterMovResRowID(aMovRes, condpair=condpair, padjThd=padjThd, valueThd=valueThd, upThd=upThd, dnThd=dnThd)
  rid=filtered$rid
  condname=filtered$condname
  heat=filtered$heat

  if (out=='FULL') { #output: [DEGENE=gene]+[DEPAC=anno]+[DEPAC/DEGENE=counts]+value, padj
    dds=NULL
    if (toupper(aMovRes@method)=='DESEQ2') {
      dds=PACds2DESeqDds(PACds, noIntergenic = TRUE)
    } else if (toupper(aMovRes@method)=='DEXSEQ') {
      dds=PACdataset2Dxd(PACds)
    }

    #allres=matrix[sample_counts, .log2FC, .padj] padj, logFC, and counts from all cond pairs.
    allres=heat
    allres$gene=rownames(aMovRes@pairwise$padj)
    allres=allres[rid, , drop=F]

    if (!is.null(dds)) { #method=deseq/dexseq
      counts=SummarizedExperiment::assays(dds)$counts

      if (toupper(aMovRes@method)=='DEXSEQ') {
        counts=counts[,which(dds@colData$exon=='this')]
        colnames(counts)=dds@colData$sample[which(dds@colData$exon=='this')]
      }
    } else { #method=chisq
      counts=PACds@counts
    }

    counts=cbind(gene=rownames(counts), .asDf(counts))
    allres=merge(counts, allres, by.x='gene',by.y='gene')

    colnames(allres)[colnames(allres)=='gene']='PA'
    if (!is.null(PACds@anno)) {
      if (length(base::intersect(rownames(PACds@anno),allres$PA))==0) { #If PA rows are differently named.
        PA=paste0(PACds@anno$gene,':',PACds@anno$coord)
      } else {
        PA=rownames(PACds@anno)
      }
      anno=cbind(PA=PA,PACds@anno)
      allres=merge(anno,allres, by.x='PA',by.y='PA')
    }

    #allres=allres[order(allres$padj), ]
    return(allres)
  } #~full

  if (out=='PV') {
    pv=.asDf(heat[rid, , drop=F])
    rownames(pv)=rownames(aMovRes@pairwise$padj)[rid]
    #pv=pv[order(pv$padj), ]
    return(pv)
  }

  if (out=='PA') {
    return(rownames(aMovRes@pairwise$padj)[rid])
  }

}
)

# Get PAs from UTR/swRes object for given genes
# @param aMovRes A movUTRTrendRes/movAPASwitchRes
# @param genes A gene vector
# @param condname cond1.cond2 or c(cond1, cond2) or c(cond2, cond1)
#        If it is not in aMovRes, then return NULL.
# @return: PA names
.getPAbyGeneFromMovUTRTrendRes<-function(aMovRes, genes, condname) {
  if (!(class(aMovRes) %in% c('movUTRTrendRes','movAPASwitchRes'))) stop("movRes class!=movUTRTrendRes/movAPASwitchRes")
  if (length(condname)==2) {
    condname=paste0(condname[1],'.',condname[2])
  }

  if (!(condname %in% names(aMovRes@fullList))) {
    condname=unlist(strsplit(condname, split='\\.'))
    condname=.getMatchCondPairs(condname, aMovRes@conds, markFlip = FALSE)
    if (is.null(condname)) {
      cat('condname not in aMovRes\n')
      return(NULL)
    }
    condname=paste0(condname[1,1],'.',condname[1,2])
  }

  sw=aMovRes@fullList[[condname]]
  #print(head(sw))
  idx=which(sw$gene %in% genes)
  cols=colnames(sw)
  if (length(idx)==0) return(NULL)
  if ('PA1' %in% cols & 'PA2' %in% cols) {
    PAs=unique(c(sw$PA1[idx],sw$PA2[idx]))
  } else if ('PAs1' %in% cols & 'PAs1' %in% cols) {
    PAs=unique(c(sw$PAs1[idx],sw$PAs2[idx])) #"PH02Gene00202:PA3=92;PH02Gene00202:PA4=42"
    PAs=unlist(strsplit(PAs,';'))
    PAs=gsub('=.*','',PAs)
  } else {
    stop("movRes class=movUTRTrendRes/movAPASwitchRes, out=PA, but PA1/2 or PAs1/2 not in fullList cols")
  }
  return(unique(PAs))
}

#' Select movRes
#'
#' movSelect(movUTRTrendRes) subset results from a movUTRTrendRes object.
#'
#' @describeIn movSelect for a movUTRTrendRes object.
#' @export
setMethod("movSelect", signature(aMovRes = "movUTRTrendRes"), def=function(aMovRes, condpair, padjThd=NULL,
                                                                           valueThd=NULL, upThd=NULL, dnThd=NULL, out) {

  out=toupper(out)
  if (!(out %in% c('GENE','PA','PV','FULL'))) stop("movSelect: error aMovRes class and out\n")

  filtered=.filterMovResRowID(aMovRes, condpair=condpair, padjThd=padjThd, valueThd=valueThd, upThd=upThd, dnThd=dnThd)
  rid=filtered$rid
  if (is.null(rid)) return(NULL)

  condname=filtered$condname
  heat=filtered$heat

  if (out=='PV') {
    pv=.asDf(heat[rid, , drop=F])
    rownames(pv)=rownames(aMovRes@pairwise$padj)[rid]
    return(pv)
  }

  if (out=='PA') {
    genes=rownames(aMovRes@pairwise$padj)[rid]
    return(.getPAbyGeneFromMovUTRTrendRes(aMovRes, genes, condname))
  }

  if (out=='GENE') {
    return(rownames(aMovRes@pairwise$padj)[rid])
  }

  if (out=='FULL') {
    genes=rownames(aMovRes@pairwise$padj)[rid]
    sw=aMovRes@fullList[[condname]]
    return(sw[sw$gene %in% genes,])
  }
}
)


#' Select movRes
#'
#' movSelect(movAPASwitchRes) subset results from a movAPASwitchRes object.
#'
#' @describeIn movSelect for a movAPASwitchRes object.
#' @export
setMethod("movSelect", signature(aMovRes = "movAPASwitchRes"), def=function(aMovRes, condpair, padjThd=NULL,
                                                                            valueThd=NULL, upThd=NULL, dnThd=NULL, out) {
  aMovRes=new("movUTRTrendRes", group=aMovRes@group, method=aMovRes@method, conds=aMovRes@conds, pairwise=aMovRes@pairwise, fullList=aMovRes@fullList)
  return(movSelect(aMovRes=aMovRes, condpair=condpair, padjThd=padjThd,valueThd=valueThd, upThd=upThd, dnThd=dnThd, out=out))
}
)


#' Select movRes
#'
#' movSelect(movIndexDiffRes) subset results from a movIndexDiffRes object.
#'
#' @describeIn movSelect for a movIndexDiffRes object.
#' @export
setMethod("movSelect", signature(aMovRes = "movIndexDiffRes"), def=function(aMovRes, condpair, padjThd=NULL,
                                                                            valueThd=NULL, upThd=NULL, dnThd=NULL, out) {
  out=toupper(out)
  if (!(out %in% c('GENE','PV','FULL'))) stop("movSelect: error aMovRes class and out\n")
  if (out=='FULL' & length(aMovRes@fullList)==0) { return(NULL) }
  aMovRes=new("movUTRTrendRes", group=aMovRes@group, method=aMovRes@method, conds=aMovRes@conds, pairwise=aMovRes@pairwise, fullList=aMovRes@fullList)
  return(movSelect(aMovRes=aMovRes, condpair=condpair, padjThd=padjThd,valueThd=valueThd, upThd=upThd, dnThd=dnThd, out=out))
}
)

# Get condpair row names from movRes, normally used to get condname from movRes$pairwise[, condpair]
# Get unified condname from movRes@conds
# @param aMovRes A movRes class
# @param condpair Can be cond1.cond2 or c(cond1, cond2) (a pair of conds)
# @param check If TRUE, then if not found in movRes, then stop
# @param markFlip If T, then return list(condname=A.B, flip=TRUE/FALSE); else return only A.B
# @return rownames(aMovRes$conds)[i]
# @examples
# getMovResPairwiseCondName(DEgene, c('0350','0140')); head(DEgene@pairwise$padj)
# getMovResPairwiseCondName(DEgene, 'embryo.anther'); head(DEgene@pairwise$padj)
# getMovResPairwiseCondName(DEgene, 'a.b', check=F)
# getMovResPairwiseCondName(DEgene, 'embryo.anther', markFlip=FALSE);
getMovResPairwiseCondName <- function (aMovRes, condpair, check=TRUE, markFlip=FALSE) {
  if (length(condpair==1)) {
    condpair=unlist(strsplit(x=condpair, split='.', fixed=TRUE))
    if (length(condpair)!=2) {
      stop(cat("getMovResPairwiseCondName: error condpair", condpair, '(not A.B or c(A,B)\n'))
    }
  }
  x=apply(aMovRes@conds, 1, function(par) return(length(base::intersect(par, condpair))))
  idx=which(x==2)
  if (length(idx)>1) stop(cat("getMovResPairwiseCondName: more than one condpair from aMovRes for",condpair,'\n'))
  if (length(idx)==0) {
    if (check) {
      stop(cat("getMovResPairwiseCondName: NO condpair from aMovRes for",condpair,'\n'))
    } else {
      return(NULL)
    }
  }
  if (markFlip) {
    return(list(condname=rownames(aMovRes@conds)[idx], flip=!(aMovRes@conds[idx,1]==condpair[1] & aMovRes@conds[idx,2]==condpair[2])))
  }
  return(rownames(aMovRes@conds)[idx])
}

# ---- movStat -----

#' heatmapResult object
#'
#' heatmapResult is an S4 class to represent a list with padj, value and colData. Normally it is from getHeatmapRes().
#'
#' The slot of padj or value is a data frame with each row being gene/PA, each column being a condpair or cond (if DEXSeq, then it is the exonUsage matrix).
#'
#' @slot padj a data frame recording the pvalue of each gene or PA in each condition.
#' @slot value a data frame with the same dimension as padj.
#' @slot colData stores padj/value matrix's annotation info. Its row names are the column names of the padj matrix.
#' Normally, there are two columns (cond1 and cond2) in colData.
#' @family movRes objects
heatmapResults <- setClass("heatmapResults", slots=list(padj='data.frame',value='data.frame', colData='data.frame'))

#' Convert movRes object to heatmapResults object
#'
#' movRes2heatmapResults converts movRes object to heatmapResults object,
#' so that then can use subsetHeatmap functions. This function is mainly used by statGeneralHeat().
#'
#' @param aMovRes a movRes object
#' @return A heatmapResults object.
#' @export
setGeneric("movRes2heatmapResults", function(aMovRes) standardGeneric("movRes2heatmapResults"))

#' Convert movRes object to heatmapResults object
#'
#' @describeIn movRes2heatmapResults for a movRes object
setMethod("movRes2heatmapResults", signature="movRes", def=function(aMovRes) {
  heat=new("heatmapResults", padj=aMovRes@pairwise$padj, value=aMovRes@pairwise$value, colData=aMovRes@conds)
  return(heat)
}
)

#' Make statistics for a object
#'
#' movStat(object) makes statistics for a movRes or PACdataset object.
#' It is a generic function to make statistics of a object which is normally of movRes or PACddataset class.
#'
#' @param object Normally a movRes or PACdataset object.
#' @param ... Arguments passing to other functions.
#' @family mov-like functions
#' @export
setGeneric("movStat", function(object, ...) standardGeneric("movStat"))

# -------- movStat movRes -------------

#' Make statistics for a movRes object
#'
#' movStat(movRes) makes statistics for a movRes object.
#'
#' Statistical results of movStat(movRes) are numbers of significant DE genes/PACs/switching events,
#' the significant gene/PAC list for each condition pair, etc.
#'
#' @param padjThd cutoff to filter rows with aMovRes@pairwise$padj<padjThd.
#' @param valueThd cutoff to filter rows with aMovRes@pairwise$value|>=valueThd.
#' @param upThd cutoff to filter up-regulated rows with aMovRes@pairwise$value>=upThd.
#' @param dnThd cutoff to filter down-regulated rows with aMovRes@pairwise$value<=dnThd.
#' @return A list of several data frames including nsig, siglist, tf01, ovp, de01, and deNum.
#' @describeIn movStat for a movRes object.
#' @export
setMethod("movStat", signature="movRes", def=function(object, padjThd=NULL,
                                                      valueThd=NULL, upThd=NULL, dnThd=NULL) {
  heat=movRes2heatmapResults(object)
  ol=statGeneralHeat(heat=heat, padjThd=padjThd, valueThd=valueThd, upThd=upThd, dnThd=dnThd, lbl='sig.num')
  return(ol)
}
)

#' Make statistics for a PACdataset object
#'
#' movStat(PACdataset) makes statistics for a PACdataset object.
#'
#' Statistical results of movStat(PACdataset) are the distributions PATs, PACs etc.
#'
#' @param minPAT min PAT number to filter valid PACs with expression level >= minPAT, can be like 5 or c(5,10)
#' @param ofilePrefix output stat results to a txt file named <ofilePrefix.patX.stat>. If it is NULL, then do not output to file.
#' @return A list of data frames representing results from different cutoffs like [PAT1, PAT5, PAT10]. Each item is a data frame recording statistical results.
#' @seealso \code{\link{plotPACdsStat}} for plot the stat results of PACdataset.
#' @examples
#' ## movStat(PACdataset) ##
#' data(PACds)
#' ## Stat a PACds, filtering PACs with cutoffs of 1, 5, and 10.
#' pstats=movStat(PACds, minPAT=c(1,5,10), ofilePrefix='statFileEach')
#'
#' ## If the PACds includes replicates, you may want to pool replicates first and then make stats.
#' pstats=movStat(subsetPACds(PACds, group='group', pool=TRUE),
#'               minPAT=c(1,5,10), ofilePrefix='statFilePool')
#' @describeIn movStat for a PACdataset object.
#' @export
setMethod("movStat", signature="PACdataset", def=function(object, minPAT=c(1,5,10), ofilePrefix=NULL) {

  ## bug: not suitable for large table

  ostats=list()
  if (!isPACdsAnnotated(object)) {
    cat("The PACdataset object is not annotated, please call annotatePACds to annotate first!\n")
    return(NULL)
  }
  for (patThd in minPAT) {

    count=.asDf(object@counts)
    count[count<patThd]=0
    count$total=rowSums(count)

    nPAC=colSums(count>0)
    nPAT=colSums(count)

    ostat=data.frame(nPAC, nPAT)
    rownames(ostat)=colnames(count)

    #stat non-igt
    idx=getNonItgFtrId(object@anno$ftr)
    count2=count
    count2$gene=object@anno$gene
    count2=count2[idx,]
    nGene=c(); nAPAgene=c(); APAextent=c()
    for (i in 1:(ncol(count2)-1)) {
      s=count2[,c(i,ncol(count2))]
      s=s[s[,1]>0, ]
      nGene=c(nGene, length(unique(s$gene)))
      nAPAgene=c(nAPAgene, length(unique(s$gene[duplicated(s$gene)])))
    }
    ostat=cbind(ostat, nGene, nAPAgene)
    ostat$APAextent=ostat$nAPAgene/ostat$nGene

    #region distribution
    count$ftr=object@anno$ftr
    count$ftr=as.factor(count$ftr)
    nPAT=aggregate(. ~ ftr, data = count, sum)
    rownames(nPAT)=nPAT$ftr
    nPAT$ftr=NULL
    nPAT=t(nPAT)
    colnames(nPAT)=paste0(colnames(nPAT),'_','nPAT')
    count[,-which(colnames(count)=='ftr')]=count[,-which(colnames(count)=='ftr')]>0
    nPAC=aggregate(. ~ ftr, data = count, sum)
    rownames(nPAC)=nPAC$ftr
    nPAC$ftr=NULL
    nPAC=t(nPAC)
    colnames(nPAC)=paste0(colnames(nPAC),'_','nPAC')
    ostat=cbind(ostat, nPAT, nPAC)

    ostats[[paste0('pat',patThd)]]=ostat

    if (!is.null(ofilePrefix)) {
      ofile=paste0(ofilePrefix,'.pat',patThd,'.stat')
      ostat=cbind(cond=rownames(ostat),ostat)
      write.table(ostat, file=ofile, row.names = FALSE, col.names = TRUE, sep="\t", quote=F)
      cat('>>>',ofile,'\n')
    }

  }

  return(ostats)
}
)



stat_plot_theme<- ggplot2::theme_bw() + ggplot2::theme(
  #text=element_text(family="Arial"),
  axis.title.x =ggplot2::element_text(color="black",size=16) ,
  axis.title.y =ggplot2::element_text(color="black",size=16) ,
  axis.text.x =ggplot2::element_text(color="black",size=14) ,
  axis.text.y =ggplot2::element_text(color="black",size=14) ,
  legend.text =ggplot2::element_text(color="black",size=14),
  legend.title=ggplot2::element_text(color="black",size=14),
  legend.background = ggplot2::element_blank(),
  #panel.border = element_blank(),
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  strip.text = ggplot2::element_text(size=14)
  )

## gene model without background
viz_geneModel_theme<- ggplot2::theme_bw() + ggplot2::theme(
  panel.border=ggplot2::element_rect(colour = "grey", fill=NA),
  legend.position = "none",
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank()
)


#' Plot movStat results
#'
#' plotPACdsStat plots the results from movStat(PACdataset).
#'
#' Output plots are barplots with facet of minPAT.
#'
#' @param pstats a list from movStat(PACdataset).
#' Each item in the list is a data frame, including columns [nPAC nPAT nGene nAPAgene APAextent 3UTR_nPAT 5UTR_nPAT
#' AMB_nPAT CDS_nPAT intron_nPAT 3UTR_nPAC 5UTR_nPAC AMB_nPAC CDS_nPAC intron_nPAC]
#' @param pdfFile If not NULL, then print to pdfFile.
#' @param minPAT min number of PATs to filter PACs with expression levels >= minPAT, values like NULL/5/c(1,5,10)
#' @param conds conditions to make statistics, values like NULL / c('0400','1200') / 'total'.
#' All conds will be in one plot.
#' @param ... Dot arguments passing to other functions.
#' @return NULL. Output to text file and/or print the plots to a device.
#' @examples
#'
#' data(PACds)
#' ## To make statistics of distributions of PACs for each sample, first we pooled replicates.
#' PACds1=subsetPACds(PACds, group='group', pool=TRUE)
#' head(PACds@counts)
#' head(PACds1@counts)
#'
#' ## Make statistics of distribution of PACs using different PAT cutoffs.
#' ## minPAT=5 means that only PACs with >=5 reads are used for statistics.
#' pstats=movStat(PACds1, minPAT=c(1, 5, 10, 20, 50, 60), ofilePrefix=NULL)
#' names(pstats)
#' pstats$pat10
#'
#' ## Plot statistical results by barplots: PAC#, PAT#, APA gene%, PAC%,
#' ## PAT% across samples and genomic regions.
#' ## Plot all statistical results with all cutoffs,
#' ## with each plot having six smaller plots corresponding to the six cutoffs.
#' plotPACdsStat(pstats, pdfFile='PACds_stat.pdf', minPAT=c(5,10))
#' ## Plot specific cutoffs and conditions.
#' plotPACdsStat(pstats, pdfFile='PACds_stat_anther_embryo.pdf',
#'              minPAT=c(5,10), conds=c('anther','embryo'))
#' ## Plot the overall distributions using pooled samples (total) and two cutoffs.
#' plotPACdsStat(pstats, pdfFile='PACds_stat_total.pdf',
#'              minPAT=c(5,10), conds=c('total'))
#' ## Plot the overall distributions using pooled samples (total) and one cutoff.
#' plotPACdsStat(pstats, pdfFile='PACds_stat_total_PAT10.pdf',
#'              minPAT=c(10), conds=c('total'))
#'
#' ## Plot to the current device, but just the last plot will be shown.
#' plotPACdsStat(pstats, pdfFile=NULL, minPAT=c(5,10))
#' @importFrom reshape2 melt
#' @family PACdataset functions
#' @seealso \code{\link{movStat}} to get statistical results.
#' @export
plotPACdsStat <- function(pstats, pdfFile=NULL, minPAT=NULL, conds=NULL, ...) {

  #library(ggplot2, verbose = FALSE)

  #Add minPAT column
  patThd=as.numeric(gsub('pat','',names(pstats)))
  if (is.null(minPAT)) minPAT=patThd
  if (!AinB(minPAT, patThd)) stop("minPAT not all in pstats")

  if (length(minPAT)>6) {
    cat("Warning: the plot may be messy, please use [minPAT] to choose <=6 cutoffs!\n")
  }


  if (!is.null(pdfFile))  {
    if (length(minPAT)<=2) {
      height=8
    } else {
      height=floor((length(minPAT)+1)/2)*4
    }
    pdf(file=pdfFile, width=10, height=height, ...)
    cat('Plot >>>> ', pdfFile,'\n')
  }


  pstat=pstats[[paste0('pat',minPAT[1])]]
  pstat$minPAT=minPAT[1]
  pstat$cond=rownames(pstat)

  if (nrow(pstat)>=5) {
    cat("Warning: the plot may be messy, please use [conds] to choose <=5 conditions!\n")
  }


  for (t in minPAT[-1]) {
    x=pstats[[paste0('pat',t)]]
    x$minPAT=t
    x$cond=rownames(x)
    pstat=rbind(pstat, x)
  }

  #subset conds
  if (is.null(conds)) {
    conds=unique(pstat$cond)
    conds=conds[conds!='total']
  }
  if (!AinB(conds, unique(pstat$cond))) stop("conds not all in pstat")
  pstat=pstat[pstat$cond %in% conds, , drop=F]

  #order by minPAT as factors, used for facet
  pstat$minPAT=pstat$minPAT[order(pstat$minPAT)]
  pstat$minPAT=factor(pstat$minPAT, levels=unique(pstat$minPAT), labels=unique(c(paste0('PAT>=',pstat$minPAT))), ordered=TRUE)

  #nPAC
  d=pstat[, c('cond','minPAT','nPAC')]
  p=ggplot(d, aes(x = cond, y=nPAC, fill=cond)) + ggplot2::geom_bar(stat="identity") +
    geom_text(aes(label=nPAC), vjust=-0.2,col="black",size=4) +
    ylab('Number of PACs')+ xlab('Condition') + ggtitle("Number of PACs") + guides(fill=FALSE) +
    facet_wrap( ~ minPAT, ncol=2) + stat_plot_theme
  print(p)
  cat('Plot Number of PACs\n')

  #nPAT
  d=pstat[, c('cond','minPAT','nPAT')]
  p=ggplot(d, aes(x = cond, y=nPAT, fill=cond)) + ggplot2::geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%d", nPAT)), vjust=-0.2, col="black",size=4) +
    ylab('Number of PATs')+ xlab('Condition') + ggtitle("Number of PATs")  + guides(fill=FALSE) +
    facet_wrap( ~ minPAT, ncol=2) + stat_plot_theme
  print(p)
  cat('Plot Number of PATs\n')

  #single, APA bar plot, and text APA%
  d=pstat[, c('cond','minPAT','nGene','nAPAgene','APAextent')]
  d$singlePA=d$nGene-d$nAPAgene
  d$APA=d$nAPAgene
  d$nAPAgene=NULL; d$nGene=NULL
  d=reshape2::melt(d, id.vars = c('cond','minPAT','APAextent'), value.name = 'nGene', variable.name = 'gene')
  d$APAextent=sprintf("%.0f%%", 100*d$APAextent)
  d[d$gene=='singlePA','APAextent']=NA
  p=ggplot(data=subset(d),  aes(x=cond, y=nGene, fill=gene)) +
    ggplot2::geom_bar(stat="identity")  +
    geom_text(aes(label=APAextent), vjust=-0.2, col="black",size=4) +
    ylab('Number of genes')+ xlab('Condition') + ggtitle("APA extent")  +
    facet_wrap( ~ minPAT, ncol=2) + stat_plot_theme
  print(p)
  cat('Plot APA extent\n')

  #PAC distribution
  d=pstat[, c('cond','minPAT', colnames(pstat)[grep('_nPAC', colnames(pstat))])]
  colnames(d)=gsub('_nPAC','',colnames(d))
  d=melt(d, id.vars = c('cond','minPAT'), value.name = 'nPAC', variable.name = 'region')
  p=ggplot(d, aes(x = region, y=nPAC, fill=cond)) + ggplot2::geom_bar(stat="identity", position=position_dodge()) +
    ylab('Number of PACs')+ xlab('Region') + ggtitle("PAC# distribution")  +
    stat_plot_theme +
    theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
    facet_wrap( ~ minPAT, ncol=2)
  print(p)
  cat('Plot PAC# distribution\n')

  #PAC%
  d=pstat[, c('cond','minPAT', colnames(pstat)[grep('_nPAC', colnames(pstat))])]
  d[,-c(1:2)]=d[,-c(1:2), drop=F]/rowSums(d[,-c(1:2), drop=F])
  colnames(d)=gsub('_nPAC','',colnames(d))
  d=melt(d, id.vars = c('cond','minPAT'), value.name = 'nPAC', variable.name = 'region')
  p=ggplot(d, aes(x = region, y=nPAC, fill=cond)) + ggplot2::geom_bar(stat="identity", position=position_dodge()) +
    ylab('PACs%')+ xlab('Region') + ggtitle("PAC% distribution")  +
    stat_plot_theme +
    theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
    facet_wrap( ~ minPAT, ncol=2)
  print(p)
  cat('Plot PAC% distribution\n')

  #PAT distribution
  d=pstat[, c('cond','minPAT', colnames(pstat)[grep('_nPAT', colnames(pstat))])]
  colnames(d)=gsub('_nPAT','',colnames(d))
  d=melt(d, id.vars = c('cond','minPAT'), value.name = 'nPAT', variable.name = 'region')
  p=ggplot(d, aes(x = region, y=nPAT, fill=cond)) + ggplot2::geom_bar(stat="identity", position=position_dodge()) +
    ylab('Number of PATs')+ xlab('Region') + ggtitle("PAT distribution")  +
    stat_plot_theme +
    theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
    facet_wrap( ~ minPAT, ncol=2)
  print(p)
  cat('Plot PAT# distribution\n')

  #PAT%
  d=pstat[, c('cond','minPAT', colnames(pstat)[grep('_nPAT', colnames(pstat))])]
  d[,-c(1:2)]=d[,-c(1:2), drop=F]/rowSums(d[,-c(1:2), drop=F])
  colnames(d)=gsub('_nPAT','',colnames(d))
  d=melt(d, id.vars = c('cond','minPAT'), value.name = 'nPAT', variable.name = 'region')
  p=ggplot(d, aes(x = region, y=nPAT, fill=cond)) + ggplot2::geom_bar(stat="identity", position=position_dodge()) +
    ylab('PATs%')+ xlab('Region') + ggtitle("PAT distribution (%)")  +
    stat_plot_theme +
    theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
    facet_wrap( ~ minPAT, ncol=2)
  print(p)
  cat('Plot PAT% distribution\n')

  if (!is.null(pdfFile)) invisible(dev.off())
}



# -------- class others -------------
# digDESeqAll..output object
# Each element, such as ref is from DESeq results' list[matrix=log2FC, padj, pvalue,...]
# In ref dataframe, each row is a gene, each column is an experiment, the value is log2FC/padj...
DEAllResults <- setClass("DEAllResults", slots=list(ref="list", pairwise="list", series="list"))

setMethod("show",
          "heatmapResults",
          function(object) {
            cat('@padj',nrow(object@padj),'rows\n')
            print(head(object@padj,3))
            cat('@value',nrow(object@value),'rows\n')
            print(head(object@value,3))
            cat('@colData','\n')
            print(object@colData)
          }
)

# Store switching results
swResults <- setClass("swResults",
                      slots=list(sw="list",
                                 type="character", method="character"))
#type=APAswitching, UTRtrend
#method=lineartrend, de, dex
#sw=list[[cond1.cond2]]


# ---------------- *** PACdataset *** ------------------


#' Make an example PACdataset
#'
#' makeExamplePACds returns an example PACdataset of Arabidopsis.
#'
#' @param n number of PACs.
#' @param seed random seed.
#' @param forceSparse whether to make the counts table as a sparse matrix, otherwise is a normal matrix.
#' @return A PACdataset of Arabidopsis with Chr1..Chr5. There are six columns in the anno slot.
#' @examples
#' makeExamplePACds(n=1000)
#' makeExamplePACds(n=1000, seed=456, forceSparse=T)
#' @name makeExamplePACds
#' @family PACdataset functions
#' @export
makeExamplePACds <- function(n = 1000, seed=123, forceSparse=F) {
  if (n==0) {
    PACds=new("PACdataset")
    return(PACds)
  }

  set.seed(seed)
  colData=.asDf(matrix(c(rep('control',3),rep('treat',3)),nrow=6,ncol=1,
                               dimnames =list(c(paste0('control',1:3),paste0('treat',1:3)),'group')))
  cts <- matrix(floor(rnorm(n*6, mean=50, sd=10)), nrow = n,
                              dimnames = list(paste0('PA',1:n),rownames(colData)))
  anno=sample(paste0('Chr',1:5), n, replace=TRUE)
  anno=.asDf(anno)
  colnames(anno)='chr'
  rownames(anno)=rownames(cts)
  anno$strand=sample(c('+','-'), n ,replace = TRUE)
  anno$UPA_start=floor(runif(n, min = 0, max = 1e6))
  anno$UPA_end=anno$UPA_start+floor(runif(n, min = 0, max = 100))
  anno$coord=floor((anno$UPA_start+anno$UPA_end)/2)
  anno$ftr=sample(c('3UTR','5UTR','intergenic','CDS','intron'), n, replace=TRUE)
  anno$ftr_start=anno$UPA_start-200
  anno$ftr_end=anno$UPA_end+200
  PACds=createPACdataset(counts=cts, colData=colData, anno=anno, forceSparse = forceSparse)
  return(PACds)
}

# - pacds: a PACdataset with @anno$chr, or df/matrix with chr column, or vector
# - obj: BSgenome, or FaFile, or data.frame/matrix or vector
# - col: the colname of chr in obj[if obj is df or matrix]
# - allin: pacds$chr should all in obj's chr, otherwise base::intersect is allowed
# - return: TRUE (consistent or base::intersect) or FALSE
# Ex.
# isChrConsistent(PACds, obj=bsgenome, allin=TRUE)
# isChrConsistent(PACds, obj=FaFile(aFastaFile), allin=TRUE)
# isChrConsistent(PACds, obj=unique(PACds@anno$chr), allin=TRUE)
# isChrConsistent(PACds, obj=c('x','a'), allin=TRUE)
# isChrConsistent(PACds, obj=d, allin=TRUE)
# isChrConsistent(PACds, obj=d, col='strand',allin=TRUE)
# isChrConsistent(c('1','2'), obj=c('chr1','chr2'),allin=TRUE)
isChrConsistent<-function(pacds, obj, col=NULL, allin=FALSE) {
  if (is.vector(pacds)) {
    pacchr=pacds
  } else if  (inherits(pacds, 'data.frame') | inherits(obj, 'matrix')) {
    pacchr=unique(as.character(pacds$chr))
  } else if (inherits(pacds, 'PACdataset') ) {
    pacchr=unique(pacds@anno$chr)
  }

  if (inherits(obj, 'BSgenome')) {
    objchr=GenomeInfoDb::seqnames(obj) ## 2023/5/24 fix bugs
  } else if (inherits(obj, 'FaFile')) {
    if (!file.exists(obj$index)) {
      cat("Indexing fa...\n")
      Rsamtools::indexFa(obj$path)
    }
    objchr=names(GenomeInfoDb::seqlengths(obj))
  } else if (inherits(obj, 'data.frame') | inherits(obj, 'matrix')) {
    if (is.null(col)) {
      if ('seqnames' %in% colnames(obj)) {
        col='seqnames'
      } else if ('chr' %in% colnames(obj)) {
        col='chr'
      }
    }
    if (!(col %in% colnames(obj))) {
      stop(cat('chrCol',col, 'not in obj\n'))
    }
    objchr=unique(obj[,col])
  } else {
    objchr=obj
  }
  if (allin) {
    return(AinB(pacchr, objchr, all=TRUE))
  } else {
    return(length(base::intersect(pacchr,objchr))!=0)
  }
}


# Given the ftr vectorget id that are not ^inter
# - ftr: character vector
# - return: idx of non-itg ftrs
getNonItgFtrId<-function(ftr) {
  id=grep('^inter',ftr)
  if (length(id)>0) {
    return((1:length(ftr))[-id])
  } else {
    return(1:length(ftr))
  }
}


#' Subset a PACdataset
#'
#' subsetPACds returns a subset of PACdataset.
#'
#' @param pacds A PACdataset object.
#' @param group a sample group, must be present in PACds@colData.
#' If NULL, then parameters of cond1, cond2, and avg are not in use.
#' @param cond1 a condition, must be present in PACds@colData.
#' @param cond2 a condition, must be present in PACds@colData.
#' If cond1 or cond2 is NULL, then both cond1 and cond2 are set NULL.
#' @param conds to subset multiple conditions.
#' @param avgPACtag if >0, then filter PACs with average PAT number across all samples >= avgPACtag, after subseting by group and conds.
#' @param avgGeneTag similar to avgPACtag, but for gene.
#' @param totPACtag filter PACs with total PAT number of all samples >=totPACtag.
#' @param noIntergenic whether to remove intergenic PACs.
#' @param avg TRUR/FALSE. If TRUE then take average of replicates of each of subset conditions.
#' @param pool TRUR/FALSE. If TRUE then take pool value of replicates of each of subset conditions.
#' @param choosePA value can e NULL/distal/proximal/APA. Only if PACds is all 3UTRAPA, then choosePA can be distal/proximal.
#'  APA means choose PACs from APA genes.
#' @param PAs to filter PACs by a PAC name list according to rownames in PAs.
#' @param genes to filter by a gene name list according to PACds@anno$gene.
#' @param chrs to filter by a chr name list according to PACds@anno$chr.
#' @param clearPAT If >0 then set @counts[<=clearPAT]=0 and remove blank lines and reset PACds@anno.
#' @param minExprConds default 1, min number of condtions (columns) expressed.
#' @param verbose T to show message during subsetting.
#' @return A subset of PACdataset with 0-count lines removed.
#' @examples
#' data(PACds)
#' PACds1=subsetPACds(PACds, group='group', pool=TRUE)
#' ## Subset two conditions.
#' pacds=subsetPACds(PACds, group='group',cond1='anther', cond2='embryo')
#' ## Subset PACs from APA genes in two conditions.
#' PACds1=subsetPACds(PACds, group='group', cond1='anther', cond2='embryo', choosePA='apa')
#' ## Subset PACs of given genes.
#' subsetPACds(PACds, genes=PACds@anno$gene[1:5], verbose=TRUE)
#' ## Subset PACs of given PA names.
#' swPAC=subsetPACds(PACds, PAs=rownames(PACds@counts)[1:5], verbose=TRUE)
#' @name subsetPACds
#' @family PACdataset functions
#' \code{\link{samplePACds}} to get sampled PACds.
#' @export
subsetPACds<-function(pacds,
                      group=NULL, cond1=NULL, cond2=NULL, conds=NULL,
                      avgPACtag=0, avgGeneTag=0, totPACtag=0,
                      choosePA=NULL, PAs=NULL, genes=NULL, chrs=NULL,
                      noIntergenic=FALSE, avg=FALSE, pool=FALSE,
                      clearPAT=0, minExprConds=1,
                      verbose=FALSE) {

  pacds@counts=asAnyMatrix(pacds@counts)

  if (avg & pool) stop('avg or pool')
  n1=nrow(pacds@counts)
  txt=c('before subsetPACds'=n1)

  if (is.null(group)) {
    cond1=NULL; cond2=NULL; avg=FALSE; conds=NULL

    if(avg | pool) stop("group is NULL but avg or pool is true, don't know how to summarize")
  }

  if (!is.null(choosePA)) {
    choosePA=tolower(choosePA)
    if (!(choosePA %in% c('distal','proximal','apa'))) {
      stop("choosePA must be distal/proximal/APA")
    }
  }

  if (noIntergenic) {
    ii=getNonItgFtrId(pacds@anno$ftr)
    pacds=pacds[ii]
    txt=c(txt,'noItg'=nrow(pacds@counts))
  }

  if (is.null(conds)) {
    if (!is.null(group) & (is.null(cond1) | is.null(cond2))) {
      conds=levels(pacds@colData[, group])
    } else if (!is.null(group)) {
      conds=c(cond1, cond2)
    }
  }


  if (length(conds)>0) {
    if (!AinB(conds, pacds@colData[,group]) ) {
      stop(cat(paste(conds[1], collapse = ','),'not all in pacds@colData\n'))
    }

    smps=rownames(pacds@colData)[which(pacds@colData[,group] %in% conds)]
    pacds@counts=pacds@counts[,smps,drop=F]
    pacds@colData=pacds@colData[smps,,drop=F]
    for (i in 1:ncol(pacds@colData)) {
      if (is.factor(pacds@colData[,i])) pacds@colData[,i]=droplevels(pacds@colData[,i])
    }
    row0=which(rowSums(pacds@counts)==0)
    if (length(row0)>0) {
      pacds=pacds[-row0]
    }
    txt=c(txt,'After filter conds'=nrow(pacds@counts))
  }

  if (avgPACtag>0) {
    pacds=pacds[rowMeans(pacds@counts)>=avgPACtag]
    n2=nrow(pacds@counts)
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]=paste0('avgPACtag>=',avgPACtag)
  }

  if (totPACtag>0) {
    tot=rowSums(pacds@counts)
    pacds=pacds[tot>=totPACtag]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]=paste0('totPACtag>=',totPACtag)
  }

  if (minExprConds>0) {
    nv=rowSums(pacds@counts>0)
    pacds=pacds[nv>=minExprConds]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]=paste0('minExprConds>=',minExprConds)
  }

  if (avgGeneTag>0) {
    pacds@anno=pacds@anno[rownames(pacds@counts),,drop=F]
    rs=rowSums(pacds@counts)
    b=aggregate(rs, list(gene=pacds@anno$gene), sum)
    gs=b$gene[b$x/ncol(pacds@counts)>=avgGeneTag]
    pacds=pacds[pacds@anno$gene %in% gs]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]=paste0('avgGeneTag>=',avgGeneTag)
  }

  if (!is.null(PAs)) {
    rn=rownames(pacds@counts)
    pacds=pacds[rn[rn %in% PAs]]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]='PAs'
  }

  if (!is.null(genes)) {
    gcol=NA
    if ('gene_id' %in% colnames(pacds@anno)) gcol='gene_id'
    if ('gene' %in% colnames(pacds@anno)) gcol='gene'
    if (is.na(gcol)) stop("subsetPACds by gene, but gene_id or gene not in pacds@anno\n")
    pacds=pacds[pacds@anno[, gcol] %in% genes]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]='genes'
  }


  if (!is.null(chrs)) {
    pacds=pacds[pacds@anno$chr %in% chrs]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]='chrs'
  }

  if (avg | pool) {
    if (avg) txt=c(txt,"averaging")
    if (pool) txt=c(txt,"pool")
    d=matrix()
    conds=unique(pacds@colData[, group])
    for (cond in conds) {
      if (avg) {
        di=rowMeans(pacds@counts[,rownames(pacds@colData)[which(pacds@colData[,group]==cond)], drop=F])
      } else {
        di=rowSums(pacds@counts[,rownames(pacds@colData)[which(pacds@colData[,group]==cond)], drop=F])
      }
      if (conds[1]==cond) {
        d=di
      } else {
        d=cbind(d,di)
      }
    }
    colnames(d)=conds
    rownames(d)=rownames(pacds@counts)
    pacds@counts=asAnyMatrix(d)
    pacds@colData=data.frame(group=conds)
    colnames(pacds@colData)=group
    rownames(pacds@colData)=conds
    #2023/12/29 not round because the PACds may be RUD-like ratio data  
    #pacds@counts=round(pacds@counts)
    for (i in 1:ncol(pacds@colData)) {
      pacds@colData[,i]=factor(pacds@colData[,i])
    }
  }

  #Note the filtering order!
  #If avg/pool=T, then do that first and then clearPAT
  #if a line is all 0, then remove the PAC.
  if (!is.null(clearPAT)) {
    if (clearPAT>0) {
      pacds@counts[pacds@counts<=clearPAT]=0
      idx=rowSums(pacds@counts)!=0
      pacds=pacds[idx]
      txt=c(txt, 'After clearPAT (set<=clearPAT=0)'=nrow(pacds@counts))
    }
  }

  idx=rowSums(pacds@counts)!=0
  if (sum(idx)!=length(pacds)) {
    pacds=pacds[idx]
    txt=c(txt,'After removing 0 lines'=nrow(pacds@counts))
  }

  #After clearPAT, then choosePA, because clearPAT will remove all-0 lines.
  if (!is.null(choosePA)) {
    if (choosePA %in% c('distal','proximal')) {
      pacds=get3UTRAPAds(pacds, sortPA=TRUE, choose2PA='PD')
      if (choosePA=='distal') {
        pacds=pacds[seq(2,nrow(pacds@counts),2)]
      } else {
        pacds=pacds[seq(1,nrow(pacds@counts),2)]
      }
    } else if (choosePA=='apa') {
      genes1=unique(pacds@anno$gene[duplicated(pacds@anno$gene)])
      pacds=pacds[pacds@anno$gene %in% genes1]
    }

    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]=choosePA
  }


  if (verbose) {
    txt=cbind(count=txt)
    print(txt)
  }

  return(pacds)
}


#' Sample a PACdataset
#'
#' samplePACds returns a random subset of PACdataset.
#'
#' @param pacds a PACdataset
#' @param N number of sampled PACs. If N>nrow(pacds), then do not sampling.
#' @param nPAT filter PACs with total tagnum>=nPAT before sampling.
#' @return A subset of PACdataset.
#' @examples
#' data(PACds)
#' pacds=samplePACds(PACds, N=1000)
#' @name samplePACds
#' @family PACdataset functions
#' \code{\link{subsetPACds}} to get a subset PACdataset.
#' @export
samplePACds<-function(pacds, N, nPAT=1) {
  if (nPAT>1) pacds=subsetPACds(pacds, totPACtag=nPAT)
  if (N>=nrow(pacds@counts)) return(pacds)
  idx=sample(1:nrow(pacds@counts), N)
  pacds=pacds[idx]
  return(pacds)
}


#' Get 3'UTR APA PACdataset
#'
#' get3UTRAPAds subset a PACdataset to get all PACs of genes with 3'UTR APA sites.
#'
#' Genes with 3'UTR APA sites are genes with multiple PACs in its 3'UTR.
#' If the column toStop is not in PACds@anno, then will calculate the 3'UTR length first.
#' First filter by paramters like avgPACtag, etc., and then by choose2PA.
#'
#' @param pacds A PACdataset object.
#' @param sortPA If TRUE, then order PACs by the respective 3'UTR length in each gene.
#' @param choose2PA specify whether and how to choose two PACs when there are >2 PACs.
#' The value can be NULL (use all PACs), PD (choose only proximal and distal sites),
#' farest (choose two PACs that are with the longest distance), most (choose two PACs with the most abundance).
#' @param avgPACtag if >0, then to filter by PAC tag num, >=avgPACtag, see subsetPACds().
#' @param avgGeneTag if >0, then to filter by PAC tag num, >=avgGeneTag, see subsetPACds().
#' @param clearPAT if >0, then to filter by clearPAT, see subsetPACds().
#' @return A PACdataset with only 3'UTR APA sites. If there is no result, return an empty PACdataset with 0-row @anno and @counts.
#' @examples
#' data(PACds)
#' ## Get all 3'UTR APA
#' pp1=get3UTRAPAds(PACds,sortPA=TRUE); summary(pp1); is3UTRAPAds(pp1); is3UTRAPAds(PACds)

#' ## Get proximal distal PA
#' pp1=get3UTRAPAds(PACds, sortPA=TRUE, choose2PA='most');
#' summary(pp1); is3UTRAPAds(pp1); is3UTRAPAds(PACds); head(pp1@anno)
#' pp1@anno[pp1@anno$gene=='PH02Gene00018',]
#' @name get3UTRAPAds
#' @family PACdataset functions
#' @export
get3UTRAPAds<-function(pacds, sortPA=TRUE, choose2PA=NULL, avgPACtag=0, avgGeneTag=0, clearPAT=0) {

  if (!is.null(choose2PA)) {
    sortPA=TRUE
    choose2PA=toupper(choose2PA)
    if (!(choose2PA %in% c('PD','MOST'))) stop("choose2PA must be PD or MOST\n")
  }

  #filter by tagnum
  pacds=subsetPACds(pacds, avgPACtag=avgPACtag, avgGeneTag=avgGeneTag, noIntergenic=TRUE, clearPAT=clearPAT, verbose=FALSE)

  #filter 3UTR APA
  if (!is3UTRAPAds(pacds)) {
    pacds@anno=pacds@anno[which(pacds@anno$ftr=='3UTR'),]
    genes=unique(pacds@anno$gene[duplicated(pacds@anno$gene)])
    if (length(genes)==0) { #no 3utr APA
      pacds@anno=pacds@anno[character(0), ]
      pacds@counts=pacds@counts[character(0), ]
      return(pacds)
    }
    pacds@anno=pacds@anno[pacds@anno$gene %in% genes,]
    if (!('toStop' %in% colnames(pacds@anno))) {
      if ('three_UTR_length' %in% colnames(pacds@anno)) {
        pacds@anno$toStop=pacds@anno$three_UTR_length
      } else {
        pacds@anno$toStop=pacds@anno$coord
        id=grep('\\+',pacds@anno$strand)
        pacds@anno$toStop[id]=pacds@anno$coord[id]-pacds@anno$ftr_start[id]+1
        id=grep('\\-',pacds@anno$strand)
        pacds@anno$toStop[id]=pacds@anno$ftr_end[id]-pacds@anno$coord[id]+1
      }
    }
    if (min(pacds@anno$toStop)<0) stop("get3UTRAPAds: error toStop (<0), please check coord/ftr_start/ftr_end\n")
    pacds@counts=pacds@counts[rownames(pacds@anno), ]
  }

  #sort by 3UTR len
  if (sortPA) {
    pacds=pacds[order(pacds@anno$gene, pacds@anno$toStop,decreasing = FALSE)]
  }

  if (!is.null(choose2PA)) {

    rs=rowSums(pacds@counts)
    d=pacds@anno[,c('gene','toStop')]
    d$tag=rs
    d$PA=rownames(pacds@anno)

    .choose2<-function(x, method) {
      if (nrow(x)==2) {
        return(x$PA)
      }
      if (method=='PD') {
        return(x$PA[c(1,nrow(x))])
      } else if (method=='MOST') {
        x=x[order(x$tag,decreasing = TRUE),]
        return(x$PA[c(1:2)])
      }
    }

    sp=split(d, d$gene, drop=T)
    PAs=unlist(lapply(sp, .choose2, method=choose2PA))
    #x=matrix(unlist(PAs), ncol=2, byrow = TRUE)
    #PAs=.asDf(cbind(gene=names(PAs), x))

    #PAs=plyr::ddply(d,.(gene), .choose2, method=choose2PA)
    #PAs=unlist(PAs[,2:3])

    pacds@anno=pacds@anno[PAs,]
    pacds@counts=pacds@counts[PAs,]
    if (sortPA) {
      pacds=pacds[order(pacds@anno$gene, pacds@anno$toStop,decreasing = FALSE)]
    }
  }

  return(pacds)
}


#' Get distance between nearby pAs
#'
#' getNearbyPAdist calcualtes distance between adjacent pAs within each gene or within a given `ftr` (e.g., 3UTR) of each gene.
#' This function can be used to determine the `minDist` and `maxDist` parameters in \code{\link{get3UTRAPApd}}.
#'
#' @param pacds a PACdataset, which should contain gene and/or ftr columns in its anno slot.
#' @param ftrs default is 3UTR, used to filter pacds. It couls be a vector like c('3UTR','intron').
#' If it is NULL, then will not filter pacds by ftr.
#' @param within should be gene (default) or ftr, which means to calculate distance within each gene or ftr (but still within each gene).
#' @param verbose TRUE to print the summary of distance.
#' @param return can be summary/mean/median/df/other any character. mean to return the mean distance value;
#' summary (default) to return summary of distance; median to return median distance; df to return the data.frame;
#' others to not return anything but just print summary.
#' @return Depends on the `return` parameter. If return='df', then a data.frame with two columns (id, distance). The id depends on the `within` parameter,
#' would be like `chrX.-.1756` (chr.strand.gene) or `chrX.-.1756.3UTR` (chr.strand.gene.ftr).
#' @examples
#' data(PACds)
#' ## summary of distance
#' s=getNearbyPAdist(PACds, ftrs='3UTR', within='gene')
#' ## mean distance
#' getNearbyPAdist(PACds, ftrs='3UTR', within='gene', return='mean')
#' ## all distance in a data frame
#' d=getNearbyPAdist(PACds, ftrs='3UTR', within='gene', return='df')
#' ## median distance, and not print summary
#' d[d$distance==max(d$distance), ]
#' getNearbyPAdist(PACds, ftrs='3UTR', within='gene', return='median', verbose=FALSE)
#' ## just print summary, not return anything
#' getNearbyPAdist(PACds, ftrs='3UTR', within='gene', return='null')
#' ## not filter ftrs, but within=ftr, this will calculate within each gene's each ftr
#' getNearbyPAdist(PACds, ftrs=NULL, within='ftr', return='null')
#' @name getNearbyPAdist
#' @seealso \code{\link{get3UTRAPApd}}
#' @export
getNearbyPAdist<-function(pacds, ftrs='3UTR', within='gene', return='summary', verbose=TRUE) {

  if (!is.null(ftrs)) {
    if (!('ftr' %in% colnames(pacds@anno))) stop("nearbyAPAdist: ftr column not in pacds@anno, cannot filter by ftr!\n")
    pacds=pacds[pacds@anno$ftr %in% ftrs]
  }

  within=tolower(within)
  if (within=='gene') {
    if (!('gene' %in% colnames(pacds@anno))) stop("nearbyAPAdist: gene column not in pacds@anno, cannot use within=gene!\n")
    d=pacds@anno
    d=pacds@anno[order(d$chr, d$strand, d$gene, d$coord), c('chr','strand', 'gene', 'coord')]
    d=d[d$gene!='character(0)', ]
    sp=split(d, d[,c('chr', 'strand', 'gene')], drop=T)
  } else if (within=='ftr') {
    if (!all(c('gene','ftr') %in% colnames(pacds@anno))) stop("nearbyAPAdist: ftr/gene columns not in pacds@anno, cannot use within=ftr!\n")
    d=pacds@anno
    d=pacds@anno[order(d$chr, d$strand, d$gene, d$ftr, d$coord), c('chr','strand', 'gene', 'ftr','coord')]
    d=d[d$gene!='character(0)', ]
    sp=split(d, d[,c('chr', 'strand', 'gene', 'ftr')], drop=T)
  } else {
    stop("nearbyAPAdist: within should be gene or ftr, which means get distance of adjacent pAs within each gene or within ftr of each gene!\n")
  }

  spdiff=lapply(sp, FUN=function(x) abs(diff(x$coord)))

  s=summary(unlist(spdiff))

  return=tolower(return)
  if (return=='mean') {
    out=s[4]
  } else if (return=='median') {
    out=s[3]
  } else if (return=='df') {
    spid=rep(names(spdiff), Biobase::listLen(spdiff))
    out=data.frame(id=spid, distance=unlist(spdiff))
    rownames(out)=NULL
  } else if (return=='summary') {
    out=s
  } else {
    out=NULL
    verbose=TRUE
  }

  if (verbose) {
    cat('Nearby pA distance summary:\n')
    print(s)
    cat('\n')
  }

  if (is.null(out)) return(invisible(NULL))
  return(out)
}


#' Get proximal and distal 3'UTR APAs
#'
#' get3UTRAPApd provides a smarter way to get proximal and distal 3'UTR APAs then get3UTRAPAds.
#' Criteria for a proximal and distal pair in a gene: first select APA pairs from all APA pairs that meet the dist and ratio criteria; then the final P/D pair is the one with the highest count sum and the closest distance.
#' Note that there is no filtering of PA expression levels here, which can be performed before get3UTRAPApd using subsetPACds().
#'
#' @param pacds a PACdataset, if it contains non-3UTR PAs, will run get3UTRAPAds first to subset 3UTR PAs.
#' @param minDist the distance of the proximal and distal PA should be in [minDist, maxDist].
#' @param maxDist the distance of the proximal and distal PA should be in [minDist, maxDist].
#' It is recommended to use \code{\link{getNearbyPAdist}} to estimate an appropriate distance.
#' @param minRatio If minRatio>0, then will remove PAs with ratio<minRatio before getting P/D PAs.
#' @param fixDistal TRUE means that the most distal PA is fixed and then only searches for the proximal PA for that distal PA.
#' @param addCols If not NULL, then four columns would be added to PACds@anno, inlucidng <addCols>Which/Score/Ratio/Dist. Among these columns, the 'which' column denotes P or D of each PA.
#' @return A PACdataset with only 3'UTR P-D APA sites. If there is no result, return an empty PACdataset with 0-row @anno and @counts.
#' @examples
#' data(PACds)
#' ## to get a glance of the distance between nearby 3UTR pAs
#' getNearbyPAdist(PACds, ftrs='3UTR', within='ftr', return='null')
#' pd1=get3UTRAPApd(PACds, minDist=50, maxDist=1000, minRatio=0.05, fixDistal=F)
#' pd2=get3UTRAPApd(PACds, minDist=50, maxDist=1000, minRatio=0.05, fixDistal=T)
#' g=PACds@anno$gene[!(pacds@anno$gene %in% pd1@anno$gene)]
#' subsetPACds(pd1, gene=g)
#' subsetPACds(pd2, gene=g)
#' @name get3UTRAPApd
#' @family PACdataset functions
#' @seealso \code{\link{subsetPACds}} \code{\link{get3UTRAPAds}} \code{\link{movAPAindex}}
#' @export
get3UTRAPApd<-function(pacds, minDist=50, maxDist=1e4, minRatio=0.05, fixDistal=FALSE, addCols='pd') {

  pacds@counts=asAnyMatrix(pacds@counts)

  if ('ftr' %in% colnames(pacds@anno)) {
    ftrs=unique(pacds@anno$ftr)
    if (length(ftrs)!=1 || ftrs!='3UTR') {
      n1=length(pacds)
      pacds=get3UTRAPAds(pacds)
      cat( sprintf("get3UTRAPApd: run get3UTRAPAds to get %d 3UTR APAs from %d PAs\n", length(pacds), n1) )
    }
  }

  if (!('gene' %in% colnames(pacds@anno))) stop("get3UTRAPApd: gene is not in pacds@anno, please use annotatePAC to annotate the pacds first!")

  nc=rowSums(pacds@counts)
  annos=cbind(pacds@anno[, c('chr','strand','coord', 'gene')], nc=nc)
  annos$PA=rownames(annos)

  gnc=tapply(annos$nc, as.factor(annos$gene), sum)
  gnp=tapply(annos$gene, as.factor(annos$gene), length)
  g=data.frame(gene=names(gnc), gnc, gnp)

  annos=merge(annos, g, by='gene')
  annos$ratio=annos$nc/annos$gnc

  if (minRatio>0) {
    ridx=(annos$ratio>minRatio)
    if (sum(ridx)>0) {
      annos=annos[ridx, ]
      gnp=tapply(annos$gene, as.factor(annos$gene), length)
      gpass=names(gnp)[gnp>=2]
      cat('get3UTRAPApd: filtering by minRatio: gene# before:', length(gnp), '; after:', length(gpass),'; remove:', length(gnp)-length(gpass), '\n')
      annos=annos[annos$gene %in% gpass, ]
    }
  }

  # gene/chr/strand/coord
  .sortAnnos<-function(annos) {
    anno1=annos[annos$strand=='+', ]
    anno1=anno1[order(anno1$gene, anno1$chr, -anno1$coord), ]
    anno2=annos[annos$strand=='-', ]
    anno2=anno2[order(anno2$gene, anno2$chr, anno2$coord), ]
    annos=rbind(anno1,anno2)
    return(annos)
  }

  annos=.sortAnnos(annos)

  ####
  .getAPApd<-function(anno) {

    # ?̶?? [1]?
    if (fixDistal) {
      anno$dist=abs(anno$coord-anno$coord[1])
      idx=which(anno$dist>0 & anno$dist<=maxDist & anno$dist>=minDist)
      if (length(idx)>0) {     mid=idx[anno$nc[idx]==max(anno$nc[idx])][1]
        return(anno$PA[c(1, mid)])
      } else {
        return(c())
      }
    }

    # not fixDistal
    d=c()
    r=c()
    p1=c()
    p2=c()
    for (i in 1:(nrow(anno)-1)) {
      for (j in (i+1):(nrow(anno))) {
        p1=c(p1, anno$PA[i]); p2=c(p2, anno$PA[j])
        d=c(d, abs(anno$coord[i]-anno$coord[j]))
        r=c(r, anno$ratio[i]+anno$ratio[j])
      }
    }
    d=data.frame(p1=p1, p2=p2, dist=d, ratio=r)
    idx1=which(d$dist>0 & d$dist<=maxDist & d$dist>=minDist)
    idx2=which(d$ratio==max(d$ratio))
    idx=base::intersect(idx1, idx2)
    if (length(idx)==1) {
      return(c(d$p1[idx], d$p2[idx]))
    } else if (length(idx)==0) {
      return(c())
    } else { #idx[which(d$dist[idx]==min(d$dist[idx]))]
      idx3=idx[which(d$dist[idx]==min(d$dist[idx]))][1]
      return(c(d$p1[idx3], d$p2[idx3]))
    }
    return(c())
  }
  ####

  sp=split(annos, annos$gene, drop=T)

  pas=unlist(lapply(sp, .getAPApd))
  cat('get3UTRAPApd: filtering pd (dist between): gene# before:', length(sp), '; after:', length(pas)/2,'; remove:', length(sp)-length(pas)/2, '\n')

  if (length(pas)==0) {
    pacds@counts=pacds@counts[-(1:nrow(pacds@counts)), ]
    pacds@anno=pacds@anno[-(1:nrow(pacds@anno)), ]
    return(pacds)
  } else {

    pacds@anno=pacds@anno[pas, ]
    pacds@anno=.sortAnnos(pacds@anno)
    pacds@counts=pacds@counts[rownames(pacds@anno), ]

    if (!is.null(addCols)) {
      newcols=paste0(addCols[1], c('Which','Score','Ratio','Dist'))
      cat("get3UTRAPApd: add four columns to pacds@anno:", toString(newcols), '\n')
      pacds@anno[, newcols[1]]=rep(c('D','P'),length.out=nrow(pacds@anno))

      nc=rowSums(pacds@counts)
      idD=seq(1, nrow(pacds@anno), by=2)
      idP=idD+1
      ncg=rep(nc[idD]+nc[idP], each=2)

      pacds@anno[, newcols[2]]=nc
      pacds@anno[, newcols[3]]=nc/ncg

      pacds@anno[, newcols[4]]=rep(abs(pacds@anno$coord[idD]-pacds@anno$coord[idP]), each=2)
    }
    return(pacds)
  }
}



# If a PACds is 3UTRAPA
# Judge by only ftr=3UTR, have toStop column, and APA
# will not order PACds
is3UTRAPAds<-function(pacds) {
  if (nrow(pacds@anno)==0) return(FALSE)
  if (sum(pacds@anno$ftr!='3UTR')>0) return(FALSE)
  if (!('toStop' %in% colnames(pacds@anno))) return(FALSE)
  if (min(table(pacds@anno$gene))<=1) return(FALSE) #not APA
  return(TRUE)
}


# Stat gene info for a PACds
# For each column in @counts
# @param avgUTRlen PACds should have only 3UTR ftr and toStop column
# @param geneTag Get nPAC and tag number of genes
# @param PAinfo For each gene in each cond, get info like coord1=tag1;coord2=tag2.
# @param ftr For each gene in each cond, get info like "3UTR, intergenic".
# @param merge: if TRUE then merge list results to one matrix (gene number must be the same)
# @return list[[avgUTRlen]]=(gene, cond1, cond2...), [[geneTag]]=... (nPAC, cond1, cond2...), [[PAinfo]]=(... /53609464=5;53609939=118/)
#       merge=TRUE --> c('gene','avgUTRlen1', 'avgUTRlen2','nPAC','geneTag1', 'geneTag2','PAs1', 'PAs2')
geneStatForPACds<-function(pacds, avgUTRlen=TRUE, geneTag=TRUE, PAinfo=TRUE, ftr=FALSE, merge=TRUE) {

  ls=list()

  pacds@counts=asAnyMatrix(pacds@counts)

  pacds@counts=round(pacds@counts)

  .statUTR<-function(x) {
    tot=sapply(x[,-c(1:2)],function(xx, tostop) {sum(xx*tostop)}, tostop=x$toStop)
    tot=tot/colSums(x[,-c(1:2)])
    tot=data.frame(matrix(tot, nrow=1))
  }

  .statTag <- function(x) {
    data.frame(matrix( c(nrow(x), sapply(x[,-1], sum) ), nrow=1))
  }

  if (avgUTRlen) {

    if (!is3UTRAPAds(pacds)) {
      stop("avgUTRlen=TRUE, but pacds is non-3UTR-APA dataset\n")
    }

    d=cbind(pacds@anno[,c('gene','toStop')], .asDf(pacds@counts))
    d=d[order(d$gene),]

    #avgUTR=plyr::ddply(d, .(gene), .statUTR )

    avgUTR=.asDf(d %>% dplyr::group_by(gene) %>% dplyr::do(.statUTR(.)))
    colnames(avgUTR)=c('gene', colnames(pacds@counts))

    ls[['avgUTRlen']]=avgUTR

    if (merge)   gs=ls$avgUTRlen
  }

  if (geneTag) {
    d=cbind(gene=pacds@anno$gene, .asDf(pacds@counts))
    d=d[order(d$gene),]

    #tag1=plyr::ddply(d,.(gene), function(x) c('nPAC'=nrow(x),sapply(x[,-1], sum) ) )
    tag=.asDf(d %>% dplyr::group_by(gene) %>% dplyr::do(.statTag(.)))
    colnames(tag)=c('gene','nPAC',colnames(d)[-1])

    ls[['geneTag']]=tag

    if (merge) {
      if (avgUTRlen) {
        if (nrow(ls$geneTag)!=nrow(ls$avgUTRlen)) {
          stop(cat('merge=TRUE but ls$geneTag and ls$avgUTRlen not equal gene number\n'))
        }
        ls$geneTag=ls$geneTag[,-1]
        colnames(ls$geneTag)[-1]=paste0('geneTag_',colnames(ls$geneTag)[-1])
        gs=cbind(gs, ls$geneTag)
      } else {
        gs=ls$geneTag
      }
    }
  }

  .getPAinfo<-function(x) {
    info=sapply(x[,-c(1:2)],function(xx, PA) {paste(PA,xx,sep='=')}, PA=x$PA)
    if (!is.vector(info)) {
      info=apply(info, 2, paste, collapse=';')
    }
    info=data.frame(matrix(info, nrow=1))
  }

  if (PAinfo) {
    d=cbind(gene=pacds@anno$gene, PA=rownames(pacds@anno), .asDf(pacds@counts))
    d=d[order(d$gene),]

    #info1=plyr::ddply(d, .(gene), .getPAinfo )

    info=.asDf(d %>% dplyr::group_by(gene) %>% dplyr::do(.getPAinfo(.)))
    colnames(info)=c('gene', colnames(pacds@counts))

    ls[['PAinfo']]=info

    if (merge) {
      if (avgUTRlen | geneTag) {
        if (nrow(ls$geneTag)!=nrow(ls$PAinfo)) {
          stop(cat('merge=TRUE but ls$geneTag and ls$PAinfo not equal gene number\n'))
        }
        ls$PAinfo=ls$PAinfo[-1]
        colnames(ls$PAinfo)[-1]=paste0('PAs_',colnames(ls$PAinfo)[-1])
        gs=cbind(gs, ls$PAinfo)
      } else {
        gs=ls$PAinfo
      }
    }
  }

  if (ftr) {
    d=.asDf(cbind(gene=pacds@anno$gene, ftr=pacds@anno$ftr))
    d=d[order(d$gene),]

    #ftr1=plyr::ddply(d,.(gene), function(x) paste(unique(sort(x$ftr)), collapse = ',') )

    ftr=.asDf(d %>% dplyr::group_by(gene) %>% dplyr::do(
        data.frame( V1=matrix (paste(unique(sort(.$ftr)), collapse = ','), nrow=1))
      )
    )

    colnames(ftr)=c('gene','ftr')
    ls[['ftr']]=ftr


    if (merge) {
      if (avgUTRlen | geneTag) {
        if (nrow(ls$geneTag)!=nrow(ls$ftr)) {
          stop(cat('merge=TRUE but ls$geneTag and ls$ftr not equal gene number\n'))
        }
        gs=cbind(gs, ftr=ls$ftr[-1])
      } else {
        gs=ls$ftr
      }
    }
  }

  if (merge) {
    return(gs)
  } else {
    return(ls)
  }
}

#Convert a PACds to dataframe [chr,strand,coord, samples...], rownames=rownames(PACds@counts)
PACds2PAdf <- function(PACds) {
  if (!AinB(c('chr','strand','coord'), colnames(PACds@anno))) stop('chr/strand/coord not in PACds@anno')
  if (inherits(x = PACds@counts, what = 'dgCMatrix')) {
    if (ncol(PACds@counts)>500 & nrow(PACds@counts)>1e4)
      warning("The PACds@counts is a big dgCMatrix[>10000 rows, >500 columns,] which may need BIG mem!",
        call. = FALSE, immediate. = TRUE)
    PACds@counts=as.matrix(PACds@counts)
  }
  df=cbind(PACds@anno[,c('chr','strand','coord')], PACds@counts, stringsAsFactors = FALSE)
  return(df)
}

# Merge multiple PACdataset
#
# mergePACds groups nearby PACs from single/multiple PACdataset objects.
#
# This function is particularlly useful for grouping nearby cleavage sites into PACs.
# If is also useful When you have multiple PA or PAC files, each file is from one sample.
# Then you need to merge these PACds into one PACds for DE or other analyses.
# But after grouping and/or merging, you may need to annotatePAC by a GFF annotation.
# @usage mergePACds(PACdsList, d=24)
# @param PACdsList a PACdataset, or a list of multiple PACdataset objects. The PACds@anno should have columns chr/strand/coord.
# If there is no colData in PACds, then will be set as groupN.
# If PACdsList is a PACdataset, then will treat it as PA and group nearby PAs into PACs.
# @param d distance to group nearby PACds, default is 24 nt.
# @return A merged PACdataset. The counts slot stores counts of merged samples.
# If sample names from different PACdataset objects are duplicated, then will be set as .x,.y.
# The colData slot stores the merge sample annotation from the first column of each @colData.
# The anno slot will have these columns: chr, strand, coord, tottag, UPA_start, UPA_end, nPA, maxtag.
# @examples
# ## group PA into PACs
# data(PACds)
# PACds@counts=rbind(PACds@counts, PACds@counts)
# PACds@anno=rbind(PACds@anno, PACds@anno)
# ds=mergePACds(PACds, d=24)
# ## merge two PACds
# pacds=new('PACdataset', counts=PACds@counts, anno=PACds@anno)
# PACdsList=list(PACds, pacds)
# ds=mergePACds(PACdsList, d=24)
mergePACds_old <- function (PACdsList, d=24) {

  ## from .r script ##
  #library(GenomicRanges, verbose = FALSE)
  #library(data.table, verbose = FALSE)

  if (inherits(PACdsList, "PACdataset")) PACdsList=list(PACdsList)

  for (i in 1:length(PACdsList)) {
    if (!AinB(c('chr','strand','coord'), colnames(PACdsList[[i]]@anno))) stop('chr/strand/coord not in PACdsList@anno')
    if (nrow(PACdsList[[i]]@colData)==0) {
      PACdsList[[i]]@colData=.asDf(matrix( rep(paste0('group',i),
                                                       ncol(PACdsList[[i]]@counts)), ncol=1,
                                                   dimnames =list(colnames(PACdsList[[i]]@counts),'group') ))
    }
  }

  allpa=PACds2PAdf(PACdsList[[1]])

  if (length(PACdsList)>=2) {
    for (j in 2:length(PACdsList)) {
      pa2=PACds2PAdf(PACdsList[[j]])
      allpa=merge(allpa,pa2, all=T, by.x=c('chr','strand','coord'), by.y=c('chr','strand','coord'))
    }
  }
  allpa[is.na(allpa)]=0
  invisible(gc())

  ## allpa group by dist
  allpa=allpa[order(allpa$chr,allpa$strand,allpa$coord),]
  invisible(gc())
  diffcrd=diff(allpa$coord) #distance
  same=(allpa$chr[1:(nrow(allpa)-1)]==allpa$chr[2:nrow(allpa)]) & (allpa$strand[1:(nrow(allpa)-1)]==allpa$strand[2:nrow(allpa)])
  #allpa[!same,c('chr','strand')]
  indist=(diffcrd<=d & same) #distance<X and the same chr-strand
  invisible(gc())
  indist=as.character(indist)
  indist[indist=='FALSE']=paste('F',1:length(indist[indist=='FALSE']),sep='')
  rle=S4Vectors::Rle(indist)
  cum=cumsum(runLength(rle))
  Tto=cum[S4Vectors::runValue(rle)=='TRUE']+1
  Tfrom=cum[which(S4Vectors::runValue(rle)=='TRUE')-1]+1
  if (length(Tfrom)==length(Tto)-1) {
    Tfrom=c(1,Tfrom)
  }
  if (length(Tfrom)!=length(Tto)) {
    stop("error: length(Tfrom)!=length(Tto)\n")
  }
  #write.table(cbind(Tfrom,Tto,Tto-Tfrom),file='tfrom_tTo.txt')
  ir=IRanges::IRanges(start=Tfrom,end=Tto)
  gap=gaps(ir,start=1,end=nrow(allpa))
  gap=.asDf(gap)
  Ffrom=unlist(apply(gap,1,function(par) return(seq(par[1],par[2]))))
  Fto=Ffrom
  groups=rbind(cbind(from=Tfrom,to=Tto),cbind(from=Ffrom,to=Fto))
  groups=groups[order(groups[,'from']),]
  gnames=paste('g',1:nrow(groups),sep='')
  pacChrStrand=allpa[groups[,'from'],c('chr','strand')]
  groups=.asDf(cbind(groups,pacChrStrand,gnames))
  groups$from=as.numeric(groups$from)
  groups$to=as.numeric(groups$to)
  gnames=rep(gnames,times=groups[,'to']-groups[,'from']+1)
  rm(rle,indist,gap,ir,cum,Tto,Tfrom,Fto,Ffrom,pacChrStrand)
  invisible(gc())
  ##groups=from to  chr strand gnames

  #cat('Calculating nPA, PACtag (slow) ...\n')
  smpcols=colnames(allpa)[!(colnames(allpa) %in% c('chr','strand','coord'))]

  if (length(smpcols)>1) {
    tottag=rowSums(allpa[,smpcols])
  } else {
    tottag=allpa[,smpcols]
  }
  allpa$tottag=tottag
  allpa$gnames=gnames
  UPA_start=allpa$coord[groups[,'from']]
  UPA_end=allpa$coord[groups[,'to']]
  nPA=groups[,'to']-groups[,'from']+1
  groups=cbind(groups,UPA_start,UPA_end,nPA)


  sgroups=groups[nPA==1,]
  mgroups=groups[nPA>1,]
  sallpa=allpa[allpa$gnames %in% sgroups$gnames,]
  mallpa=allpa[allpa$gnames %in% mgroups$gnames,]
  rm(groups,allpa)

  if (nrow(mallpa)>0) {
    smpcols=c(smpcols,'tottag')
    x=gc(verbose =F)
    mpactag=aggregate(mallpa[,smpcols], list(groups=mallpa$gnames), sum)
    x=gc(verbose =F)

    maxtag=aggregate(mallpa[,'tottag'], list(groups=mallpa$gnames), max)
    colnames(maxtag)=c('groups','tottag')
    x=gc(verbose =F)
    maxtag=data.table::data.table(maxtag)
    mallpa=data.table::data.table(mallpa)
    m=merge(maxtag,mallpa[,c('gnames','coord','tottag'),with=FALSE],by.x=c('groups','tottag'),by.y=c('gnames','tottag'))
    x=gc(verbose =F)
    m=m[order(m$groups),]
    rle=S4Vectors::Rle(m$groups)
    cum=cumsum(runLength(rle))
    m=m[cum,]
    colnames(m)=c('groups','maxtag','coord')
    mgroups=merge(mgroups,m,by.x='gnames',by.y='groups')
    mgroups=merge(mgroups,mpactag,by.x='gnames',by.y='groups')
  }

  sallpa$maxtag=sallpa$tottag
  sgroups=sgroups[,!(colnames(sgroups) %in% c('chr','strand'))]
  sgroups=merge(sgroups,sallpa,by.x='gnames',by.y='gnames')
  if (nrow(mgroups)>0)  {
    sgroups=sgroups[,colnames(mgroups)]
    groups=rbind(sgroups,mgroups)
  } else {
    groups=sgroups
  }

  smpcols=smpcols[smpcols!='tottag']
  cols=c('chr','strand','coord','tottag','UPA_start','UPA_end','nPA','maxtag',smpcols)
  groups=groups[,cols]

  counts=groups[,smpcols]
  colData=.asDf(matrix(unlist(lapply(PACdsList, function(ds) return(as.character(ds@colData[,1])))),
                               ncol=1, dimnames = list(colnames(counts),'group')))
  anno=groups[,c('chr','strand','coord','tottag','UPA_start','UPA_end','nPA','maxtag')]

  d=new("PACdataset",counts=counts, colData=colData, anno=anno)
  return(d)
}



#' Merge multiple PACdatasets for RNA-seq
#'
#' mergePACds_v0 groups nearby PACs from single/multiple PACdataset objects.
#' This function is the mergePACds in previous version of movAPA, which is not recommanded to use now because it is not very suitable for (large) single-cell data.
#' This function is particularly useful for grouping nearby cleavage sites into PACs.
#' It is also useful When you have multiple PA or PAC files, each file is from one sample.
#' Then you need to merge these PACds into one PACds for DE or other analyses.
#' But after grouping and/or merging, you may need call annotatePAC to annotate the merged PACs by a GFF annotation.
#' @param PACdsList is a PACdataset, or a list of multiple PACdataset objects. The PACds@anno should have columns chr/strand/coord.
#' If there is no colData in PACds, then the sample label will be set as groupN.
#' If PACdsList is a PACdataset, then will treat it as PA and group nearby PAs into PACs.
#' @param d distance to group nearby PACds, default is 24 nt.
#' @return A merged PACdataset. The counts slot stores counts of merged samples.
#' If sample names from different PACdataset objects are any duplicated, then the sample name will be added a suffix .i for each sample in PACds[[i]].
#' The colData slot stores the merge sample annotation from the first column of each @colData.
#' The anno slot contains these columns: chr, strand, coord, tottag, UPA_start, UPA_end, nPA, maxtag.
#' @examples
#' ## Group PA into PACs
#' data(PACds)
#' PACds@counts=rbind(PACds@counts, PACds@counts)
#' PACds@anno=rbind(PACds@anno, PACds@anno)
#' ds=mergePACds_v0(PACds, d=24)
#' ## merge two PACds
#' ds1=makeExamplePACds()
#' ds2=makeExamplePACds()
#' ds=mergePACds_v0(list(ds1, ds2), d=24)
#' @name mergePACds_v0
#' @family PACdataset functions
#' \code{\link{annotatePAC}} to annotate a PACdataset; [rbind()] to combine multiple PACdatasets of the same format.
#' @export
mergePACds_v0 <- function (PACdsList, d=24) {

  if (inherits(PACdsList, "PACdataset")) PACdsList=list(PACdsList)

  for (i in 1:length(PACdsList)) {
    if (!AinB(c('chr','strand','coord'), colnames(PACdsList[[i]]@anno))) stop('chr/strand/coord not in PACdsList@anno')
    if (nrow(PACdsList[[i]]@colData)==0) {
      PACdsList[[i]]@colData=.asDf(matrix( rep(paste0('group',i),
                                                       ncol(PACdsList[[i]]@counts)), ncol=1,
                                                   dimnames =list(colnames(PACdsList[[i]]@counts),'group') ))
    }
  }

  ### format colData to avoid duplicated sample names
  if (length(PACdsList)>1) {
    sn=unlist(lapply(PACdsList, function(ds) return(rownames(ds@colData))))
    snu=unique(sn)
    if (length(sn)!=length(snu)) {
      cat(sprintf("mergePACds_v0: there are %d duplicated sample names in the PACdsList, will add .N to sample names of each PACds\n", length(sn)-length(snu)))
      #sufx=c('x','y','z',letters[1:23])
      for (i in 1:length(PACdsList)) {
        colnames(PACdsList[[i]]@counts)=paste0(colnames(PACdsList[[i]]@counts), '.', i)
        rownames(PACdsList[[i]]@colData)=paste0(rownames(PACdsList[[i]]@colData), '.', i)
      }
    }
  }

  #colnames used for restore order of sample columns
  cnames=unlist(lapply(PACdsList, function(d) return(colnames(d@counts))))

  allpa=PACds2PAdf(PACdsList[[1]])

  if (length(PACdsList)>=2) {
    for (j in 2:length(PACdsList)) {
      pa2=PACds2PAdf(PACdsList[[j]])
      allpa=merge(allpa,pa2, all=T, by.x=c('chr','strand','coord'), by.y=c('chr','strand','coord'))
    }
  }
  allpa[is.na(allpa)]=0
  invisible(gc())


  gr <- with(allpa, GRanges(seqnames = chr,
                                   ranges =IRanges::IRanges(start=coord,
                                                   end=coord),
                                   strand = strand) )

  #resize width as d+1
  gr=GenomicRanges::resize(gr, width=d+1, fix="start", use.names=TRUE, ignore.strand=FALSE)

  itv=GenomicRanges::reduce(gr, drop.empty.ranges=TRUE)

  cat("mergePACds_v0: group PA to PACs\n")
  ov = GenomicRanges::findOverlaps(gr, itv,
                    maxgap=-1L, minoverlap=1L,
                    type=c("any"), select='all',
                    ignore.strand=FALSE)
  ov=.asDf(ov)
  allpa$idx=1:nrow(allpa)
  allpa=merge(allpa, ov, by.x='idx', by.y='queryHits')

  allpa$idx=NULL
  smpcols=colnames(allpa)[!(colnames(allpa) %in% c('chr','strand','coord','subjectHits'))]

  #sum tag per interval
  cat('count tot tag for each sample within each PAC\n')
  allpa$tottag=rowSums(allpa[, smpcols, drop=F])
  byItv <- dplyr::group_by(allpa, subjectHits)
  dots <- sapply(smpcols ,function(x) substitute(sum(x), list(x=as.name(x))))
  dots[['tottag']]=substitute(sum(x),list(x=as.name('tottag')))
  pac=do.call(dplyr::summarise, c(list(.data=byItv), dots))

  cat('mergePACds_v0: annotate the range of each PAC\n')
  #get interval range...
  pacAnno=byItv  %>% dplyr::summarise(nPA=n(), UPA_start=min(coord), UPA_end=max(coord), maxtag=max(tottag), coord=coord[which.max(tottag)], chr=chr[1], strand=strand[1])

  pac=merge(pac, pacAnno, by.x='subjectHits', by.y='subjectHits')

  #pacds
  smpcols=smpcols[smpcols!='tottag']
  counts=pac[, smpcols, drop=F]
  anno=pac[,c('chr','strand','coord','tottag','UPA_start','UPA_end','nPA','maxtag')]

  colData=.asDf(matrix(unlist(lapply(PACdsList, function(ds) return(as.character(ds@colData[,1])))),
                               ncol=1, dimnames = list(colnames(counts),'group')))

  if(all(cnames %in% colnames(counts))) {
    counts=counts[, cnames , drop=F]
    colData=colData[cnames,  , drop=F]
  }

  d=createPACdataset(counts=counts, colData=colData, anno=anno)
  return(d)
}


## check columns of annos
## by: coord or T(coord)/F
## ifnotstop, otherwise return NULL
.checkPACdsBy<-function(PACdsList, by='coord', verbose=TRUE) {

  if (is.logical(by))
    bycoord=by
  else
    bycoord=ifelse(tolower(by)=='coord', TRUE, FALSE)

  if (inherits(PACdsList, "PACdataset")) PACdsList=list(PACdsList)

  for (i in 1:length(PACdsList)) {
    if (bycoord) {
      if (!AinB(c('chr','strand','coord'), colnames(PACdsList[[i]]@anno))) stop('mergePACds: by=coord, but chr/strand/coord not in PACdsList@anno')
    }
    else {
      if (!AinB(c('chr','strand','UPA_start','UPA_end'), colnames(PACdsList[[i]]@anno))) stop('mergePACds: by=range, but chr/strand/UPA_start/UPA_end not in PACdsList@anno')
    }
  }
  if (verbose & !bycoord) cat("by is not coord, then use UPA_start~UPA_end to merge PACds.\n")
  return(invisible(NULL))
}

### PA ranges table: PA/chr/strand/..start/end
.getPAranges<-function(p, bycoord=TRUE, pn='') {
  if (bycoord) {
    p=p@anno[, c('chr','strand', 'coord')]
    p$UPA_start=p$coord
    p$UPA_end=p$coord
    p$coord=NULL
  } else {
    p=p@anno[, c('chr','strand', 'UPA_start', 'UPA_end')]
  }
  if (!is.null(pn)) p$PA=paste0(pn, 1:nrow(p))
  rownames(p)=NULL
  return(p)
}

### counts tables: id1/sample/count
.getPAcounts<-function(p, pn='') {
  dat=p@counts
  if (is.data.frame(dat)) dat=as(as.matrix(dat), 'dgCMatrix')
  if (is.matrix(dat)) dat=as(dat, 'dgCMatrix')
  rid=paste0(pn, 1:nrow(dat))
  cid=colnames(dat)

  df <- .asDf(summary(dat))
  colnames(df)=c('i','j','count')
  df$id1 <- rid[df$i]
  df$sample <-cid[df$j]
  df=transform(df, i=NULL, j=NULL)
  #df=reshape2::melt(d, value.name='count', variable.name='sample', id.vars='PA')
  return(df)
}


#' Merge multiple PACdatasets
#'
#' mergePACds groups nearby PACs from single/multiple PACdataset objects, with or without reference PACds.
#' If first reduces all pA regions by d-nt to a single GRanges, and then count each sample in the merged ranges.
#' For reducing GRanges, two ways are provided: with refPACds or without.
#' This function is suitable for both bulk and single-cell data.
#' This function is particularly useful for grouping nearby cleavage sites into PACs.
#' It is also useful When you have multiple PA or PAC files, each file is from one sample.
#' Then you need to merge these PACds into one PACds for DE or other analyses.
#' But after grouping and/or merging, you may need call annotatePAC to annotate the merged PACs by a GFF annotation.
#' @param PACdsList a PACdataset, or a list of multiple PACdataset objects. The PACds@anno should have columns chr/strand/coord.
#' If there is no colData in PACds, then the sample label will be set as groupN.
#' If PACdsList is a PACdataset, then will treat it as PA and group nearby PAs into PACs.
#' @param d distance to group nearby PACds, default is 24 nt.
#' @param by a charater of "coord" or other string. If coord then use the PA's coord for merging,
#' otherwise use UPA_start and UPA_end in the anno slot for merging PAs.
#' @param refPACds a reference PACds for merging PACdsList in a smarter way.
#' Providing refPACds for merging is useful when there are multiple large pA lists to be merged, which can prevent generating pAs with a very wide range.
#' If there is reference pAs from 3'seq, it is recommended to use it.
#' Meanwhile, it is also recommended to call \code{\link{buildRefPACdsAnno}} with both refPACds and PACdsList (set min.cells/min.counts/max.width) to get a high-confident pA list as reference.
#' However, if there is no refPACds from 3'seq, it is still encouraged to call \code{\link{buildRefPACdsAnno}} with PACdsList alone (and set either PACds as the ref) to obtain high-quality pAs as ref.
#' @return A merged PACdataset. The counts slot stores counts of merged samples.
#' If sample names from different PACdataset objects are any duplicated, then the sample name will be added a suffix .i for each sample in PACds[[i]].
#' The colData slot stores the merge sample annotation from the first column of each @colData.
#' The anno slot contains these columns: chr, strand, coord, UPA_start, UPA_end. Note: Three columns in previous version of movAPA (or current function mergePACds_v0) -- tottag, nPA, maxtag columns -- are not output here anymore.
#' @examples
#' ## make example pacds
#' ds1=makeExamplePACds(seed=123)
#' ## change a bit for ds2
#' ds2=ds1
#' ds2@anno$coord[1:3]=ds2@anno$coord[1:3]-1
#' ds2@anno$UPA_start[1:3]=ds2@anno$UPA_start[1:3]-10
#' ds2@anno$UPA_end[1:3]=ds2@anno$UPA_end[1:3]+5
#'
#' newc=ds2@anno[1:3, ]$coord
#' oldc=ds1@anno[1:3, ]$coord
#'
#' ### merge two pacds by coord
#' ## without using refPACds
#' p1=mergePACds(list(ds1, ds2), d=0, by='coord')
#'
#' ## merge with refPACds ds2
#' ## ps uses ds2 as ref,
#' ## so the final output will use ds2's info if a merged pA can from both ds1 and ds2
#' p2=mergePACds(PACdsList=list(ds1, ds2), refPACds=ds2, d=0, by='coord')
#'
#' ## use ds1 as ref
#' p3=mergePACds(PACdsList=list(ds1, ds2), refPACds=ds1, d=0, by='coord')
#'
#' summary(p1)
#' summary(p2)
#' summary(p3)
#'
#' ## the width of the final PA ranges
#' summary(p1@anno$UPA_end-p1@anno$UPA_start)
#' summary(p2@anno$UPA_end-p2@anno$UPA_start)
#' summary(p3@anno$UPA_end-p3@anno$UPA_start)
#'
#' ## number of reads will not change after merging
#' sum(ds1@counts)+sum(ds2@counts)
#' sum(p1@counts)
#' sum(p2@counts)
#' sum(p3@counts)
#'
#' ## all TRUE, because p2 uses ds2 as reference
#' newc %in% p2@anno$coord
#' oldc %in% p2@anno$coord
#' ## all TRUE, because p3 uses ds2 as reference
#' oldc %in% p3@anno$coord
#' newc %in% p3@anno$coord
#'
#' #### merge two pacds by ranges
#' p1=mergePACds(list(ds1, ds2), d=0, by='range')
#' p2=mergePACds(PACdsList=list(ds1, ds2), refPACds=ds2, d=0, by='range')
#' p3=mergePACds(PACdsList=list(ds1, ds2), refPACds=ds1, d=0, by='range')
#'
#' newc=ds2@anno[1:3, ]$UPA_start
#' oldc=ds1@anno[1:3, ]$UPA_start
#'
#' ## all TRUE, because p2 uses ds2 as reference
#' newc %in% p2@anno$UPA_start
#' oldc %in% p2@anno$UPA_start
#' ## all TRUE, because p3 uses ds2 as reference
#' oldc %in% p3@anno$UPA_start
#' newc %in% p3@anno$UPA_start
#'
#' @name mergePACds
#' @family PACdataset functions
#' @seealso \code{\link{mergePACds_v0}} old version of mergePACds; \code{\link{buildRefPACdsAnno}} to build a reference PACds; rbind() to combine multiple PACdatasets of the same format.
#' @export
mergePACds <- function (PACdsList, d=24, by='coord', refPACds=NULL) {

  .isRefPACds<-function(refPACds) {
    return('source' %in% colnames(refPACds@anno))
  }

  useRef=FALSE

  if (!is.null(refPACds)) {
    if (!.isRefPACds(refPACds)) {
      cat("Warning: it seems that refPACds is not from buildRefPACdsAnno, it is recommend call buildRefPACdsAnno(refPACds, PACdsList, ...) first to build a more abundant reference!\n")
    }
    useRef=TRUE
  }

  bycoord=ifelse(tolower(by)=='coord', TRUE, FALSE)

  if (inherits(PACdsList, "PACdataset")) PACdsList=list(PACdsList)

  .checkPACdsBy(PACdsList, by=bycoord)

  for (i in 1:length(PACdsList)) {
    ## save original sample groups
    if (nrow(PACdsList[[i]]@colData)==0) {
      PACdsList[[i]]@colData=.asDf(matrix( rep(paste0('group',i),
                                                       ncol(PACdsList[[i]]@counts)), ncol=1,
                                                   dimnames =list(colnames(PACdsList[[i]]@counts),'group') ))
    }
  }

  ### format colData to avoid duplicated sample names
  if (length(PACdsList)>1) {
    sn=unlist(lapply(PACdsList, function(ds) return(rownames(ds@colData))))
    snu=unique(sn)
    if (length(sn)!=length(snu)) {
      cat(sprintf("mergePACds: there are %d duplicated sample names in the PACdsList, will add .N to sample names of each PACds\n", length(sn)-length(snu)))
      #sufx=c('x','y','z',letters[1:23])
      for (i in 1:length(PACdsList)) {
        colnames(PACdsList[[i]]@counts)=paste0(colnames(PACdsList[[i]]@counts), '.', i)
        rownames(PACdsList[[i]]@colData)=paste0(rownames(PACdsList[[i]]@colData), '.', i)
      }
    }
  }

  #colnames used for restore order of sample columns
  cnames=unlist(lapply(PACdsList, function(d) return(colnames(d@counts))))

  ### colData
  sng=unlist(lapply(PACdsList, function(ds) return(as.character(ds@colData[, 1]))))
  sn=unlist(lapply(PACdsList, function(ds) return(rownames(ds@colData)) ))
  colData=data.frame(group=sng)
  rownames(colData)=sn


  allAnno=.getPAranges(PACdsList[[1]], bycoord=bycoord, pn='P1_')

  if (length(PACdsList)>=2) {
    for (j in 2:length(PACdsList)) {
      pa2=.getPAranges(PACdsList[[j]], bycoord=bycoord, pn=paste0('P',j,'_'))
      allAnno=rbind(allAnno, pa2)
    }
  }
  invisible(gc())

  cat(sprintf("mergePACds: total %d redundant PACs from %d PACds to merge\n", nrow(allAnno), length(PACdsList)))

  allAnno <- with(allAnno, GenomicRanges::GRanges(seqnames = chr,
                            ranges =IRanges::IRanges(start=UPA_start,
                                            end=UPA_end),
                            strand = strand, id1=PA) )


  # expand 3'end, resize width as d+1
  if (d>0) allAnno=GenomicRanges::resize(allAnno, width=d+IRanges::width(allAnno),
                                         fix="start", use.names=TRUE, ignore.strand=FALSE)

  ### merge Annos ###
  if (useRef) { ## smarter way to merge annos with refPACds
    ## use PACdsList and ref to get merged Anno
    mergedAnnoRef=buildRefPACdsAnno(refPACds=refPACds, PACdsList=PACdsList, by=by, d=d,
                                 min.counts = NULL, min.smps=NULL, max.width=NULL,
                                 verbose=FALSE)
    mergedAnno <- with(mergedAnnoRef@anno, GenomicRanges::GRanges(seqnames = chr,
                                                          ranges =IRanges::IRanges(start=merged_start, ###!!! not UPA_start
                                                                                   end=merged_end),
                                                          strand = strand, id2=paste0('M', 1:length(mergedAnnoRef))) )

    cat(sprintf("mergePACds with refPACds: %d merged PACs (d=%dnt)\n",
                length(mergedAnno), d))

  } else { ## merge all annos without refPACds
    mergedAnno=GenomicRanges::reduce(allAnno, drop.empty.ranges=TRUE)
    mcols(mergedAnno)$id2=paste0('M', 1:length(mergedAnno))
    cat(sprintf("mergePACds without refPACds: %d separate PACs reduce to %d PACs (d=%dnt)\n",
                length(allAnno), length(mergedAnno), d))
  }


  ### note: after buildRef, some range or coord of pAs will be priorily set as ref.
  ## which means that some pAs in `allAnno` cannot be overlapped with ref anymore (because the region may be narrowed by refPACds).
  ## so after findOverlaps, the not-over-lapping ones should be related to the nearest ref.

  ### map oldIDs -- mergedIDs
  #cat("mergePACds: Mapping old PAids to merged PAids\n")
  ov = GenomicRanges::findOverlaps(allAnno, mergedAnno,
                    maxgap=-1L, minoverlap=1L,
                    type=c("any"), select='first',
                    ignore.strand=FALSE)

  ## some pAs are lost, should recover those to nearese Merge-id
  if (sum(is.na(ov))>0) {
    cat(sprintf('mergePACds: %d pAs (total %d) cannot be mapped back to the merged pA ranges, which may because the refPACds is used.\n',
                sum(is.na(ov)), length(ov)))
    allAnno.lost=allAnno[ is.na(ov) ]
    ov.lost = GenomicRanges::findOverlaps(allAnno.lost, mergedAnno,
                                     maxgap=5000L, minoverlap=0L,
                                     type=c("any"), select='first',
                                     ignore.strand=FALSE)
    if (sum(is.na(ov.lost))>0) {
      cat(sprintf('mergePACds: %d pAs from %d unmapped pAs still cannot be mapped to nearest merged pAs (within 5000 bp); these pAs are deleted!\n',
                  sum(is.na(ov.lost)), length(ov.lost)))
    } else {
      cat(sprintf('mergePACds: all these %d pAs can be successfully mapped to nearest merged pAs (within 1000 bp)!\n',
                  length(ov.lost)))
    }
    ov[is.na(ov)]=ov.lost
  }


  ## restore mergedAnno to smaller UPA_start/end, after have done findoverlaps
  if (useRef) { ## smarter way to merge annos with refPACds
    mergedAnno <- with(mergedAnnoRef@anno, GenomicRanges::GRanges(seqnames = chr,
                                                               ranges =IRanges::IRanges(start=UPA_start, ###!!! back to UPA_start
                                                                                        end=UPA_end),
                                                               strand = strand, id2=paste0('M', 1:length(mergedAnnoRef))) )
  }

  ov=data.frame(id1=1:length(ov), id2=ov)
  ov=ov[!is.na(ov$id2), ]

  ov$id1=mcols(allAnno)$id1[ov$id1]
  ov$id2=mcols(mergedAnno)$id2[ov$id2]

  allcnt=.getPAcounts(PACdsList[[1]],  pn='P1_')

  if (length(PACdsList)>=2) {
    for (j in 2:length(PACdsList)) {
      pa2=.getPAcounts(PACdsList[[j]], pn=paste0('P',j,'_'))
      allcnt=rbind(allcnt, pa2)
    }
  }
  invisible(gc())
  cat(sprintf('mergePACds: melted all counts tables, total %d triplet rows\n', nrow(allcnt)))

  PACdsList=NULL

  ### add new id to allcnt: id1/count/sample/id2
  cat(sprintf('mergePACds: link %d old PA IDs to %d new PA IDs by merge\n',
              length(unique(ov$id1)), length(unique(ov$id2))))
  allcnt=data.table::data.table(allcnt)
  ov=data.table::data.table(ov)
  allcnt=merge(allcnt, ov, by='id1')

  #group by id2+sample to get sum counts (not need!)
  #allcnt=allcnt %>% dplyr::group_by(id2, sample)  %>% dplyr::summarise(count=sum(count))

  cat(sprintf('mergePACds: convert %d triplets to dgCMatrix\n', nrow(allcnt)))
  ### change [i,j,value] to id-sample dgCMatrix
  rid=factor(allcnt$id2)
  cid=factor(allcnt$sample)
  allcnt = transform(allcnt,
                     id2 = rid,
                     sample = cid)
  allcnt = Matrix::sparseMatrix(as.integer(allcnt$id2), as.integer(allcnt$sample), x = allcnt$count)
  colnames(allcnt) = levels(cid)
  rownames(allcnt) = levels(rid)
  invisible(gc())
  cat(sprintf('mergePACds: construct Matrix[PA, sample], [%d, %d]\n', nrow(allcnt), ncol(allcnt)))

  if (ncol(allcnt)<=50) allcnt=as.matrix(allcnt) # use standard matrix

  ### reorder mergedAnno and allcnt
  #head(mergedAnno)
  #head(allcnt)
  ## shrink the resize-PA-range, narrow the 3'end
  if (d>0 & !useRef) mergedAnno=GenomicRanges::resize(mergedAnno, width=IRanges::width(mergedAnno)-d,
                                            fix="start", use.names=TRUE, ignore.strand=FALSE)

  mergedAnno=.asDf(mergedAnno)
  colnames(mergedAnno)=c('chr','UPA_start','UPA_end','width','strand','PAid')

  ## Add coord
  mergedAnno$coord=mergedAnno$UPA_end
  mergedAnno$coord[mergedAnno$strand=='-']=mergedAnno$UPA_start[mergedAnno$strand=='-']
  rownames(mergedAnno)=mergedAnno$PAid
  mergedAnno=mergedAnno[, c('chr','strand','coord','UPA_start','UPA_end')]

  mergedAnno$chr=as.character(mergedAnno$chr)
  mergedAnno$strand=as.character(mergedAnno$strand)

  id1=sort(rownames(mergedAnno))
  id2=sort(rownames(allcnt))
  if (length(id1)!=length(id2) | !identical(id1, id2)) {
    cat(sprintf("Warning: in mergePACds, the number of rows (PAs) or rownames after merging is not the same between anno (%d) and counts (%d)\nThis may be because refPACds is used. ", length(id1), length(id2)))
    cmid=base::intersect(id1, id2)
    cat(sprintf("Will use common PAs between anno and counts (%d PAs)\n", length(cmid)))
    mergedAnno=mergedAnno[cmid, ]
    allcnt=allcnt[cmid, ]
  } else {
    allcnt=allcnt[rownames(mergedAnno), ]
  }

  # a few samples (cells) may be 0 after merging (dont know why, but it happens!)
  if(!all(cnames %in% colnames(allcnt))) {
    notin=which(!(cnames %in% colnames(allcnt)))
    warning(sprintf("mergePACds: %d samples are removed after merging\n", length(notin)))
    cnames=cnames[-notin]
    colData=colData[cnames, , drop=F]
  }

  # restore colname order
  if(all(cnames %in% colnames(allcnt))) {
    allcnt=allcnt[, cnames, drop=F]
    colData=colData[cnames, , drop=F]
  }

  d=createPACdataset(counts=allcnt, anno=mergedAnno, colData=colData)
  return(d)
}



#' Build a reference annotation of PACdataset
#'
#' buildRefPACdsAnno builds a reference annotation of PACdataset which can be used for `mergePACdsByRef`.
#' First high-quality pAs meet min.counts/min.smps/max.width are filtered from PACdsList.
#' Then PACdsList and the refPACds are merged into a single PACds, based on `d` nt and by=coord/range.
#' For each merged pA, if it is combined from multiple sources,
#' then the final pA range is determined by priority: source=ref >> max.counts >> min.width >> most distant one.
#'
#' @param refPACds a reference PACds, usually is a PACds with high-confidence, such as polyA sites from 3'seq data.
#' refPACds is used as the basic reference to include more high-confident pAs from other sources provided by PACdsList.
#' @param PACdsList default is NULL, or can be a PACdataset, or a list of multiple PACdataset objects. The PACds@anno should have columns chr/strand/coord.
#' If it is NULL, then will only reduce refPACds by `d` nt. min.counts, min.smps, and max.width are used only when PACdsList is not NULL.
#' The name of the list would be used for the `source` column in the final built refPACds@anno, which represent the source of that pA.
#' @param d distance to group nearby PACds, default is 24 nt.
#' @param by a charater of "coord" or other string. If coord then use the PA's coord for merging,
#' otherwise use UPA_start and UPA_end in the anno slot for merging PAs.
#' @param min.counts pAs with total tagnum >=min.counts will be filtered from PACdsList for building the reference.
#' It can be set as NULL to not filter by counts.
#' @param min.smps pAs expressed in >=min.smps samples (columns of counts) will be filtered from PACdsList for building the reference.
#' It can be set as NULL to not filter by number of samples.
#' @param max.width pAs spanning a region <=max.width will be filtered from PACdsList for building the reference.
#' It can be set as NULL to not filter by width of pA.
#' If all the above three parametes are set NULL, then all pAs in PACdsList will be used for building reference.
#' This is useful when just want to merge all pAs' ranges in a smarter way than calling \code{\link{mergePACds}}.
#' @param verbose TRUE to show message on screen.
#' @return A PACdataset representing the reference pAs.
#' The counts slot stores counts of the finally used samples, which is just for reference.
#' The colData slot is simply set as `group`, which is meanningless.
#' The anno slot contains these columns: chr, strand, coord, UPA_start, UPA_end, source, counts, merged_start, merged_end.
#' The `source` column is the data source of the reference pA, `ref` for the refPACds, and others for other sources in PACdsList.
#' The merged_start and merged_end are the bigger range of the pA that may be combined from multiple sources,
#' while the UPA_start and UPA_end are the smaller range only retieved from one source (with ref in the first priority) .
#'
#' @examples
#' ## make example pacds
#' ds1=makeExamplePACds(seed=123)
#' ## only use ds1 to build ref (which means only reduce ranges by d=24nt)
#' ref=buildRefPACdsAnno(refPACds=ds1, by='coord', d=24)
#' table(ref@anno$source)
#'
#' ## change a bit
#' ds2=ds1
#' ds2@anno$coord[1:3]=ds2@anno$coord[1:3]-1000
#' ds2@anno$UPA_start[1:3]=ds2@anno$coord[1:3]-10
#' ds2@anno$UPA_end[1:3]=ds2@anno$coord[1:3]+5
#'
#' ## using refPACds and high-quality pAs from PACdsList to build ref
#' ref=buildRefPACdsAnno(refPACds=ds1, PACdsList=list(ds1=ds1, ds2=ds2),
#'                       by='coord', d=24,
#'                       min.counts = 50, min.smps=3, max.width=100,
#'                       verbose=TRUE)
#' table(ref@anno$source)
#'
#'
#' ## no high-quality from PACdsList, so only pAs in reduced refPACds were output
#' ref2=buildRefPACdsAnno(refPACds=ds2, PACdsList=list(ds1, ds2),
#'                        by='range', d=0,
#'                        min.counts = 5000, min.smps=3, max.width=50,
#'                        verbose=TRUE)
#' table(ref2@anno$source)
#' @name buildRefPACdsAnno
#' @seealso \code{\link{mergePACds}}
#' @export
buildRefPACdsAnno <- function (refPACds, PACdsList=NULL, d=24, by='coord',
                           min.counts=10, min.smps=10, max.width=100, verbose=TRUE) {

  bycoord=ifelse(tolower(by)=='coord', TRUE, FALSE)

  if (!is.null(PACdsList)) {
    if (inherits(PACdsList, "PACdataset")) PACdsList=list(PACdsList)

    if (is.null(names(PACdsList))) names(PACdsList)=paste0('PACds', 1:length(PACdsList))
    if (any(names(PACdsList)=='')) names(PACdsList)[which(names(PACdsList)=='')]=paste0('PACds', which(names(PACdsList)==''))

    .checkPACdsBy(PACdsList, by=bycoord)
  }


  ## filter high-quality PACds
  .getGoodPACds<-function(PACdsList, verbose=TRUE) {

    if (is.null(min.counts) & is.null(min.smps) & is.null(max.width)) return(PACdsList)

    for (i in 1:length(PACdsList)) {
      p=PACdsList[[i]]

      nidc=1:length(p); nidv=1:length(p); nidw=1:length(p)

      if (!is.null(min.counts)) nidc=which(rowSums(p@counts, na.rm = TRUE)>=min.counts)
      if (!is.null(min.smps)) nidv=which(rowSums(p@counts>0, na.rm = TRUE)>=min.smps)
      if (!is.null(max.width)) nidw=which(p@anno$UPA_end-p@anno$UPA_start<=max.width)

      nid=intersect(nidc, nidv)
      nid=intersect(nid, nidw)

      PACdsList[[i]]=p[nid]
      if (verbose & !is.null(min.counts) & !is.null(min.smps) & !is.null(max.width))
        cat(sprintf("%s: total=%d --> good=%d [>=min.counts(%d) & >=min.smps(%d)) & <=max.width(%d)]\n",
              names(PACdsList)[i], length(p), length(nid), min.counts, min.smps, max.width) )
    }
    return(PACdsList)
  }

  if (!is.null(PACdsList)) {
    PACdsList=.getGoodPACds(PACdsList, verbose=verbose)

    nQC=sum(unlist(lapply(PACdsList, length)))

    if (verbose)
      cat(sprintf( "buildRefPACds: %d highQuality PACs (filtered by min.counts=%s, min.smps=%s, max.width=%s)\n",
                   nQC, min.counts,  min.smps, max.width) )
  }

  ## get a combined GR for all
  ## chr/strand/start/end + mcols=counts, source
  .getGRwithNC<-function(PACdsList, bycoord=TRUE) {
    gr=NULL
    for (i in 1:length(PACdsList)) {
      pacds=PACdsList[[i]]
      source=names(PACdsList)[i]
      if (length(pacds)==0) next

      if (bycoord) {
        pacds@anno$UPA_start=pacds@anno$coord
        pacds@anno$UPA_end=pacds@anno$coord
      }
      gr1 <- GenomicRanges::GRanges(seqnames = pacds@anno$chr,
                                 ranges =IRanges::IRanges(start=pacds@anno$UPA_start,
                                                 end=pacds@anno$UPA_end),
                                 strand = pacds@anno$strand,
                                 source=source, counts=rowSums(pacds@counts, na.rm=TRUE))
      if (is.null(gr))
        gr=gr1
      else
        gr=c(gr, gr1)
    }
    return(gr)
  }

  # seqnames        ranges strand |      source    counts
  #  Chr2   89566-89615      - |         ref       272
  if (!is.null(PACdsList)) {
    allAnno=.getGRwithNC(  c(list(ref=refPACds), PACdsList), bycoord )
  } else {
    allAnno=.getGRwithNC(  list(ref=refPACds), bycoord )
  }

  nref=length(refPACds)

  ## reduce allGR by d
  # expand 3'end, resize width as d+1
  if (d>0) allAnno=GenomicRanges::resize(allAnno, width=d+IRanges::width(allAnno),
                                         fix="start", use.names=TRUE, ignore.strand=FALSE)

  mergedAnno=GenomicRanges::reduce(allAnno, drop.empty.ranges=TRUE)

  if (verbose) {
    if (!is.null(PACdsList)) {
      cat(sprintf("buildRefPACds: %d [%d ref + %d highQuality PACs] reduce to %d PACs (d=%dnt)\n",
                  length(allAnno), nref, nQC, length(mergedAnno), d))
    } else {
      cat(sprintf("buildRefPACds: %d [ref] reduce to %d PACs (d=%dnt)\n",
                  length(allAnno), length(mergedAnno), d))
    }
  }


  ### map allAnno -- mergedAnno, to refine final refs
  ov = GenomicRanges::findOverlaps(allAnno, mergedAnno,
                    maxgap=-1L, minoverlap=1L,
                    type=c("any"), select='first',
                    ignore.strand=FALSE)

  ov=data.frame(ridS=1:length(ov), ridM=ov)
  colnames(ov)=c('ridS','ridM') #single, merge
  ov=ov[!is.na(ov$ridM), , drop=F]


  ## restore the PA-range, narrow the 3'end
  if (d>0) allAnno=GenomicRanges::resize(allAnno, width=IRanges::width(allAnno)-d,
                                         fix="start", use.names=TRUE, ignore.strand=FALSE)
  allAnno=.asDf(allAnno)
  allAnno$ridS=1:nrow(allAnno)
  allAnno=data.table::as.data.table(allAnno)
  ov=data.table::as.data.table(ov)

  n1=nrow(allAnno)

  allAnno=merge(allAnno, ov)
  n2=nrow(allAnno)

  # ridS seqnames  start    end width strand source counts ridM
  #  1     Chr2  89566  89615    50      -    ref    272  106
  if (n1!=n2) cat( sprintf("Warnning in buildRefPACds: (%d) pAs (total %d) were lost after reduce pAs\n", n1-n2, n1))

  ## add merged ranges to allAnno
  mergedAnno=.asDf(mergedAnno)
  mergedAnno$ridM=1:nrow(mergedAnno)
  mergedAnno$merged_start=mergedAnno$start
  mergedAnno$merged_end=mergedAnno$end
  mergedAnno=data.table::as.data.table(mergedAnno[, c('ridM', 'merged_start', 'merged_end')])
  allAnno=merge(allAnno, mergedAnno, by='ridM')

  allAnno$source=as.factor(allAnno$source)
  allAnno$source=relevel(allAnno$source, "ref")

  allAnno$pos=allAnno$end
  allAnno$pos[allAnno$strand=='-']=(-1)*allAnno$end[allAnno$strand=='-']

  ## order of a merged PA's sources by priority
  allAnno=allAnno[order(allAnno$ridM, allAnno$source, -allAnno$counts, allAnno$width, -allAnno$pos, decreasing = FALSE), ]

  ## get the first row of the same Mid
  sid=(1:nrow(allAnno))[!duplicated(allAnno$ridM)]
  allAnno=allAnno[sid, , drop=F]

  allAnno=allAnno[, c('seqnames','strand', 'start','end', 'width','source','counts', 'merged_start', 'merged_end')]
  allAnno$coord=allAnno$end
  allAnno$coord[allAnno$strand=='-']=allAnno$start[allAnno$strand=='-']
  colnames(allAnno)[c(1,3, 4)]=c('chr', 'UPA_start','UPA_end')
  allAnno$chr=as.character(allAnno$chr)
  allAnno$strand=as.character(allAnno$strand)
  allAnno$source=as.character(allAnno$source)

  if (verbose)
    cat( sprintf("buildRefPACds: %d reference PACs returned\n", nrow(allAnno)) )

  d=createPACdataset(counts=matrix(allAnno$counts, ncol=1), anno=allAnno, verbose=FALSE)
  return(d)
}


#' Count overlapping between PACdatasets
#'
#' countOvpPACds counts overlapping between PACdatasets.
#'
#' This function is used to compare PA coordinates of two PACdataset objects.
#' @param pacds1 a PACdataset.
#' @param pacds2 another PACdataset.
#' @param mode1 use coord or UPA_start/UPA_end for pacds1, default is point (using coord). If mode1=region, then use UPA_start/UPA_end.
#' @param mode2 same as mode1 except for pacds2.
#' @param d distance to count overlapping, default is 0 (no gap between qry and sbj).
#' @param pats1 cutoffs to filter pacds1 by total tag number of all samples. If pats1<=1, then do not filter.
#' Multiple cutoffs can be set, like pats1=c(10, 20, 30).
#' @param pats2 same as pats1 but to filter pacds2.
#' @return A dataframe with 8 columns: PAT1, PAT2, Total1, Total2, Ovp1, Ovp2, Ovp1Pct, Ovp2Pct.
#' Row number is the combination number of pats1 and pats2.
#' @examples
#' data(PACds)
#' pacds1=PACds
#' pacds2=PACds
#' pacds2@anno$coord=pacds2@anno$coord+sample(-10:10, 1)
#' ## Count exactly the same positions between two pacds without filtering
#' countOvpPACds(pacds1, pacds2, d=0)
#' ## count overlapping with different cutoffs
#' countOvpPACds(pacds1, pacds2, d=50, pats1=1, pats2=seq(0,20,10))
#' ## count overlapping with different cutoffs, comparing ds1's coord ~ ds2's region
#' countOvpPACds(pacds1, pacds2, mode1='point', mode2='region', d=50, pats1=1, pats2=seq(0,20,10))
#' @name countOvpPACds
#' @family PACdataset functions
#' \code{\link{findOvpPACds}} to get an overlapped PACdataset.
#' @export
countOvpPACds <- function (pacds1, pacds2, mode1='point', mode2='point', d=0, pats1=1, pats2=1) {


  mode1=tolower(mode1)
  mode2=tolower(mode2)
  if (!(mode1 %in% c('point', 'region'))) stop("mode1 must be point or region")
  if (!(mode2 %in% c('point', 'region'))) stop("mode2 must be point or region")

  if (mode1=='point') {
  gr1 <- with(pacds1@anno, GRanges(seqnames = chr,
                              ranges =IRanges::IRanges(start=coord,
                                              end=coord),
                              strand = strand) )
  } else if (mode1=='region') {
    gr1 <- with(pacds1@anno, GRanges(seqnames = chr,
                                       ranges =IRanges::IRanges(start=UPA_start,
                                                       end=UPA_end),
                                       strand = strand) )
  }

  if (mode2=='point') {
    gr2 <- with(pacds2@anno, GRanges(seqnames = chr,
                                     ranges =IRanges::IRanges(start=coord,
                                                     end=coord),
                                     strand = strand) )
  } else if (mode2=='region') {
    gr2 <- with(pacds2@anno, GRanges(seqnames = chr,
                                     ranges =IRanges::IRanges(start=UPA_start,
                                                     end=UPA_end),
                                     strand = strand) )
  }

  if (is.null(pats1)) pats1=1
  if (is.null(pats2)) pats2=1

  if (!(length(pats1)==1 & pats1[1]==1)) tot1=rowSums(pacds1@counts)
  if (!(length(pats2)==1 & pats2[1]==1)) tot2=rowSums(pacds2@counts)

  ostat=matrix(data=NA, nrow=0, ncol=6, dimnames = list(c(NULL),c('PAT1','PAT2','Total1','Total2','Ovp1','Ovp2')))

  for (p1 in pats1) {
    gr11=gr1
    if (p1>1) {
      idx1=which(tot1>=p1)
      gr11=gr1[idx1, ]
    }

    for (p2 in pats2) {
      gr22=gr2
      if (p2>1) {
        idx2=which(tot2>=p2)
        gr22=gr2[idx2, ]
      }

      ov = GenomicRanges::findOverlaps(gr11, gr22,
                        maxgap=d-1, minoverlap=0L,
                        type=c("any"), select='all',
                        ignore.strand=FALSE)
      ov=.asDf(ov)
      n1=length(gr11)
      n2=length(gr22)
      ov1=length(unique(ov$queryHits))
      ov2=length(unique(ov$subjectHits))

      ostat=rbind(ostat, c(p1, p2, n1, n2, ov1, ov2))
    }
  }

  ostat=.asDf(ostat)
  ostat$Ovp1Pct=sprintf("%d%s",round(ostat$Ovp1/ostat$Total1*100), '%')
  ostat$Ovp2Pct=sprintf("%d%s",round(ostat$Ovp2/ostat$Total2*100), '%')
  return(ostat)
}


#' Find overlapping PACdataset
#'
#' findOvpPACds returns overlapping PACdataset.
#'
#' This function is used to subset one PACds with another PACds.
#' @param qryPACds a query PACdataset, will subset and return the qryPACds.
#' @param sbjPACds a subject PACdataset.
#' @param qryMode use coord or UPA_start/UPA_end for qryPACds. Default is point (using coord), if qryMode=region, then use UPA_start/UPA_end.
#' @param sbjMode same as qryMode except for sbjPACds.
#' @param d distance to count overlapping, default is 0 (no gap between qry and sbj).
#' @param qryPAT cutoff to filter qryPACds. If sbjPAT<=1, then do not filter.
#' @param sbjPAT same as qryPAT but to filter sbjPACds.
#' @param returnOvp default is TRUE; if TRUE then return overlapping qryPACds.
#' @param returnNonOvp default is FALSE; if TRUE then return nonOverlapping qryPACds.
#' If both returnOvp=T and returnNonOvp=TRUE, then return a list(ovp, nonOvp)
#' @param verbose print information.
#' @return A subset PACdataset of qryPACds.
#' @examples
#' data(PACds)
#' pacds1=PACds
#' pacds2=PACds
#' pacds2@anno$coord=pacds2@anno$coord+sample(-10:10, 1)
#' ## return overlapped qryPACds without filtering, comparing coord~coord
#' pd=findOvpPACds(qryPACds=pacds1, sbjPACds=pacds2, d=0, qryPAT=1, sbjPAT=1)
#' ## return overlapped qryPACds with filtering, comparing coord~coord
#' pd=findOvpPACds(qryPACds=pacds1, sbjPACds=pacds2, d=50, qryPAT=1, sbjPAT=10)
#' ## return overlapped qryPACds comparing qry's coord ~ sbj's UPA_start/end
#' pd=findOvpPACds(qryPACds=pacds1, sbjPACds=pacds2, d=0, qryMode='point', sbjMode='region')
#' @name findOvpPACds
#' @family PACdataset functions
#' \code{\link{countOvpPACds}} to count overlapping between two PACdatasets.
#' @export
findOvpPACds <- function (qryPACds, sbjPACds, qryMode='point',sbjMode='point', d=0, qryPAT=1, sbjPAT=1, returnOvp=TRUE, returnNonOvp=FALSE, verbose=TRUE) {

  qryMode=tolower(qryMode)
  sbjMode=tolower(sbjMode)
  if (!(qryMode %in% c('point', 'region'))) stop("qryMode must be point or region")
  if (!(sbjMode %in% c('point', 'region'))) stop("sbjMode must be point or region")
  if (!returnNonOvp & !returnOvp) stop("At least returnNonOvp=T or returnOvp=T!")

  if(qryPAT>1) {
    n1=nrow(qryPACds@counts)
    qryPACds=subsetPACds(qryPACds, totPACtag = qryPAT)
    if (verbose) {
      n2=nrow(qryPACds@counts)
      cat(sprintf("qryPACds after (qryPAT>=%d): %d --> %d", qryPAT, n1, n2), '\n')
    }
  }

  if(sbjPAT>1) {
    n1=nrow(sbjPACds@counts)
    sbjPACds=subsetPACds(sbjPACds, totPACtag = sbjPAT)
    if (verbose) {
      n2=nrow(sbjPACds@counts)
      cat(sprintf("sbjPACds after (sbjPAT>=%d): %d --> %d", sbjPAT, n1, n2), '\n')
    }
  }

  if (qryMode=='point') {
    gr1 <- with(qryPACds@anno, GRanges(seqnames = chr,
                                     ranges =IRanges::IRanges(start=coord,
                                                     end=coord),
                                     strand = strand) )
  } else if (qryMode=='region') {
    gr1 <- with(qryPACds@anno, GRanges(seqnames = chr,
                                       ranges =IRanges::IRanges(start=UPA_start,
                                                       end=UPA_end),
                                       strand = strand) )
  }

  if (sbjMode=='point') {
    gr2 <- with(sbjPACds@anno, GRanges(seqnames = chr,
                                     ranges =IRanges::IRanges(start=coord,
                                                     end=coord),
                                     strand = strand) )
  } else if (sbjMode=='region') {
    gr2 <- with(sbjPACds@anno, GRanges(seqnames = chr,
                                       ranges =IRanges::IRanges(start=UPA_start,
                                                       end=UPA_end),
                                       strand = strand) )
  }

  ov = GenomicRanges::findOverlaps(gr1, gr2,
                     maxgap=d-1, minoverlap=0L,
                     type=c("any"), select='all',
                     ignore.strand=FALSE)
  ov=.asDf(ov)

  qryIdx=unique(ov$queryHits)

  n1=length(gr1)

  if (returnNonOvp) {
    qryPACdsNot=qryPACds
    qryPACdsNot@counts=qryPACds@counts[-qryIdx, , drop=F]
    qryPACdsNot@anno=qryPACds@anno[-qryIdx, , drop=F]
    if (verbose) {
      n2=nrow(qryPACdsNot@counts)
      cat(sprintf("qryPACds NotOverlapping: %d --> %d", n1, n2), '\n')
    }
  }

  if (returnOvp) {
    qryPACds@counts=qryPACds@counts[qryIdx, , drop=F]
    qryPACds@anno=qryPACds@anno[qryIdx, , drop=F]
    if (verbose) {
      n2=nrow(qryPACds@counts)
      cat(sprintf("qryPACds overlapping: %d --> %d", n1, n2), '\n')
    }
  }

  if (returnOvp & returnNonOvp) {
    return(list(ovp=qryPACds, nonOvp=qryPACdsNot))
  } else if (returnOvp) {
    return(qryPACds)
  } else if (returnNonOvp) {
    return(qryPACdsNot)
  }

}

#' Annotate a PACdataset with known PACdatasets
#'
#' annotateByKnownPAC returns overlapping status and min distance with given PACdatasets.
#'
#' This function is used to validate a PACdataset with other known PAC datasets and calculates the min distance between two PACs.
#' @param pacds a query PACdataset.
#' @param knownPACdss a PACdataset or a list of PACdatasets.
#' @param labels a character vector storing column names added to pacds.
#' Each PACdataset in knownPACdss will add two columns ({labels[i]}_ovp, _dist). The length of labels must be 1 or the same as knowPACdss (if it is a list).
#' @param d Distance (nt) allowed from known PACs to pacds.
#' @param verbose print information.
#' @return A PACdataset with columns ({labels[i]}_ovp, _dist) added.
#'
#' The _ovp column is 0/1 denoting whether a PAC in pacds is overlapping with knowPACdss.
#'
#' The _dist column is NA (no overlapping) or a integer denoting the min distance between a PAC in pacds and the overlapped PAC in knownPACdss.
#' @examples
#' data(PACds); pacds=PACds
#' knownPACdss=list(pacds[1:500], pacds[500:1000])
#' p2=annotateByKnownPAC(pacds, knownPACdss, labels=c('test1','test2'))
#' summary(p2@anno$test2_dist)
#' @name annotateByKnownPAC
#' @family PACdataset functions
#' \code{\link{findOvpPACds}} to find overlapping between two PACdatasets.
#' @export
annotateByKnownPAC <- function (pacds, knownPACdss, labels, d=50, verbose=TRUE) {

  #library(GenomicRanges, verbose = FALSE)


  gr1 <- with(pacds@anno, GRanges(seqnames = chr,
                                     ranges =IRanges::IRanges(start=coord,
                                                     end=coord),
                                     strand = strand) )

  if (inherits(knownPACdss, "PACdataset")) {
    if (length(labels)!=1) stop("knownPACdss is PACdataset, but more than one elements in labels\n")
    knownPACdss=list(knownPACdss)
    names(knownPACdss)=labels
  } else if (!is.list(knownPACdss)) {
    stop("knownPACdss is neither a PACdataset nor a list!\n")
  } else if (is.list(knownPACdss)) {
    if (length(labels)!=length(knownPACdss)) stop("Not equal element# in knownPACdss and labels\n")
    names(knownPACdss)=labels
  }

 for (i in 1:length(knownPACdss)) {
   gr2 <- with(knownPACdss[[i]]@anno, GRanges(seqnames = chr,
                                      ranges =IRanges::IRanges(start=coord,
                                                      end=coord),
                                      strand = strand) )
   lbl=names(knownPACdss)[i]

   ov = GenomicRanges::findOverlaps(gr1, gr2,
                     maxgap=d-1, minoverlap=0L,
                     type=c("any"), select='first',
                     ignore.strand=FALSE)
   isOvp=ov
   im=which(!is.na(ov))
   isOvp[-im]=0
   isOvp[im]=1
   minD=rep(NA, length(ov))
   minD[im]=abs(knownPACdss[[i]]@anno$coord[ov[im]]-pacds@anno$coord[im])

   pacds@anno[, paste0(lbl,'_ovp')]=isOvp
   pacds@anno[, paste0(lbl,'_dist')]=minD

   if (verbose) {
     cat('query:',length(gr1),'\n')
     cat('known:',length(gr2),'(',lbl,')\n')
     cat(length(im), 'query PACs overlapping known PACs\n')
   }
 }
 return(pacds)
}

# -------------- *** Normalization **** ---------------

#' Normalize a PACdataset
#'
#' normalizePACds returns a normalized PACds with @counts normalized.
#'
#' This function is used for normalization with different methods.
#' @param pacds a PACdataset.
#' @param method normalization method, can be TPM/CPM/DEseq/DESeq2/EdgeR.
#'               CPM=TPM (after normalization, libsize of all columns are 1e6), DESeq=DESeq2.
#' @return A normalized PACdataset with pacds@counts normalized.
#' @examples
#' data(PACds)
#' PACdsNorm=normalizePACds(PACds, method='TPM')
#' @name normalizePACds
#' @family PACdataset functions
#' @export
normalizePACds <- function (pacds, method='TPM') {

  method=toupper(method)
  if (!(method %in% c('TPM','CPM','DESEQ','DESEQ2','EDGER','TMM'))) {
    stop("normalizePACds: method not TPM/CPM/DESEQ/DESEQ2/EDGER")
  }

  if (method=='TPM') method='EDGER'

  sf=getSizeFactor(pacds@counts, method=method)
  dat=normBySizeFactor(pacds@counts, sizeFactor=sf)
  if (!.isAnyMatrix(dat)) dat=as.matrix(dat)
  pacds@counts=dat
  return(pacds)
}


## test sizefactors of DESeq2
## my funcs:
# sf=getSizeFactor(pacds@counts, method='deseq2')
# head(normBySizeFactor(pacds@counts, sizeFactor=sf))

## DEseq2
# dds=PACds2DESeqDds(pacds, noIntergenic = FALSE)
# sf2=sizeFactors(estimateSizeFactors(dds)) #same as sf
# dds=estimateSizeFactors(dds)
# head(counts(dds, normalized=TRUE)) #same as above

# getSizeFactor(matrix,'EDGER'/'DESEQ'/'TPM')
# return: size factors or libsizes (colsum/f)
# note: to use PAC/factor for normalization
# examples
# #get size factors
# sf=getSizeFactor(matrix,'EDGER');
# #get libsize
# libSize=getSizeFactor(matrix,'EDGER',oLibsize=T)
# #x% TPM method=TPMx, thd=0.95
# libSize=getSizeFactor(matrix,'TPMx',oLibsize=T,thd=0.95)
getSizeFactor<-function(dat, method='DESEQ', oLibsize=FALSE, thd=0.95) {

  method=toupper(method)
  if (is.data.frame(dat)) dat=as.matrix(dat)
  if (method!='EDGER' && method!='DESEQ' && method!='TPM' && method!='TPMX') {
    method='DESEQ'
  }
  stopifnot (method %in% c('EDGER','DESEQ','DESEQ2','CPM','TPM', 'TPMX'))
  if (method=='EDGER') {
    f <- edgeR::calcNormFactors(dat)
    #f <- f/exp(mean(log(f)))
  } else if (method %in% c('DESEQ','DESEQ2')) {
    ds=DESeq2::DESeqDataSetFromMatrix(countData = dat,
                              colData = S4Vectors::DataFrame(1:ncol(dat)),
                              design = ~ 1)
    f=DESeq2::sizeFactors(DESeq2::estimateSizeFactors(ds))
  } else if (method=='TPM') {
    f=colSums(dat)/1000000
  }else if (method=='TPMX') {
    rowsums=rowSums(dat);
    uprow=quantile(rowsums,probs=thd)
    rest=dat[rowsums<uprow,]
    libs=colSums(rest);
    minlib=min(libs)
    f=libs/minlib
  }
  if (oLibsize){
    f=as.integer(colSums(dat)/f)
  }
  return(f)
}

#  normBySizeFactor(matrix, sizeFactor)
#  dat: a matrix
#  useage:
#  sf=getSizeFactor(matrix,'EDGER'); dat.norm=normBySizeFactor(dat,sf)
normBySizeFactor<-function(dat, sizeFactor) {
  tmp=matrix(rep(sizeFactor,nrow(dat)), nrow=nrow(dat), byrow=T)
  return(round(dat/tmp))
}


# -------------- *** Annotation **** ---------------

#' Get actual grams for an abbreviation
#'
#' getVarGrams returns a vector of grams given an abbreviation.
#'
#' This function is used to get common grams.
#' @param grams a character string (e.g., V1, v1, MOUSE, mouse, mm10, mm) or a string vector.
#' V1 means AATAAA and its variants; MOUSE means hexamers of polyA signals in mouse.
#' Mouse signals are obtained from https://polyasite.unibas.ch/atlas#2, which reside in a region of 60 nt upstream to 10 nt downstream of one of the poly(A) sites of a cluster.
#' @return A vector of grams.
#' @examples
#' ## get AATAAA and its variants
#' getVarGrams('v1')
#' ## get mouse hexamers
#' getVarGrams('mm')
#' ## not doing anything
#' getVarGrams(c('AATAAA','AATT'))
#' @name getVarGrams
#' @family APA signal functions
#' @export
getVarGrams<-function(grams) {
  if (!is.null(grams)) {
    grams=toupper(grams)
    if (length(grams)==1 & grams[1]=='V1')
      grams=c('AATAAA','TATAAA','CATAAA','GATAAA','ATTAAA','ACTAAA','AGTAAA','AAAAAA','AACAAA','AAGAAA','AATTAA','AATCAA','AATGAA','AATATA','AATACA','AATAGA','AATAAT','AATAAC','AATAAG')
    if (length(grams)==1 & (grams[1] %in% c('MM','MOUSE','MM10'))) {
      grams=c('AATAAA','ATTAAA','TATAAA','AGTAAA','AATACA','CATAAA','AATATA','GATAAA','AATGAA','AATAAT','AAGAAA','ACTAAA','AATAGA','ATTACA','AACAAA','ATTATA','AACAAG','AATAAG')
    }
  }
  return(grams)
}

#' Annotate a PACdataset with polyA signals
#'
#' annotateByPAS returns distance-to-PAC of give pattern(s) within a given range of PAC.
#'
#' This function is used to get polyA signals around PACs and calculates the min distance between the signal to PAC.
#' @param pacds a query PACdataset.
#' @param bsgenome chrmosome fasta files, an object of BSgenome or FaFile, see faFromPACds().
#' @param grams a character vector to specify a gram like AATAAA, or v1 (AATAAA's variants), or multiple grams. grams can be not equal length, like c('AATAAA','AAATTT','CCCT')
#' @param from to specify the range near PACs, PAC is the 0 position.
#' e.g., from=-50, to=-1 to subset 50 nt (PAC is the 0 or 51st position, upstream 50nt of PAC), see faFromPACds().
#' @param to similar to from.
#' @param label a character to specify output column name. pacds will be added one or two columns ({label}_gram, _dist).
#' If only one element in grams, then label could be NULL, then label=gram. Only {label}_dist will be ouput.
#' If multiple elements in grams, then label should be specified. Then two columns ({label}_gram, _dist) will be added.
#' The _gram column gives the gram that is closest to the PAC.
#' _dist is the start position of a gram to a PAC.
#' @param priority a numeric vector to set the priority and subgroups of grams if grams has multiple elements, default is NULL.
#'               For example, if grams=c('AATAAA','ATTAAA','AAAAAA','TTTAT), priority=(1,2,3,3), then will first search for AATAAA, if not exists, then for ATTAAA, then the remaining AAAAAA/TTTAT.
#'               If priority=NULL, then will treat all elements in grams as the same group (no priority).
#' @param chrCheck if TRUE, then all chr in PACds should be in bsgenome, otherwise will ignore those non-consistent chr rows in PACds.
#' @return A PACdataset with columns (label_gram (if multiple grams), _dist) added.
#'
#' The _gram column is NA (no signal) or the gram closest to the PAC.
#'
#' The _dist column is NA (no signal) or a integer denoting the min distance between a PAC and grams.
#'
#' For example, given * is the PAC, then ...AATAAA*, dist=6; ...*AATAAA, dist=1.
#' @examples
#' \dontrun{
#' ## First, load the reference genome sequences that are already represented as a BSgenome object.
#' library("BSgenome.Oryza.ENSEMBL.IRGSP1")
#' bsgenome <- BSgenome.Oryza.ENSEMBL.IRGSP1
#' data(PACds); pacds=PACds
#' ## scan AATAAA upstream 50bp of PACs
#' test=annotateByPAS(pacds, bsgenome, grams='AATAAA', from=-50, to=-1, label=NULL)
#' summary(test@anno$AATAAA_dist)
#' ## scan AATAAA's 1nt variants
#' test=annotateByPAS(pacds, bsgenome, grams='V1', from=-50, to=-1, label=NULL)
#' table(test@anno$V1_gram)
#' ## scan custom grams
#' test=annotateByPAS(pacds, bsgenome, grams=c('AATAAA','ATTAAA','GATAAA','AAAA'),
#'                                      from=-50, to=-1, label='GRAM')
#' table(test@anno$GRAM_gram)
#' ## scan with priority, scan AATAAA first,
#' ## if no, then scan ATTAAA, if no, then scan the remaining grams.
#' test2=annotateByPAS(pacds, bsgenome, grams=c('AATAAA','ATTAAA','GATAAA','AAAA'),
#'                                       priority=c(1,2,3,3), from=-50, to=-1, label='GRAM')
#' table(test2@anno$GRAM_gram)
#' ## only scan AATAAA, the number will be the same as test2's AATAAA
#' test3=annotateByPAS(pacds, bsgenome, grams=c('AATAAA'),
#'                                       priority=NULL, from=-50, to=-1, label='GRAM')
#' sum(!is.na(test3@anno$GRAM_dist))
#' }
#' @name annotateByPAS
#' @family PACdataset functions
#' @family APA signal functions
#' @export
annotateByPAS <- function (pacds, bsgenome, grams, from, to, priority=NULL, label=NULL, chrCheck=TRUE) {

  grams=toupper(grams)

  if (length(grams)!=1) {
    if (is.null(label)) stop("grams has >1 elements, label should not be NULL!\n")
  }

  if (length(grams)==1) {
    if (!is.null(priority)) stop("grams has only 1 element, priority should be NULL!\n")
  }

  if (is.null(priority)) priority=rep(1, length(grams))
  if (length(priority)!=length(grams)) stop("priority must be the same length as grams!\n")

  if (is.null(label)) label=grams

  grams=getVarGrams(grams)

  if (!isChrConsistent(pacds, bsgenome, allin=TRUE)) {
    cat('PACds chr not all in seqnames of bsgenome\n')
    if (chrCheck) {
      stop('Please check chr names!\n')
    } else {
      cat("chrCheck=FALSE, Remove rows in PACds with chr name not in bsgenome\n")
      pacds=subsetPACds(pacds, chrs=chrnames(bsgenome), verbose=TRUE)
    }
  }

  #seq=faFromPACds(pacds[1:5], bsgenome, what='updn', fapre=NULL, up=-5, dn=5)
  seq=faFromPACds(pacds, bsgenome, what='updn', fapre=NULL, up=from, dn=to)

  if (length(seq)!=nrow(pacds@counts)) {
    cat(length(seq),'fasta seqs, but',nrow(pacds@counts),'PACs, not matched!\n')
    stop("This is probably because that the coord/from/to of pacds is out of range of chromosomes\n")
  }

  #grams=c('GCGA','TATG','CC')
  dd=list()
  for (gram in grams) {
    vm=Biostrings::vmatchPattern(pattern=gram, subject=seq, fixed=TRUE)
    #CCCTTCCTTTATAACTAGTGTCGCAACAATAAAATTTGAGCTTTGATCA
    # vm$start=28, 27chars before AATAAA

    nr=elementNROWS(vm)
    hasGram=nr
    hasGram[nr!=0]=1
    idx=which(hasGram!=0)
    svm=lapply(vm[idx], function(par) {Biostrings::start(par)})

    #convert to 0-base position
    svm2=lapply(svm, function(par) {
      min(abs(par+from-1))
      })

    gramDist=rep(NA, length(hasGram))
    gramDist[idx]=unlist(svm2)

    if (length(grams)==1) {
      pacds@anno[, paste0(label, '_dist')]=gramDist
    } else {
      dd[[gram]]=gramDist
    }
  }

  if (length(grams)==1) return(pacds)

  ## dd is a list, like dd$AATAAA=[NA, 5, NA, NA,... nseq]

  # if there are multiple grams, get the gram with min distance to PAC
  dd2=as.matrix(.asDf(dd)) #row=PAs, col=grams
  pacds@anno[, paste0(label, '_gram')]=NA
  pacds@anno[, paste0(label, '_dist')]=NA

  levels=sort(unique(priority))
  for (l in levels) {

    # deal with the first priority grams
    dd=dd2[, priority==l, drop=F]

    # get gram of the min dist
    mi=apply(dd, 1, function(par) {
      if (sum(!is.na(par))==0) return(NA)
      m=min(par, na.rm=TRUE)
      mi=which(par==m)[1]
      return(mi)
    })

    # get min dist
    mind=apply(dd, 1, function(par) {
      if (sum(!is.na(par))==0) return(NA)
      m=min(par, na.rm=TRUE)
      return(m)
    })

    minGram=mi
    minGram[!is.na(mi)]=grams[priority==l][mi[!is.na(mi)]]

    idNA=is.na(pacds@anno[, paste0(label, '_gram')]) # if no prior gram was found, then fill the PAS
    pacds@anno[idNA, paste0(label, '_gram')]=minGram[idNA]
    pacds@anno[idNA, paste0(label, '_dist')]=mind[idNA]
  }

  return(pacds)
}



# If a PACds is annotated.
# @param PACds a PACdataset
# return TRUE or FALSE
isPACdsAnnotated <- function (PACds) {
  if (sum(!(c('chr','strand','gene','ftr','coord','ftr_start','ftr_end') %in% colnames(PACds@anno)))!=0) {
    #stop("chr,strand,gene,ftr,coord,ftr_start,ftr_end not all in PACds@anno")
    return(FALSE)
  }
  return(TRUE)
}

#' Annotate PACs
#'
#' annotatePAC annotates PACs with a given GFF annotation.
#'
#' If after annotation, the PAC number is changed, it will raise a warning but not error.
#' In such case, you may need to check aGFF.
#' PAs in gene named `character(0)` will also be removed.
#'
#' @param pac a data frame with [chr/coord/strand] or a PACdataset.
#' @param aGFF specify a genome annotation, see parseGenomeAnnotation().
#' @param verbose TRUE to show message.
#' @return  A PACdataset with annotation if pac is PACdataset, or a data frame if pac is a data frame.
#'          If a PACdataset is returned, the original @anno columns in pac are remained (duplicated annotation columns are removed);
#'          duplicated rows with the same chr/strand/coord will be removed; will add @supp$stopCodon to store the stopCodon of all transcripts, which can be used for ext3UTRPACds().
#'
#'   The following columns are added for annotation:
#'   "ftr", "ftr_start","ftr_end","gene","biotype",","gene_start","gene_end","gene_stop_codon",
#'   "upstream_id","upstream_start","upstream_end","downstream_id","downstram_start",
#'   "downstream_end","three_UTR_length","three_extend".
#'   \itemize{
#'  \item{ftr}: the type of feature about coord of PAC, including 3UTR, 5UTR, CDS, intron, exon, and intergenic;
#'  \item{ftr_start}: start position of the "ftr";
#'  \item{ftr_end}: end position of the "ftr";
#'  \item{gene}: the gene name of the feature;
#'  \item{gene_start}: start position of the gene;
#'  \item{gene_end}: end position of the gene;
#'  \item{gene_type}: the classification of the gene, such as protein_coding,
#'    long non-conding RNA (lncRNA), non-coding RNA (ncRNA), tRNA and so on;
#'  \item{gene_stop_codon}: the end position of stop codon for protein coding gene or transcript,
#'    is the end postion of gene for ncRNA;
#'  \item{upstream_id}: the upstream gene name of poly(A) sites that are located in intergenic;
#'  \item{upstream_start}: the start position of `upstream_id`;
#'  \item{upstream_end}: the end position of `upstream_id`;
#'  \item{downstream_id} same as upstream_id expected for downstream gene;
#'  \item{downstream_start} see upstream_start;
#'  \item{downstream_end} see upstream_end;
#'  \item{three_UTR_length}: the length of 3'UTR, which is equal to
#'    poly(A) sites minus the end position of stop codon of gene. This is only for poly(A) site
#'    located in 3'UTR or intergenic;
#'  \item{three_extend}: used to identify the poly(A) site in 3'UTR extension, which is equal to
#'     poly(A) sites minus end positon of its upstream gene. This is only for poly(A) site
#'     located in intergenic;
#'  }
#' @examples
#' \dontrun{
#' load('Oryza_sativa.IRGSP-1.0.42.gff.rda')
#' data(PACds)
#' ## Because the demo data already contain the annotation,
#' ## here we removed the annotation columns first.
#' PACds1=PACds
#' PACds1@anno[,c('gene','ftr','gene_type','ftr_start','ftr_end')]=NULL
#' PACds1=annotatePAC(PACds1, gff)
#' ## annotate PAC data with a gff3 file
#' newpac=annotatePAC(pac, aGFF='mm10.gff3')
#' ## Annotate a PACdataset
#' newpac=annotatePAC(pacds, aGFF=TxDb.Mmusculus.UCSC.mm10.ensGene)
#' ## Annotate a PAC data frame with existing gff.rda from parseGenomeAnnotation()
#' newpac=annotatePAC(pac, 'mm10.rda')
#' }
#'
#' @name annotatePAC
#' @family PACdataset functions
#' @seealso \code{\link{parseGenomeAnnotation}} to get an gff annotation object.
#' @export
annotatePAC <- function(pac, aGFF, verbose=FALSE){

  pd=NULL
  if (inherits(pac, "PACdataset")) {
    if (!identical(rownames(pac@anno),rownames(pac@counts))) {
      rownames(pac@anno)=paste0('PA',1:nrow(pac@anno))
      rownames(pac@counts)=rownames(pac@anno)
    }
    pd=pac
    pac=pac@anno[,c('chr','strand','coord')]
  }

  N=nrow(pac)

  #Loading PAC information
  pac$chr <- as.character(pac$chr)
  pac$coord<- as.integer(pac$coord)
  pac$strand<- as.character(pac$strand)

  g=parseGenomeAnnotation(aGFF)
  anno.need=g$anno.need
  anno.rna=g$anno.rna
  anno.frame=g$anno.frame
  rm(g); invisible(gc())

  gffchr=unique(anno.need$seqnames)
  pacchr=unique(pac$chr)
  if (!isChrConsistent(pacchr, gffchr, allin=FALSE)) stop(cat("Chr names of GFF and PAC not the same, please check!\n"))


  #Overlap
  #find \exon\cds\3utr\5utr\intergenic\exon\intron region
  gr.pac <- GRanges(seqnames =as.character(pac$chr) ,
                    ranges =IRanges::IRanges(start=as.integer(pac$coord) ,
                                    end=as.integer(pac$coord)),
                    strand =as.character(pac$strand))

  gr.exon <- GRanges(seqnames =as.character(anno.need$seqnames) ,
                     ranges =IRanges::IRanges(start=as.integer(anno.need$start) ,
                                     end=as.integer(anno.need$end)),
                     strand =as.character(anno.need$strand))

  ov <- GenomicRanges::findOverlaps(gr.pac, gr.exon)
  ov <- .asDf(ov)
  colnames(ov)<-c("pacID","annoID")
  ov$type <- as.character(anno.need$type[ov$annoID])
  ov$gene_id <- as.character(anno.need$gene_id[ov$annoID])
  ov$biotype <- as.character(anno.need$biotype[ov$annoID])
  ov$strand <- as.character(anno.need$strand[ov$annoID])
  ov$start <- as.integer(anno.need$start[ov$annoID])
  ov$end <- as.integer(anno.need$end[ov$annoID])
  ov$labe <- paste0(ov$biotype,"_",ov$type)
  #print("###all RNA biotype and type information")
  #print(table(ov$labe))
  level.names <- unique(ov$labe)
  high.levels <- c("protein_coding_three_prime_UTR","NA_three_prime_UTR",
                   "protein_coding_CDS","nontranslating_CDS_exon", "NA_CDS",
                   "protein_coding_intron","nontranslating_CDS_intron","NA_intron",
                   "protein_coding_five_prime_UTR" ,"NA_five_prime_UTR",
                   "protein_coding_exon", "NA_exon",
                   "lncRNA_exon" , "ncRNA_exon",
                   "lncRNA_intron","ncRNA_intron")
  add.levels <- base::setdiff(level.names,high.levels)
  #levels=level.names[c(4,3,2,6,7,1,9,5,8)]
  #add.levels <-c("ncRNA_three_prime_UTR","nothing_intron")

  three.id <- grep("three_prime_UTR$",add.levels,ignore.case = FALSE)
  exon.id <- grep("exon$",add.levels,ignore.case = FALSE)
  intron.id <- grep("intron$",add.levels,ignore.case = FALSE)
  five.id <- grep("five_prime_UTR",add.levels,ignore.case=FALSE)
  id <- c(three.id,intron.id,five.id,exon.id)
  add.levels.first <- add.levels[id]
  add.levels.rest <- add.levels[-id]

  ov$labe<- factor(ov$labe,
                   levels=c(high.levels,add.levels.first,add.levels.rest))

  ov.order <- ov[order(ov$pacID,ov$labe),]
  ov.unique <- ov.order[!duplicated(ov.order$pacID),]

  #check result
  gr.rna <- GRanges(seqnames =as.character(anno.rna$seqnames) ,
                    ranges =IRanges::IRanges(start=as.integer(anno.rna$start) ,
                                    end=as.integer(anno.rna$end)),
                    strand =as.character(anno.rna$strand))
  ov.rna <- GenomicRanges::findOverlaps(gr.pac,gr.rna)
  if(length(base::setdiff(ov.rna@from,ov.unique$pacID))){
    print("warning this annotation have some problem")
    print(length(base::setdiff(ov.rna@from,ov.unique$pacID)))
  }
  if(length(base::setdiff(ov.unique$pacID,ov.rna@from))){
    print("warning this annotation have some problem")
    print(length(base::setdiff(ov.unique$pacID,ov.rna@from)))
  }
  #ov.unique:for mRNA and lncRNA annotation

  pac$ftr <-  "intergenic"
  pac$biotype <- NA
  pac$ftr_strand <- NA
  pac$ftr_start <- NA
  pac$ftr_end <- NA
  pac$gene_id <- NA
  pac$gene_start <- NA
  pac$gene_end <- NA
  pac$gene_stop_codon <- NA
  pac$ftr[ov.unique$pacID] <- ov.unique$type
  pac$gene_id[ov.unique$pacID]  <- ov.unique$gene_id
  pac$gene_start[ov.unique$pacID]  <- anno.rna$start[match(ov.unique$gene_id,anno.rna$Parent)]
  pac$gene_end[ov.unique$pacID]  <- anno.rna$end[match(ov.unique$gene_id,anno.rna$Parent)]


  pac$biotype[ov.unique$pacID]  <-ov.unique$biotype
  pac$ftr_strand[ov.unique$pacID]  <- ov.unique$strand
  pac$ftr_start[ov.unique$pacID]  <- ov.unique$start
  pac$ftr_end[ov.unique$pacID]  <- ov.unique$end

  three.utr <- anno.frame[grep("three_prime_UTR",anno.frame$type,ignore.case = TRUE),]
  anti.id <- which(three.utr$strand=="-")
  three.utr.anti <- three.utr[anti.id,]
  three.utr.anti <- three.utr.anti[order(three.utr.anti$Parent,three.utr.anti$start,decreasing=TRUE),]
  three.utr.anti <- three.utr.anti[!duplicated(three.utr.anti$Parent),]
  three.utr.anti$stop_codon <- as.integer(three.utr.anti$end+1)

  three.utr.sen <- three.utr[-anti.id,]
  three.utr.sen <- three.utr.sen[order(three.utr.sen$Parent,three.utr.sen$start),]
  three.utr.sen <-  three.utr.sen[!duplicated(three.utr.sen$Parent),]
  three.utr.sen$stop_codon <- as.integer(three.utr.sen$start -1)
  three.utr <- rbind(three.utr.sen,three.utr.anti)

  anno.rna$stop_codon <- anno.rna$end
  anno.rna$stop_codon[which(anno.rna$strand=="-")] <- anno.rna$start[which(anno.rna$strand=="-")]

  index <- match(three.utr$Parent, anno.rna$transcript_id)

  if (sum(is.na(index))>0) stop("3UTR's transcript not in anno.rna$transcript_id")

  anno.rna$stop_codon[index] <- three.utr$stop_codon
  #annotation of intergeneic
  #Upstream
  pac$upstream_id <- NA
  pac$upstream_start <- NA
  pac$upstream_end <- NA
  #pac$upstream_stop_codon <- NA
  #ov.int <- follow(gr.pac,gr.rna)
  ov.int <- follow(gr.pac, gr.rna)
  query.id <-which(!is.na(ov.int))
  inter.id <- which(pac$ftr=="intergenic")
  #length(query.id);length(inter.id);length(base::intersect(inter.id,query.id))
  #unknow.id <- base::setdiff(inter.id,query.id)
  inter.id.sub <- base::intersect(inter.id,query.id)
  subject.id <- ov.int[inter.id.sub]

  pac$gene_id[inter.id.sub] <- as.character(anno.rna$Parent[subject.id])
  pac$biotype[inter.id.sub] <- as.character(anno.rna$biotype[subject.id])
  pac$ftr_strand[inter.id.sub] <-as.character(anno.rna$strand[subject.id])
  pac$ftr_start[inter.id.sub] <- as.integer(anno.rna$start[subject.id])
  pac$ftr_end[inter.id.sub] <- as.integer(anno.rna$end[subject.id])

  pac$gene_start[inter.id.sub] <- as.integer(anno.rna$start[subject.id])
  pac$gene_end[inter.id.sub] <- as.integer(anno.rna$end[subject.id])

  pac$upstream_id[inter.id.sub] <- as.character(anno.rna$Parent[subject.id])
  pac$upstream_start[inter.id.sub] <- as.integer(anno.rna$start[subject.id])
  pac$upstream_end[inter.id.sub] <- as.integer(anno.rna$end[subject.id])

  index <- match(pac$gene_id,anno.rna$Parent)
  pac$gene_stop_codon <- anno.rna$stop_codon[index]

  #downstream
  pac$downstream_id <- NA
  pac$downstream_start <- NA
  pac$downstream_end <- NA
  #ov.int <- follow(gr.pac,gr.rna)
  ov.int <- precede(gr.pac,gr.rna)
  query.id <-which(!is.na(ov.int))
  inter.id <- which(pac$ftr=="intergenic")
  #length(query.id);length(inter.id);length(base::intersect(inter.id,query.id))
  #unknow.id <- base::setdiff(inter.id,query.id)
  inter.id.sub <- base::intersect(inter.id,query.id)
  subject.id <- ov.int[inter.id.sub]
  pac$downstream_id[inter.id.sub] <- as.character(anno.rna$Parent[subject.id])
  pac$downstream_start[inter.id.sub] <- as.integer(anno.rna$start[subject.id])
  pac$downstream_end[inter.id.sub] <- as.integer(anno.rna$end[subject.id])
  pac$ftr_strand <- NULL
  inter.id <- which(pac$ftr=="intergenic")
  pac$ftr_start[inter.id] <- pac$upstream_end[inter.id]
  pac$ftr_end[inter.id] <- pac$downstream_start[inter.id]
  neg.id <- which(pac$strand=="-")
  inter.neg.id <- base::intersect(inter.id,neg.id)
  pac$ftr_start[inter.neg.id] <- pac$upstream_start[inter.neg.id]
  pac$ftr_end[inter.neg.id] <- pac$downstream_end[inter.neg.id]

  pac$three_UTR_length <- NA
  index <- which(pac$ftr %in% c("three_prime_UTR","intergenic"))
  pac$three_UTR_length[index] <- abs(pac$coord[index]-pac$gene_stop_codon[index])

  pac$three_extend <- NA
  index <- which(!is.na(pac$upstream_id))
  pac$three_extend[index] <- pac$coord[index] - pac$upstream_end[index]
  index.anti <- which(pac$strand=="-")
  index <- base::intersect(index.anti,index)
  pac$three_extend[index] <- pac$upstream_start[index] - pac$coord[index]

  if (nrow(pac)!=N) {
    cat("Warning: PA number changed after annotation",N,'-->',nrow(pac),"\n")
  }

  invisible(gc())

  colnames(pac)[colnames(pac)=='biotype']='gene_type'
  colnames(pac)[colnames(pac)=='gene_id']='gene'
  pac$ftr[pac$ftr=='three_prime_UTR']='3UTR'
  pac$ftr[pac$ftr=='five_prime_UTR']='5UTR'

  if(!is.null(pd)) { #PACds
    #remove duplicates if exists
    dup=duplicated(pd@anno[,c('chr','strand','coord')])
    if (sum(dup)>0) {
      cat("Warning: pac is PACdataset, duplicate rows (same chr/stand/coord) exist in pac\n")
      cat("Warning: will remove duplicates", sum(dup),"rows\n")
      pd=pd[!dup]
    }

    dup=duplicated(pac[,c('chr','strand','coord')])
    if (sum(dup)>0) {
      pac=pac[!dup, ]
    }

    if (nrow(pac)!=nrow(pd@counts)) {
      cat("Warning: PA number changed after annotation and removing duplicates", nrow(pd@counts),'-->',nrow(pac),"\n")
    }

    cm=base::intersect(colnames(pd@anno), colnames(pac))
    cm=cm[!(cm %in% c('chr','strand','coord'))]
    if (length(cm)>0) pd@anno[, cm]=NULL

    pd@anno$XXID=rownames(pd@anno)
    pd@anno=merge(pd@anno, pac, by=c('chr','strand','coord'))
    rownames(pd@anno)=pd@anno$XXID; pd@anno$XXID=NULL
    pd@counts=pd@counts[rownames(pd@anno), , drop=F]

    stopCodon=anno.rna[,c('gene_id','transcript_id','stop_codon')]
    colnames(stopCodon)[colnames(stopCodon)=='gene_id']='gene'
    pd@supp$stopCodon=stopCodon

    # remove character(0) genes
    gid=which(pd@anno$gene=='character(0)')
    if (length(gid)>0) {
      cat("Warning: remove PAs in [character(0)] genes", length(gid), "rows\n")
      pd=pd[-gid, ]
    }

    return(pd)
  }

  return(pac)
}

#' Extend 3'UTR
#'
#' ext3UTRPACds extends annotated 3'UTR by a give length for a PACdataset object.
#'
#' PACs in extended 3UTR region (@anno$ftr=intergenic) are set as 3UTR, then the 3UTR length (toStop) is calculated for all 3'UTR PACs.
#'
#' @param pacds a PACdataset object, which should already be annotated by annotatePAC().
#' @param ext3UTRlen an integer denoting the extended length, like 500bp.
#' @param extFtr set the ftr column of PACs in extended 3UTR as extFtr.
#' @return A PACdataset with @anno$ftr column changed and $toStop column added.
#' PACs that are located in intergenic region but within ext3UTRlen bp of the upstream gene end are considered as 3'UTR PACs.
#' For extended 3'UTR PACs, the columns of ftr_start, ftr_end, gene_start, gene_end, gene_stop_codon are set as NA.
#' Note: 3'UTR extension is defined based on the three_extend column, so even for non-coding genes, there may also be PACs in the extended 3'UTR.
#' @examples
#' data(PACds)
#' ext3UTRPACds(PACds, 500)
#' @name ext3UTRPACds
#' @family PACdataset functions
#' \code{\link{annotatePAC}} to annotate a PACdataset.
#' @export
ext3UTRPACds <- function(pacds, ext3UTRlen, extFtr='3UTR'){
  if (!('three_extend' %in% colnames(pacds@anno))) stop("three_extend not in pacds@anno, please call annotatePAC first!")
  extid=which(pacds@anno$ftr=='intergenic' & pacds@anno$three_extend<=ext3UTRlen)
  cat(length(extid), 'PACs in extended 3UTR (ftr=intergenic >> ftr=3UTR)\n')
  cat('Get 3UTR length (anno@toStop) for 3UTR/extended 3UTR PACs\n')

  pacds@anno$ftr[extid]=extFtr
  pacds@anno$gene[extid]=pacds@anno$upstream_id[extid]
  pacds@anno[extid,c('ftr_start','ftr_end','gene_start','gene_end','gene_stop_codon')]=NA
  pacds@anno$toStop=NA
  pacds@anno$toStop[pacds@anno$ftr=='3UTR']=pacds@anno$three_UTR_length[pacds@anno$ftr=='3UTR']
  if (extFtr!='3UTR')
    pacds@anno$toStop[pacds@anno$ftr==extFtr]=pacds@anno$three_UTR_length[pacds@anno$ftr==extFtr]

  # id=pacds@anno$strand=='+' & pacds@anno$ftr=='3UTR'
  # pacds@anno$toStop[id]=pacds@anno$coord[id]-pacds@anno$ftr_start[id]+1
  # id=pacds@anno$strand=='-' & pacds@anno$ftr=='3UTR'
  # pacds@anno$toStop[id]=pacds@anno$ftr_end[id]-pacds@anno$coord[id]+1

  if (min(pacds@anno$toStop[!is.na(pacds@anno$toStop)])<0) stop("ext3UTRPACds: error ftr=3UTR but toStop (<0), please check coord/ftr_start/ftr_end/three_UTR_length\n")

  return(pacds)
}

#' Test extended 3'UTR length
#'
#' testExt3UTR test the number of 3UTR PA changes after extending a given length.
#'
#' @param pacds a PACdataset object, which should already be annotated by annotatePAC().
#' @param ext3UTRlen an integer vector denoting the extended lengthes, like c(1000, 2000, 3000)
#' @param toPlot whether to plot a barplot for the number changes.
#' @return A data frame with two columns: ext3UTRLen, addedPAnum denoting numbers of recruited 3UTR PAs for each extended length.
#' @examples
#' data(PACds)
#' testExt3UTR(PACds, seq(1000, 10000, by=1000))
#' @name testExt3UTR
#' @family PACdataset functions
#' \code{\link{ext3UTRPACds}} to extend3UTR for a PACdataset.
#' @export
testExt3UTR <- function(pacds, ext3UTRlen=seq(1000, 10000, by=1000), toPlot=TRUE){

  if (!('three_extend' %in% colnames(pacds@anno))) stop("three_extend not in pacds@anno, please call annotatePAC first!")

  nums=unlist(lapply(ext3UTRlen, function(par) sum(pacds@anno$ftr=='intergenic' & pacds@anno$three_extend<=par, na.rm=TRUE)))
  nums=data.frame(ext3UTRLen=ext3UTRlen, addedPAnum=nums)

  if (toPlot) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mfrow=c(1, 2))
    barplot(addedPAnum~ext3UTRLen, data=nums)
    if (nrow(nums)>1) {
      diffnums=diff(nums$addedPAnum)
      diffnums=data.frame(ext3UTRLen=ext3UTRlen[-1], diffPAnum=diffnums)
      barplot(diffPAnum~ext3UTRLen, data=diffnums)
    }
  }

  return(nums)
}


# -------------- *** Data I/O **** ---------------
#' Create a PACdataset from counts data
#'
#' createPACdataset create a PACdataset from counts and anno matrix/data.frame.
#' This function is particularly suitable for reading PACdatasets from single-cell data with sparse counts table and a separate meta data.
#'
#' @param counts PAC expression or ratio data, can be data.frame/matrix/dgCMatrix.
#'        If the matrix is in character mode, then will be converted into numeric mode.
#' @param anno detailed annotation information for PAs, with the same row number as the counts slot.
#' @param colData sample defination. Each row is one experiment, and the column is the sample group.
#' @param supp a list storing additional data, like stopCodon.
#' @param forceSparse whether to force the counts matrix as sparse matrix. By default the function determines whether the matrix is sparse (if it contains more than 50 percentage of 0s and ncol is greater than 10). But set forceSparse=TRUE can skip the sparse check.
#' @param verbose TRUE to print message.
#' @return A \code{\link{PACdataset}} object
#' @name createPACdataset
#' @examples
#' data(PACds)
#' p1=createPACdataset(counts=PACds@counts, anno=PACds@anno)
#' anno=PACds@anno; rownames(anno)=NULL
#' p2=createPACdataset(counts=PACds@counts, anno=anno)
#' p3=createPACdataset(counts=PACds@counts, anno=anno, forceSparse=TRUE)
#' @family PACdataset functions
#' @seealso \code{\link{readPACds}} to read PACdataset from file or data.frame.
#' @export
createPACdataset <- function (counts, anno, colData, supp, forceSparse=F, verbose=TRUE) {

  if (missing(x = counts) || missing(x = anno)) {
    stop("Must provide both 'counts' and 'anno'")
  }

  counts=asAnyMatrix(counts, forceSparse=forceSparse)

  if (!.isAnyMatrix(counts)) stop("error counts type, must be dgCMatrix/matrix/data.frame!")

  if (!is.data.frame(anno) & !is.matrix(anno)) stop("error anno type, must be data.frame/matrix!")
  anno=.asDf(anno)
  if (!all(c('chr','strand','coord') %in% colnames(anno))) stop("chr/strand/coord must be in anno to define PA coordinates!")
  anno$coord=as.numeric(anno$coord)

  if (nrow(counts)!=nrow(anno)) stop("row number of counts and anno must be the same, each row is one PA!")
  if (!identical(rownames(counts), rownames(anno))) {
    if (verbose)
      warning(
      "Not identical rownames present in counts and anno, making identical as 'PAi'",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(counts)=paste0('PA', 1:nrow(counts))
    rownames(anno)=paste0('PA', 1:nrow(counts))
  }

  if (is.null(colnames(counts))) {
    if (verbose)
      warning(
        "colnames of counts is NULL, making as 'sample<i>'",
        call. = FALSE,
        immediate. = TRUE
      )
    colnames(counts)=paste0('sample', 1:ncol(counts))
  }


  if (!missing(x = colData)) {
    if (!is.data.frame(colData) & !is.matrix(colData)) stop("error colData type, must be data.frame/matrix!")
    colData=.asDf(colData)
    if (nrow(colData)!=ncol(counts)) stop("row number of colData should be equal to column number of counts, each row/column is one sample!")
    if (!all(rownames(colData) %in% colnames(counts))) stop("row names of colData should be identical to column names of counts, each row/column is one sample!")
  } else {
    colData=.asDf(matrix( rep('group1', ncol(counts)), ncol=1, dimnames =list(colnames(counts), 'group')) )
  }

  for (i in 1:ncol(colData)) {
    colData[,i]=factor(colData[,i])
  }
  #colData=colData[order(rownames(colData)), , drop=F]

  if (!missing(x = supp)) {
  } else {
    supp=list()
  }

  PACds=new("PACdataset", counts=counts, colData=colData, anno=anno, supp=supp)
  return(PACds)
}


#' Read a PACdataset
#'
#' readPACds reads PAC counts and sample annotation into a PACdataset.
#'
#' @param pacFile a file name or a data frame/matrix/dgCMatrix. If it is a file, it should have header, with at least (chr, strand, coord) columns.
#' This file could have other columns, including gff cols (gene/gene_type/ftr/ftr_start/ftr_end/UPA_start/UPA_end) and user-defined sample columns.
#' If there are at least one non-numeric columns other than above gff cols, then all remaining columns are considered as annotation columns.
#' If all remaining columns are numeric, then they are all treated as sample columns.
#' Use annotatePAC() first if need genome annotation of coordinates.
#' @param colDataFile a file name or a data frame. If it is a file, then it is an annotation file of samples with header,
#' rownames are samples (must be all in pacFile), columns names are sample groups.
#' There could be single or multiple columns to define the groups of samples.
#' When colDataFile=NULL, then readPACds will automately retreive sample columns and gff columns (if any) from pacFile.
#' If there is no sample columns, then will set colData as a data frame with 1 column (=group) and 1 row (=tag), and element=group1.
#' If pacfile or colDataFile is a character, then it is a file name, so readPACds will read data from file.
#' @param noIntergenic TRUE/FALSE. If TRUE, then will remove PACs in intergenic (ftr='^inter')
#' @param PAname specify how to set the name (rowname) of PACs.
#' PAname=PA, the PA name is set as 'gene:PAN'; PAname=coord, then 'gene:coord'.
#' @return A PACdataset object, with @anno being a data frame with at least three columns chr/strand/coord.
#' If there is no sample column, then will add one sample named tag in group1.
#' @examples
#' data(PACds)
#' ## read simple PACfile that only has columns chr/strand/coord
#' pacFile=PACds@anno[,c('chr','strand','coord')]
#' colDataFile=NULL
#' p=readPACds(pacFile, colDataFile)

#' ## read PACfile that has columns chr/strand/coord and sample columns
#' pacFile=PACds@anno[,c('chr','strand','coord')]
#' pacFile=cbind(pacFile, PACds@counts[,c('anther1','embryo1','anther2')])
#' colDataFile=NULL
#' p=readPACds(pacFile, colDataFile)
#' p@colData; head(p@counts)

#' ## read PACfile that has columns chr/strand/coord, sample columns, and gff cols
#' ## like gene/gene_type/ftr/ftr_start/ftr_end/UPA_start/UPA_end
#' pacFile=PACds@anno
#' pacFile=cbind(pacFile, PACds@counts[,c('anther1','embryo1','anther2')])
#' colDataFile=NULL
#' p=readPACds(pacFile, colDataFile)
#' p@colData; head(p@counts); head(p@anno)

## read from data frame of PACfile and colDataFile
#' pacFile=PACds@anno
#' smps=c('anther1','embryo1','anther2')
#' pacFile=cbind(pacFile, PACds@counts[,smps])
#' colDataFile=.asDf(matrix(c('group1','group2','group1'),
#'              ncol=1, dimnames=list(smps, 'group')))
#' p=readPACds(pacFile, colDataFile)
#' p@colData; head(p@counts); head(p@anno)

## read from file names of PACfile and colDataFile
#' write.table(pacFile, file='pacFile', row.names=FALSE)
#' write.table(colDataFile, file='colDataFile', row.names=TRUE)
#' p=readPACds(pacFile='pacFile',
#'             colDataFile='colDataFile', noIntergenic=TRUE, PAname='PA')
#'
#' @name readPACds
#' @family PACdataset functions
#' @seealso \code{\link{createPACdataset}} to create PACdataset from counts and anno tables.
#' @export
readPACds<-function(pacFile, colDataFile=NULL, noIntergenic=TRUE, PAname='PA') {

  if (!(PAname %in% c('PA','coord'))) {
    stop("PAname must be PA or coord")
  }

  if (inherits(pacFile, what='character')) {
    d=read.table(pacFile, header=T)
  } else {
    d=.asDf(pacFile)
    if (!inherits(d, what='data.frame'))
      stop("Unknown type of pacFile, should be filename/data.frame/matrix/dgCMatrix!")
  }

  gffcols=c('gene','gene_type','ftr','ftr_start','ftr_end','UPA_start','UPA_end')

  if (sum(!(c('chr','strand','coord') %in% colnames(d)))!=0) {
    stop("chr,strand,coord not all in header of pacfile")
  }

  d$coord=as.numeric(d$coord)

  cat(nrow(d),'PACs\n')

  if ('ftr' %in% colnames(d)) {
    if (noIntergenic) {
      d=d[getNonItgFtrId(d$ftr),]
      cat(nrow(d),'No intergenic PACs\n')
    }
    #WB's data 20190620
    d$ftr[d$ftr=='three_prime_UTR']='3UTR'
    d$ftr[d$ftr=='five_prime_UTR']='5UTR'
  }

  #d[d=='unkown']=NA

  if ('ftr_start' %in% colnames(d)) {
    if (!is.numeric(d$ftr_start)) {d$ftr_start=as.numeric(d$ftr_start)}
    if (!is.numeric(d$ftr_end)) {d$ftr_end=as.numeric(d$ftr_end)}
  }

  if ('gene' %in% colnames(d)) {
    idx=which(d$gene=='NULL')
    if (length(idx)>0) {
      cat(length(idx),'gene name is NULL, change to chrStrand')
      d$gene[idx]=paste0(d$chr[idx],d$strand[idx])
    }

    #order 5' to 3'
    d1=d[d$strand=='+',]
    d1=d1[order(d1$gene,d1$coord,decreasing = FALSE),]
    d2=d[d$strand=='-',]
    d2=d2[order(d2$gene,d2$coord,decreasing = TRUE),]
    d=rbind(d1,d2)

    if (PAname=='coord') {
      paid=paste0(d$gene,':',d$coord)
    } else if (PAname=='PA') {
      rg=rle(d$gene)
      paid=paste0(d$gene,':PA',unlist(lapply(rg$lengths,seq)))
    }
  } else {
    paid=paste0('PA',1:nrow(d))
  }

  if (anyDuplicated(paid)) {
    dupsID=paid[duplicated(paid)]
    cat(sprintf("%d PA IDs are duplicated (probably because one gene in both strands), please remove duplicated rows before readPACds!\n", length(dupsID)))
    cat("Duplicated IDs in genes:", toString(paid[duplicated(paid)]), '\n')
    stop("readPACds error: duplicated rows!\n")
  }

  rownames(d)=paid

  allcols=colnames(d)


  # group defination of sample columns
  if (!is.null(colDataFile)) {
    if (is.vector(colDataFile)) {
      colData=read.table(colDataFile, colClasses="character")
    } else {
      colData=colDataFile
    }

    if (sum(rownames(colData) %in% colnames(d))!=nrow(colData)) {
      stop("rownames of annofile not all in columns of pacfile")
    }

  } else {
    #remove gffcolid and chr/strand/coord, other columns are sample columns, and the sample group is 'group', groupname is 'group1'
    smpcols=which(allcols %in% c(gffcols,'chr','strand','coord'))
    smpcols=allcols[-smpcols]
    #no tag columns, then add tag=1 columns
    if (length(smpcols)==0) {
      smpcols='tag'
      d$tag=1
    } else { #one columns is chr, then they are all annotations
      for (i in smpcols) {
        if (!(is.numeric(d[,i]))) {
          smpcols='tag'
          d$tag=1
          break
        }
      }
    }
    colData=.asDf(matrix( rep('group1',length(smpcols)), ncol=1, dimnames =list(smpcols,'group') ))
  }

  for (i in 1:ncol(colData)) {
    colData[,i]=factor(colData[,i])
  }
  colData=colData[order(rownames(colData)), , drop=F]
  anno=d[,-which(colnames(d) %in% rownames(colData))]
  anno[anno=='unkown']=NA

  counts=d[, rownames(colData), drop=F]
  PACds=createPACdataset(counts=counts, colData=colData, anno=anno)
  return(PACds)
}


#' Write a PACdataset into a file
#'
#' writePACds writes a PACdataset object into a file.
#'
#' @param PACds a PACdataset object.
#' @param file a file name.
#' @param colDataFile a file name. If not NULL, then output PACds@colData into this file.
#' @return NULL.
#' @examples
#' data(PACds)
#' writePACds(PACds, file='pacds.txt', colDataFile='coldata.txt')
#' readPACds(pacFile='pacds.txt', colDataFile='coldata.txt')
#' @name writePACds
#' @family PACdataset functions
#' @export
writePACds<-function(PACds, file, colDataFile=NULL) {
  d=cbind(PACds@anno, PACds@counts)
  write.table(d, file=file, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
  if (!is.null(colDataFile)) write.table(PACds@colData, file=colDataFile)
}

#' Read a DESeqDataset
#'
#' readDESeqDs reads file of PACs or genes into a DESeqDataset.
#'
#' @param datafile a gene/PA file.
#' @param annofile annotation for samples.
#' @param filetype specify the file type, gene/PA.
#' @param noIntergenic valid when filetype=PA, specifying whether to remove intergenic PACs.
#' @param PAname valid when filetype=PA, specifying how to set the name (rowname) of PACs.
#' PAname=PA, the PA name is set as 'gene.PAN';
#' PAname=coord, then 'gene.coord'.
#' @return A DESeqDataset with each row a gene or PAC.
#' @examples
#' data(PACds)
#' writePACds(PACds, file='pacds.txt', colDataFile='coldata.txt')
#'
#' ## Read PACs to DESeqDataset
#' dds=readDESeqDs(datafile='pacds.txt', annofile='coldata.txt',
#'      filetype='PA',noIntergenic=TRUE, PAname='PA')
#' head(mcols(dds)); colData(dds); head(assays(dds)$counts)
#' @name readDESeqDs
#' @export
readDESeqDs <- function (datafile, annofile, filetype, noIntergenic=TRUE, PAname=NULL) {
  #library(DESeq2, verbose = FALSE)
  if (!is.null(PAname)) {
    if (!(PAname %in% c('PA','coord'))) {
      stop("PAname must be PA or coord")
    }
  }

  if (filetype=='PA') {
    ds=readPACds(datafile, annofile, noIntergenic=TRUE, PAname='PA')

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = ds@counts,
                                  colData = ds@colData,
                                  design = ~ 1)
    mcols(dds) <- S4Vectors::DataFrame(mcols(dds), ds@anno)
    return(dds)
  }


  d=read.table(datafile,header=T)
  if (sum(!(c('gene') %in% colnames(d)))!=0) {
    stop("gene not in header of pacfile")
  }

  cat(nrow(d),'rows\n')

  if (noIntergenic) {
    id=grep('igt$',d$gene)
    if (length(id)>0) d=d[-id,]
    id=which(d$gene=="NULL")
    if (length(id)>0) d=d[-id,]
    cat(nrow(d),'No intergenic (xx.igt) or NULL genes\n')
  }

  rownames(d)=d$gene

  colData=read.table(annofile, colClasses="character")
  if (sum(rownames(colData) %in% colnames(d))!=nrow(colData)) {
    stop("rownames of annofile not all in columns of pacfile")
  }
  for (i in 1:ncol(colData)) {
    colData[,i]=factor(colData[,i])
  }
  colData=colData[order(rownames(colData)),]
  anno=d[,-which(colnames(d) %in% rownames(colData))]
  cts=d[,rownames(colData)]

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                colData = colData,
                                design = ~ 1)
  mcols(dds) <- S4Vectors::DataFrame(mcols(dds), anno)
  return(dds)
}

# Collapse PACds to geneDs (DESeqDataSet)
# @param pacds A PACdataset, with @anno should include "gene" column.
# @return A DESeqDataSet.
# Note: Will not output other anno columns like chr/strand/gene_type. (Because DESeq2 do not allow other columns in the anno) (TODO)
PACds2geneDs <- function(pacds, noIntergenic=TRUE) {
  #library(DESeq2, verbose = FALSE)

  if (noIntergenic) pacds=subsetPACds(pacds, noIntergenic=TRUE)
  d=cbind(gene=pacds@anno$gene, .asDf(pacds@counts))
  conds=colnames(pacds@counts)

  #geneds=.asDf(d %>% group_by(gene) %>% summarise_all(funs(sum)))
  geneds <- aggregate(do.call(cbind, mget(colnames(d)[-1]))~gene, d, sum)
  rownames(geneds)=geneds$gene
  geneds$gene=NULL

  # if (sum(c('chr','strand','gene_type') %in% colnames(pacds@anno))==3) {
  #   anno=pacds@anno[pacds@anno$gene %in% geneds$gene, c('gene','chr','strand','gene_type')]
  #   anno=unique(anno)
  #   rownames(anno)=anno$gene
  #   anno$gene<-NULL
  #   anno=anno[geneds$gene,]
  #   geneds$gene=NULL
  # }

  geneds=DESeq2::DESeqDataSetFromMatrix(countData = geneds,
                                colData = pacds@colData,
                                design = ~ 1)
  #Otherwise DESeq2 will error
  mcols(geneds)=NULL
  #<- DataFrame(mcols(geneds), anno)
  return(geneds)
}

# Covert pacds to DESeqDds, just change the format.
# Each PA will be considered as a gene.
PACds2DESeqDds<-function(pacds, noIntergenic=FALSE) {
  #library(DESeq2, verbose = FALSE)

  if (noIntergenic) pacds=subsetPACds(pacds, noIntergenic=noIntergenic)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = pacds@counts,
                                colData = pacds@colData,
                                design = ~ 1)
  #mcols(dds) <- DataFrame(mcols(dds), pacds@anno)
  mcols(dds) <- NULL
  return(dds)
}


# ---------------- *** DESeq2 LRT two factors *** ------------------
# Get factor and levels from ddsRd
# If factor1.. levels1 are provided, then will validate first.
getAndValidateFactorNamesFromDdsRd<-function (ddsRd, afactor1=NULL, afactor2=NULL, aref1=NULL, aref2=NULL) {
  infos=list()
  #Get refrerence level and factor name
  coefNames=DESeq2::resultsNames(ddsRd)
  factorNames=colnames(colData(ddsRd)) #len pos

  vsid=grep('_vs_',coefNames)
  refNames=strsplit(coefNames[vsid],split='_vs_')
  refNames=unique(unlist(lapply(refNames,'[',2))) #dn 0140

  ref1=refNames[1]; ref2=refNames[2]
  for (f in colnames(colData(ddsRd))) {
    if (ref1 %in% colData(ddsRd)[,f]) {
      factor1=f #pos
    }
    if (ref2 %in% colData(ddsRd)[,f]) {
      factor2=f #len
    }
  }

  levels1=levels(colData(ddsRd)[,factor1])
  levels2=levels(colData(ddsRd)[,factor2])
  infos[['factor1']]=factor1
  infos[['factor2']]=factor2
  infos[['levels1']]=levels1
  infos[['levels2']]=levels2

  if (ref1!=levels1[1] | ref2!=levels2[1]) {
    stop(cat('ref1!=levels[1]\n'))
  }
  cat('ref1=',ref1,'\n')
  cat('ref2=',ref2,'\n')

  if (!is.null(afactor1)) {
    if (factor1!=afactor1) {
      stop(cat('afactor1=',afactor1, 'not equal',factor1,'\n'))
    }
  }
  if (!is.null(afactor2)) {
    if (factor2!=afactor2) {
      stop(cat('afactor2=',afactor2, 'not equal',factor2,'\n'))
    }
  }
  if (!is.null(aref1)) {
    if (ref1!=aref1) {
      stop(cat('aref1=',aref1, 'not equal',ref1,'\n'))
    }
  }
  if (!is.null(aref2)) {
    if (ref2!=aref2) {
      stop(cat('aref2=',aref2, 'not equal',ref2,'\n'))
    }
  }
  return(infos)
}


# Get all information from a ddsRd like results(condB~condA@posX)
# @param ddsRd A reduced model from LRT
# @param afactor1..aref1.. For validation
# @return
# Except for betas_xxx is a list, other each element is a results(DDS),
# denoting the change of two levels (conditions) under the level of another factor (group).
# > names(rl) ref=0140 dn
# [1] "posmd_vs_dn.len0140"   "posmd_vs_dn.len0350"   "posmd_vs_dn.len0700"   "posmd_vs_dn.len1200"   "posup_vs_dn.len0140"
# [6] "posup_vs_dn.len0350"   "posup_vs_dn.len0700"   "posup_vs_dn.len1200"   "len0350_vs_0140.posdn" "len0350_vs_0140.posmd"
# [11] "len0350_vs_0140.posup" "len0700_vs_0140.posdn" "len0700_vs_0140.posmd" "len0700_vs_0140.posup" "len1200_vs_0140.posdn"
# [16] "len1200_vs_0140.posmd" "len1200_vs_0140.posup" "betas_posdn.len0140"
# > colnames(rl[['betas_posdn.len0140']]$betas) -- change change matrix
# [1] "pos_md_vs_dn.len0140"        "pos_up_vs_dn.len0140"        "len_0350_vs_0140.posdn"      "len_0700_vs_0140.posdn"
# [5] "len_1200_vs_0140.posdn"      "posmd_vs_dn.len0350_vs_0140" "posup_vs_dn.len0350_vs_0140" "posmd_vs_dn.len0700_vs_0140"
# [9] "posup_vs_dn.len0700_vs_0140" "posmd_vs_dn.len1200_vs_0140" "posup_vs_dn.len1200_vs_0140" "padj"

# > names(rl2) ref=1200 dn
# [1] "posmd_vs_dn.len1200"   "posmd_vs_dn.len0140"   "posmd_vs_dn.len0350"   "posmd_vs_dn.len0700"   "posup_vs_dn.len1200"
# [6] "posup_vs_dn.len0140"   "posup_vs_dn.len0350"   "posup_vs_dn.len0700"   "len0140_vs_1200.posdn" "len0140_vs_1200.posmd"
# [11] "len0140_vs_1200.posup" "len0350_vs_1200.posdn" "len0350_vs_1200.posmd" "len0350_vs_1200.posup" "len0700_vs_1200.posdn"
# [16] "len0700_vs_1200.posmd" "len0700_vs_1200.posup" "betas_posdn.len1200"
# @examples
# ddsRd <- DESeq(dds0, test="LRT", reduced = ~ pos + len)
# rl=digReduceDDS(ddsRd)
digReduceDDS<-function(ddsRd, afactor1=NULL, afactor2=NULL, aref1=NULL, aref2=NULL) {
  #library(DEXSeq, verbose = FALSE)

  #[1] "Intercept"        "pos_md_vs_dn"<-factor1    "pos_up_vs_dn"     "len_0350_vs_0140"<-factor2 "len_0700_vs_0140" "len_1200_vs_0140"
  #[7] "posmd.len0350"    "posup.len0350"    "posmd.len0700"    "posup.len0700"    "posmd.len1200"    "posup.len1200"

  # Get diff of all levels of factor1~reference level under a level from factor2
  # Column names: Factor1levelX_vs_ref1.Factor2level2
  # For example, posmd_vs_dn.len0140 is posChange@len140, i.e., DE value of md-dn under len140
  # @return A resultList, each element like len0700_vs_0140.posmd.
  # @example
  # ##posChange@len
  # resList1=.getResOfF1ChangeAtF2(ddsRd,factor1='pos',factor2='len',levels1=levels(colData(ddsRd)[,'pos']),levels2)
  .getResOfF1ChangeAtF2<-function(ddsRd,factor1,factor2,levels1,levels2) {
    resList=list()
    ref1=levels1[1]
    ref2=levels2[1]
    coefNames=DESeq2::resultsNames(ddsRd)
    for (l1 in levels1[-1]) { #pos
      for (l2 in levels2) { #len
        if (l2==ref2) { #140
          name=paste(factor1,l1,'vs',ref1,sep='_') #"pos_md_vs_dn"
          res=DESeq2::results(ddsRd, name=name, test="Wald")

        } else {

          iname=paste(factor1,l1,'.',factor2,l2,sep='')
          if (!(iname %in% coefNames)) {
            iname=paste(factor2,l2,'.',factor1,l1,sep='')
          }
          name=c(paste(factor1,l1,'vs',ref1,sep='_'),iname) #c("pos_md_vs_dn","posup.len0350")
          res=DESeq2::results(ddsRd, contrast=list(name), test="Wald")

        }
        cname=paste0(factor1,l1,'_vs_',ref1,'.',factor2,l2)
        resList[[cname]]=res
        cat(cname,'\n')
      }
    }
    return(resList)
  }



  infos=getAndValidateFactorNamesFromDdsRd(ddsRd,afactor1,afactor2,aref1,aref2)
  factor1=infos$factor1
  factor2=infos$factor2
  levels1=infos$levels1
  levels2=infos$levels2

  ref1=levels1[1]
  ref2=levels2[1]

  #posChange@len
  cat(factor1,'change @',factor2,'\n')
  resList1=.getResOfF1ChangeAtF2(ddsRd,factor1,factor2,levels1,levels2)

  #lenChange@pos
  cat(factor2,'change @',factor1,'\n')
  resList2=.getResOfF1ChangeAtF2(ddsRd,factor1=factor2,factor2=factor1,levels1=levels2,levels2=levels1)

  #change change matrix
  cat('change change betas @',factor1, factor2, '\n')
  betas=digReduceDDSIntBeta(ddsRd, ostats=c('pvalue','padj'), factor1,factor2,ref1,ref2)
  cat(paste0(factor1,ref1,'.',factor2,ref2),'\n')

  resList=c(resList1, resList2)
  resList[[paste0('betas_',factor1,ref1,'.',factor2,ref2)]]=betas
  return(resList)
}

# Get change change matrix
# @param ddsRd The result from DESeqLRT.
# @param afactor1,afactor2,aref1,aref2 Only for validation
# @return list[betas(log2+pajd), pvalue, padj, each ostats is a matrix]
# betas: first few columns are refcol of factor1, ref of factor 2, changeFC of factor1Xfactor2, overall padj
# if factor1=len, factor2=pos, then cols is FACOTR1len_0350_vs_0140.posdn ..., FACTOR2pos_md_vs_dn.len0140, ..., F1XF2 len_0350_vs_0140.pos_md_vs_dn
# @examples
# res=digReduceDDSIntBeta(ddsRd)
# > names(res)
# [1] "betas"  "pvalue" "padj"
# > colnames(res$betas )
# [1] "len_0350_vs_0140.posdn"      "len_0700_vs_0140.posdn"      "len_1200_vs_0140.posdn"
# [4] "pos_md_vs_dn.len0140"        "pos_up_vs_dn.len0140"        "len0350_vs_0140.posmd_vs_dn"
# [7] "len0700_vs_0140.posmd_vs_dn" "len1200_vs_0140.posmd_vs_dn" "len0350_vs_0140.posup_vs_dn"
# [10] "len0700_vs_0140.posup_vs_dn" "len1200_vs_0140.posup_vs_dn" "padj"
digReduceDDSIntBeta<-function(ddsRd, ostats=c('pvalue','padj'), afactor1=NULL, afactor2=NULL, aref1=NULL, aref2=NULL) {

  #library(DEXSeq, verbose = FALSE)

  infos=getAndValidateFactorNamesFromDdsRd(ddsRd,afactor1,afactor2,aref1,aref2)
  factor1=infos$factor1
  factor2=infos$factor2
  levels1=infos$levels1
  levels2=infos$levels2

  betasRd <- coef(ddsRd)
  betasRd=betasRd[,-1] #-intercept
  cnames=colnames(betasRd)
  ref1=levels1[1]
  ref2=levels2[1]
  i1=grep(paste0('^',factor1,'.*','_vs_',ref1,'$'),cnames) #pos_md/up_vs_dn
  i2=grep(paste0('^',factor2,'.*','_vs_',ref2,'$'),cnames) #len_0350/.../_vs_0140
  i3=(1:length(cnames))[-c(i1,i2)] #posmd.len0350 or len0350.posmd

  cols1=paste0(cnames[i1],'.',factor2,ref2) #pos_md_vs_dn.len0140
  cols2=paste0(cnames[i2],'.',factor1,ref1) #len_0350_vs_0140.posdn
  cols3=strsplit(cnames[i3],split='\\.') #posmd.len0350
  l1=unlist(lapply(cols3,'[',1))
  l2=unlist(lapply(cols3,'[',2))
  if (!(l1[1] %in% paste0(factor1,levels1)) & !(l1[1] %in%  paste0(factor2,levels2))) {
    stop("l1 not in levels1 or levels2!")
  }
  if (l1[1] %in% paste0(factor1,levels1)) {
    l1=paste0(l1,'_vs_',ref1)
    l2=paste0(l2,'_vs_',ref2)
    cols3=paste0(l1,'.',l2)
  } else {
    l1=paste0(l1,'_vs_',ref2)
    l2=paste0(l2,'_vs_',ref1)
    cols3=paste0(l2,'.',l1)
  }

  #get ostats
  resList=list()
  for (name in cnames[c(i1,i2,i3)]) {
    res=DESeq2::results(ddsRd, name=name, test="Wald")
    resList[[name]]=res
  }
  resList=.getStatsFromResList(resList,cols=ostats)
  cnames[i1]=cols1
  cnames[i2]=cols2
  cnames[i3]=cols3
  colnames(betasRd)=cnames
  colnames(resList$pvalue)=cnames
  colnames(resList$padj)=cnames
  betasRd=betasRd[,c(i1,i2,i3)]
  betasRd=.asDf(betasRd)
  betasRd$padj=DESeq2::results(ddsRd)$padj
  resList=c(list(betas=betasRd),resList)
  return(resList)
}


# Given a DESeqdata, contrast all levels of factor 1 and 2, and get LRT results
# @param dds0 DESeqDataSet
# @param factor1, factor2: pos, len
# @param relevel1: If TRUE then relevel all, FALSE relevel the first.
# @param relevel2: If TRUE then relevel all, FALSE relevel the first.
#        The order of relevel of the two factors are irrelevent.
# @param ref1, ref2: If provided, then relevel respective factor to ref1 and ref2, then will not perform the for loop.
#                    Then relevel1 or relevel2 is not considered, this is the same as LRT and then digDESeq.
# @return A list, each element is relevel result from DEseq and then digReduce.
# @examples
# Only use 1200 and dn as reference.
# resAll=DESeqLRTAll(dds0=dds0, factor1='len', factor2='pos', ref1='1200',ref2='dn',relevel1=TRUE, relevel2=TRUE)
# Use all levels as reference one by one.
# resAll=DESeqLRTAll(dds0=dds0, factor1='len', factor2='pos', relevel1=TRUE, relevel2=TRUE)
# > names(resAll) ref_len0140.posdn is the reuslt under the reference 140/dn, perform DESeqLRT and digReduceDDS.
# [1] "ref_len0140.posdn" "ref_len0140.posmd" "ref_len0140.posup" "ref_len0350.posdn" "ref_len0350.posmd"
# [6] "ref_len0350.posup" "ref_len0700.posdn" "ref_len0700.posmd" "ref_len0700.posup" "ref_len1200.posdn"
# [11] "ref_len1200.posmd" "ref_len1200.posup"

# > names(resAll[[1]])
# [1] "len0350_vs_0140.posdn" "len0350_vs_0140.posmd" "len0350_vs_0140.posup" "len0700_vs_0140.posdn"
# [5] "len0700_vs_0140.posmd" "len0700_vs_0140.posup" "len1200_vs_0140.posdn" "len1200_vs_0140.posmd"
# [9] "len1200_vs_0140.posup" "posmd_vs_dn.len0140"   "posmd_vs_dn.len0350"   "posmd_vs_dn.len0700"
# [13] "posmd_vs_dn.len1200"   "posup_vs_dn.len0140"   "posup_vs_dn.len0350"   "posup_vs_dn.len0700"
# [17] "posup_vs_dn.len1200"   "betas_len0140.posdn"
DESeqLRTAll<-function(dds, factor1, factor2, ref1=NULL, ref2=NULL, relevel1=TRUE, relevel2=FALSE) {
  #library(DESeq2, verbose = FALSE)

  resAll=list()
  levels1=levels(colData(dds)[,factor1])
  levels2=levels(colData(dds)[,factor2])

  ddsi=dds
  design(ddsi)= as.formula(paste0('~',factor1,'*',factor2))  # ~pos*len
  if (relevel1) {
    levs1=levels1
  } else {
    levs1=levels1[1]
  }
  if (relevel2) {
    levs2=levels2
  } else {
    levs2=levels2[1]
  }

  if (!is.null(ref1)) {
    levs1=ref1
  }

  if (!is.null(ref2)) {
    levs2=ref2
  }

  for (l1 in levs1) { #Use "0140" "0350" "0700" as ref, respectively.
    for (l2 in levs2) {
      #the LRT model after changing the ref
      cat('relevel',factor1,'to',l1,',',factor2,'to',l2,'and call DESeq LRT\n')
      colData(ddsi)[,factor1]=relevel(colData(dds)[,factor1], ref=l1)
      colData(ddsi)[,factor2]=relevel(colData(dds)[,factor2], ref=l2)
      cname=paste0('ref_',factor1,l1,'.',factor2,l2)
      if (cname %in% names(resAll)) {
        cat('already done!\n')
        next
      }
      ddsRdi <- DESeq2::DESeq(ddsi, test="LRT", reduced = as.formula(paste0('~',factor1,'+',factor2)) )
      invisible(gc())
      resAll[[cname]]=digReduceDDS(ddsRdi)
    }
  }

  return(resAll)
}


# Dig DESeqLRTAll, to get matrix (ref) which have a reference level, matrix (pairwise) which includes all pairwise, and matirx (series) which has ordered levels.
# All results are under the change of factor1 the result of a level in factor2, betas is not included.
# @param resAll Results from DESeqLRTAll
# @param  factor1, factor2, ref1, ref2 Used to get matrix (ref) and to set the order of levels.
# @param pairwise, series Whether output the corresponding matrices.
# @return DEAllResults, with log2FC, pvalue, padj matrices, list[[ref=list(log2FC, padj), pairwise=..., series=...]]
# @examples
# ## get matrices (ref, pairwise and sereies)
# resLRT_heats=digDESeqLRTAll(resLRT_All, factor1='len', factor2='pos', ref1='0700', ref2='md', pairwise=TRUE, series=TRUE)
# class(resLRT_heats) #DEAllResults
# # ref/series/pairwise are a @slot with three dataframes recording log2FC, pvalue, padj of each condpair.
# names(resLRT_heats@series)
# # [1] "log2FC_len.pos" "pvalue_len.pos" "padj_len.pos"
# head(resLRT_heats@series$log2FC_len.pos)
# # len0140_vs_0700.posmd len0350_vs_0140.posmd len1200_vs_0350.posmd len0140_vs_0700.posdn len0350_vs_0140.posdn
# # PH02Gene00003            -1.5801450            -0.6438561             1.1395508            -1.6520768             0.5000734
# # PH02Gene00004            -5.6264305             0.2630228            -2.4809073            -2.9386030            -0.5849607
# pheatmap(resLRT_heats@series$log2FC_len.pos)
digDESeqLRTAll<-function(resAll, factor1, factor2, ref1, ref2, pairwise=TRUE, series=TRUE, ostats=c('log2FC'='log2FoldChange','pvalue'='pvalue','padj'='padj')) {
  #library(DESeq2, verbose = FALSE)

  # Get name of resAll under a given ref, e.g., ref_len1200.posup
  # @param resAll A list from DESeqLRTAll
  .getRefIdxFromResAll <- function (resAll,factor1, factor2, ref1, ref2) {
    allnames=names(resAll)
    aname=paste0('ref_',factor1,ref1,'.',factor2, ref2)
    if (aname %in% allnames) {
      return(aname)
    }
    aname=paste0('ref_',factor2,ref2,'.',factor1, ref1)
    if (aname %in% allnames) {
      return(aname)
    }
    stop(cat(aname, 'not in resAll\n'))
  }

  # @param digRes A list, from digREduceDDS or an element of resAll
  # @return list
  # $`len`
  # [1] "1200" "0140" "0350" "0700"
  # $pos
  # [1] "up" "dn" "md"
  .getFactorNamesFromDigRes <- function(digRes) {
    cnames=names(digRes)
    betaname=cnames[grep('^betas_',cnames)] #betas_len1200.posup
    f1=unlist(strsplit(betaname,'\\.'))
    f1[1]=gsub('betas_','',f1[1]) #"len1200" "posup"

    vsnames=gsub('\\..*','',cnames[-which(cnames==betaname)])
    vs=strsplit(vsnames,'_vs_')
    vs=unique(unlist(lapply(vs,'[',2))) #"1200" "up"
    ref1=vs[1]; ref2=vs[2]
    if (grepl(ref1,f1[1])) {
      factor1=gsub(ref1,'',f1[1])
      factor2=gsub(ref2,'',f1[2])
    } else {
      factor1=gsub(ref2,'',f1[2])
      factor2=gsub(ref1,'',f1[1])
    }
    vsnames=gsub(paste0('_vs_',ref1),'',vsnames)
    vsnames=gsub(paste0('_vs_',ref2),'',vsnames)
    levels1=c(ref1,unique(gsub(factor1,'',vsnames[grep(factor1,vsnames)])))
    levels2=c(ref2,unique(gsub(factor2,'',vsnames[grep(factor2,vsnames)])))
    ls=list(levels1, levels2)
    names(ls)=c(factor1,factor2)
    return(ls)
  }

  #MtxF1: under different factor2 levels, the change of factor1 (ref1 as control). E.g., under up/md/dn, the 140-1200, 350-1200...
  # [1] "len0140_vs_1200.posup" "len0350_vs_1200.posup" "len0700_vs_1200.posup" "len0140_vs_1200.posdn"
  # [5] "len0350_vs_1200.posdn" "len0700_vs_1200.posdn" "len0140_vs_1200.posmd" "len0350_vs_1200.posmd"
  # [9] "len0700_vs_1200.posmd"
  getMtxF1ChangeAtF2FromResAll<-function(resAll,factor1, factor2, ref1, ref2, ostats) {
    aname=.getRefIdxFromResAll(resAll,factor1, factor2, ref1, ref2) #ref_len1200.posup
    res=resAll[[aname]]
    annos=.getFactorNamesFromDigRes(res)
    levels1=annos[[factor1]] #"1200" "0140" "0350" "0700"
    levels2=annos[[factor2]] #"up" "dn" "md"
    if (levels1[1]!=ref1 | levels2[1]!=ref2) {
      stop('reflevel from digRes not equal user provided ref\n')
    }
    # > names(res)
    # [1] "len0140_vs_1200.posup" "len0140_vs_1200.posdn" "len0140_vs_1200.posmd" "len0350_vs_1200.posup"
    # [5] "len0350_vs_1200.posdn" "len0350_vs_1200.posmd" "len0700_vs_1200.posup" "len0700_vs_1200.posdn"
    # [9] "len0700_vs_1200.posmd" "posdn_vs_up.len1200"   "posdn_vs_up.len0140"   "posdn_vs_up.len0350"
    # [13] "posdn_vs_up.len0700"   "posmd_vs_up.len1200"   "posmd_vs_up.len0140"   "posmd_vs_up.len0350"
    # [17] "posmd_vs_up.len0700"   "betas_len1200.posup"
    at=paste0(factor2,levels2) #len0140_vs_1200.[posup]
    pre=paste0(factor1,levels1[-1],'_vs_',ref1) #[len0140_vs_1200].posup
    cnames=paste(pre,rep(at,each=length(pre)),sep='.')

    ls=.getStatsFromResList(res[cnames],cols=ostats)
    names(ls)=paste0(names(ostats),'_',factor1,'.',factor2)
    return(ls)
  }


  # Mtx2: pairwise of factor1 under each of factor2
  # i.e., from results using different levels of factor1 as ref, get len0140_vs_1200.posdn len1200_vs_len0140.posdn...
  # And then remove duplicated names
  # [1] "len0140_vs_1200.posup" "len0350_vs_1200.posup" "len0700_vs_1200.posup" "len0350_vs_0140.posup" "len0700_vs_0140.posup"
  # [6] "len0700_vs_0350.posup" "len0140_vs_1200.posdn" "len0350_vs_1200.posdn" "len0700_vs_1200.posdn" "len0350_vs_0140.posdn"
  # [11] "len0700_vs_0140.posdn" "len0700_vs_0350.posdn" "len0140_vs_1200.posmd" "len0350_vs_1200.posmd" "len0700_vs_1200.posmd"
  # [16] "len0350_vs_0140.posmd" "len0700_vs_0140.posmd" "len0700_vs_0350.posmd"
  getMtxF1PairChangeAtF2FromResAll<-function(resAll,factor1, factor2, ref1, ref2, ostats) {
    aname=.getRefIdxFromResAll(resAll,factor1, factor2, ref1, ref2) #ref_len1200.posup
    annos=.getFactorNamesFromDigRes(resAll[[aname]])
    f1PairAtf2=list()
    for (l1 in annos[[factor1]]) {
      fi=getMtxF1ChangeAtF2FromResAll(resAll,factor1, factor2, ref1=l1, ref2=ref2, ostats=ostats)
      if (l1==annos[[factor1]][1]) {
        f1PairAtf2=fi
      } else {
        for (i in 1:length(ostats)) {
          f1PairAtf2[[i]]=cbind(f1PairAtf2[[i]],fi[[i]])
        }
      }
    }
    cnames=c()
    for (i in 1:(length(annos[[factor1]])-1)) {
      for (j in (i+1):(length(annos[[factor1]]))) {
        cnames=c(cnames,paste0(factor1,annos[[factor1]][j],'_vs_',annos[[factor1]][i]))  #len0350_vs_0140.posup
      }
    }
    cnames=paste(cnames,rep(paste0(factor2,annos[[factor2]]),each=length(cnames)),sep='.')
    for (i in 1:length(ostats)) {
      f1PairAtf2[[i]]=f1PairAtf2[[i]][,cnames]
    }
    return(f1PairAtf2)
  }

  ls=list()

  f1Atf2=getMtxF1ChangeAtF2FromResAll(resAll,factor1, factor2, ref1, ref2, ostats)
  ls[['ref']]=f1Atf2

  if (pairwise | series) {
    f1PairAtf2=getMtxF1PairChangeAtF2FromResAll(resAll,factor1, factor2, ref1, ref2, ostats)
    #colnames(f1PairAtf2[[1]])

    if (pairwise) {
      ls[['pairwise']]=f1PairAtf2
    }

    # Mtx3: series matrix
    # the next level of factor1~ previous level of  factor 1 under different levels of factor2.
    # Here ref1=1200, so the order is 1200 140 350 700,  and the first is 140~1200, 350~140...
    # [1] "len0140_vs_1200.posup" "len0350_vs_0140.posup" "len0700_vs_0350.posup" "len0140_vs_1200.posdn" "len0350_vs_0140.posdn"
    # [6] "len0700_vs_0350.posdn" "len0140_vs_1200.posmd" "len0350_vs_0140.posmd" "len0700_vs_0350.posmd"
    if (series) {
      aname=.getRefIdxFromResAll(resAll,factor1, factor2, ref1, ref2) #ref_len1200.posup
      annos=.getFactorNamesFromDigRes(resAll[[aname]])
      cnames=paste0(annos[[factor1]][2:length(annos[[factor1]])],'_vs_',annos[[factor1]][1:(length(annos[[factor1]])-1)])
      cnames=paste0(factor1,cnames)
      cnames=paste(cnames,rep(paste0(factor2,annos[[factor2]]),each=length(cnames)),sep='.')
      f1SeriesAtf2=f1PairAtf2
      for (i in 1:length(ostats)) {
        f1SeriesAtf2[[i]]=f1SeriesAtf2[[i]][,cnames]
      }
      ls[['series']]=f1SeriesAtf2
    }

  }#~pair

  DEAllRes <- new("DEAllResults",ref=ls$ref, pairwise=ls[['pairwise']], series=ls[['series']])
  return(DEAllRes)
}


# ---------------- *** DESeq2 one/two factors *** ------------------

# For condpair names (normally is resultsNames from dds), get all levels of a given factor.
# And set given ref1 as the first level, other levels are sorted by character order.
# @param rn Such as len_0350_vs_0140
# @param factor1 Group name like len
# @parma ref1 If provided, then sort levels vs ref1
# @return levels
# @examples
#rn=resultsNames(ddsRes)[-1]
#levels1=.getLevelsFromResultNames(rn=c("len_0350_vs_0140","len_0700_vs_0140","len_1200_vs_0140"), factor1='len', ref1=NULL)
#levels1=.getLevelsFromResultNames(rn=c("len_0350_vs_0140","len_0700_vs_0140","len_1200_vs_0140"), factor1='len', ref1='1200')
.getLevelsFromResultNames<-function(rn, factor1, ref1=NULL) {
  f1rn=unlist(strsplit(rn[grep(factor1,rn)],'_vs_'))
  levels1=unique(gsub(paste0(factor1,'_'),'',f1rn))
  if (!is.null(ref1)) {
    if (!(ref1 %in% levels1)) {
      stop(cat('ref1 not in results(ddsRes)\n'))
    }
  }
  levels1=sort(levels1)
  if (!is.null(ref1)) {
    levels1=levels1[-which(levels1==ref1)]
    levels1=c(ref1,levels1)
  }
  return(levels1)
}

# normalize DESeq2 dds
# @param dds: DESeqdataset or DEXSeqDataSet
# @param norm: FALSE. If TRUE, then call estimateSizeFactors(no matter whether DE has been called before)
#              If there is no sizefactor in dds and norm=FALSE, then will not do normalization (i.e., set factor=1)
# @param force1 To force normalize factor=1.
normalizeDDS <- function (dds, norm, force1=FALSE) {
  if (inherits(dds, "DEXSeqDataSet")) {
    nsmp=ncol((SummarizedExperiment::assays(dds)$counts))/2
  } else {
    nsmp=ncol((SummarizedExperiment::assays(dds)$counts))
  }

  if (force1) {
    sizeFactors(dds)=rep(1,nsmp)
    return(dds)
  }

  if (norm) {
    #cat("Call estimateSizeFactors\n")
    dds=DESeq2::estimateSizeFactors(dds)
  }  else {
    if (is.null(sizeFactors(dds)) & is.null(DESeq2::normalizationFactors(dds))) { #
      sizeFactors(dds)=rep(1,nsmp)
      #cat('Using normfactor=1\n')
      # gene-specific size factors
      # normFactors <- matrix(rep(1,nrow(dds)*ncol(dds)),
      #                       ncol=ncol(dds),nrow=nrow(dds),
      #                       dimnames=list(1:nrow(dds),1:ncol(dds)))
      # normalizationFactors(dds) <- normFactors
    }
  }
  return(dds)
}


# Perform DEseq without interaction
# @param dds DESeqDataSet
# @param factor1 len
# @param factor2 If NULL, then only have ~factor1
# @param norm: FALSE. If TRUE, then call estimateSizeFactors(no matter whether DE has been called before)
#              If there is no sizefactor in dds and norm=FALSE, then will not do normalization (i.e., set factor=1)
# @return A DESeq result. The only diffrence from DESeqLRTAll is this is without interaction, then one DESeq can get all pairwise results.
# @examples
# Ex. two factors
# ddsRd=DESeqNoInteractionAll(dds0,factor1='len',factor2='pos')
# > resultsNames(ddsRd)
# [1] "Intercept" "len_0350_vs_0140" "len_0700_vs_0140" "len_1200_vs_0140" "pos_md_vs_dn" "pos_up_vs_dn"
# Ex. one factor
# ddsRd=DESeqNoInteractionAll(dds0,factor1='len')
# > resultsNames(ddsRd)
# [1] "Intercept"        "len_0350_vs_0140" "len_0700_vs_0140" "len_1200_vs_0140"
DESeqNoInteractionAll<-function(dds, factor1, factor2=NULL, norm=FALSE) {
  #library(DESeq2, verbose = FALSE)
  dds=normalizeDDS(dds, norm=norm, force1=FALSE)

  mcols(dds)=NULL #TODO otherwise will error
  if (!is.null(factor2)) {
    levels1=levels(colData(dds)[,factor1])
    levels2=levels(colData(dds)[,factor2])
    ddsi=dds
    design(ddsi)= as.formula(paste0('~',factor1,'+',factor2))  # ~pos+len
  } else {
    levels1=levels(colData(dds)[,factor1])
    ddsi=dds
    design(ddsi)= as.formula(paste0('~',factor1))  # ~pos
  }

  ddsRdi <- DESeq2::DESeq(ddsi)
  return(ddsRdi)
}


# Get all results from DESeqNoInteractionAll, see digDESeqAll()

# Ex. ~ len+pos (there is no meaning for factor2 in the result,
# just that in [[ref]], the name will be added a sufix of .factor2, e.g., log2FC_len.pos.
# ddsRd=DESeqNoInteractionAll(dds0,factor1='len',factor2='pos')

# mtxAll=digDESeqNoInteractionAll(ddsRd, factor1='len', factor2='pos', ref1='1200', ref2='up', pairwise=TRUE, series=TRUE, levels1=NULL, levels2=NULL)
# > colnames(mtxAll@ref$log2FC_len.pos)
# [1] "len0140_vs_1200" "len0350_vs_1200" "len0700_vs_1200"
# > colnames(mtxAll@pairwise$log2FC_len.pos)
# [1] "len0140_vs_1200" "len0350_vs_1200" "len0700_vs_1200" "len0350_vs_0140" "len0700_vs_0140" "len0700_vs_0350"
# > colnames(mtxAll@series$log2FC_len.pos)
# [1] "len0140_vs_1200" "len0350_vs_0140" "len0700_vs_0350"
#
# mtxAll=digDESeqNoInteractionAll(ddsRd, factor1='len', factor2='pos', ref1='1200', ref2='up', pairwise=TRUE, series=TRUE, levels1=c('1200','0700','0350','0140'), levels2=NULL)
# > colnames(mtxAll@ref$log2FC_len)
# [1] "len0700_vs_1200" "len0350_vs_1200" "len0140_vs_1200"
# > colnames(mtxAll@series$log2FC_len)
# [1] "len0700_vs_1200" "len0350_vs_0700" "len0140_vs_0350"

# mtxAll=digDESeqNoInteractionAll(ddsRd, factor1='pos', factor2='len', ref1='dn', ref2='1200', pairwise=TRUE, series=TRUE)
# > colnames(mtxAll@ref$log2FC_pos)
# [1] "posmd_vs_dn" "posup_vs_dn"
# > colnames(mtxAll@pairwise$log2FC_pos)
# [1] "posmd_vs_dn" "posup_vs_dn" "posup_vs_md"
# > colnames(mtxAll@series$log2FC_pos)
# [1] "posmd_vs_dn" "posup_vs_md"

# ~len
# ddsRd=DESeqNoInteractionAll(dds0,factor1='len')
# mtxAll=digDESeqNoInteractionAll(ddsRd, factor1='len', factor2=NULL, ref1='1200', ref2=NULL, pairwise=TRUE, series=TRUE)
# > colnames(mtxAll@ref$log2FC_len)
# [1] "len0140_vs_1200" "len0350_vs_1200" "len0700_vs_1200"
digDESeqNoInteractionAll<-function(ddsRes, factor1, factor2=NULL, ref1=NULL, ref2=NULL, pairwise=TRUE, series=TRUE,
                                   levels1=NULL, levels2=NULL, ostats=c('log2FC'='log2FoldChange','pvalue'='pvalue','padj'='padj'),
                                   verbose=FALSE) {

  #ddsRes
  #[1] "Intercept"        "pos_dn_vs_md"     "pos_up_vs_md"     "len_0350_vs_0140" "len_0700_vs_0140" "len_1200_vs_0140"

  if (verbose) cat('get ref/pairwise/series change matrix from non-interaction DDS results\n')
  rn=DESeq2::resultsNames(ddsRes)[-1]

  .getLevels<-function(levels1=NULL, rn, factor1, ref1=NULL) {
    if (is.null(levels1)) {
      levels1=.getLevelsFromResultNames(rn,factor1,ref1)
    }

    if (!is.null(ref1)) {
      if (!(ref1 %in% levels1)) {
        stop(cat('ref1',ref1,'not in levels1',levels1[1],'\n'))
      }

      if (ref1!=levels1[1]) {
        if (levels1[1]!=ref1) {
          levels1=levels1[-which(levels1==ref1)]
          levels1=c(ref1,levels1)
        }
      }
    }
    return(levels1)
  }

  levels1=.getLevels(levels1, rn, factor1, ref1)
  if (!is.null(factor2))  {
    levels2=.getLevels(levels2, rn, factor2, ref2)
    cat('levels1=',levels1,'\nlevels2=',levels2,'\n')
  } else {
    cat('levels=',levels1,'\n')
  }

  ls=list()

  #ls[[ref]]=list(log2FC_len, padj_len)
  #Get change of factor1(xx~ref1) log2FC_len=matrix[len0140_vs_0350 len0700_vs_0350 len1200_vs_0350]
  ls[['ref']]=list()
  if (!is.null(ref1)) {
    resList=list()
    for (lv in levels1[-1]) {
      res=DESeq2::results(ddsRes, contrast=c(factor1, lv, ref1))
      resList[[paste0(factor1,lv,'_vs_',ref1)]]=res
    }

    ls[['ref']]=.getStatsFromResList(resList,cols=ostats)
    if (!is.null(factor2)) {
      names(ls[['ref']])=paste0(names(ostats),'_',factor1,'.',factor2)
    } else {
      names(ls[['ref']])=paste0(names(ostats),'_',factor1)
    }
  }

  ls[['pairwise']]=list()
  if (pairwise) {
    resList=list()
    for (i in 1:(length(levels1)-1)) {
      for (j in (i+1):length(levels1)) {
        res=DESeq2::results(ddsRes, contrast=c(factor1, levels1[j], levels1[i]))
        resList[[paste0(factor1,levels1[j],'_vs_',levels1[i])]]=res
      }
    }
    ls[['pairwise']]=.getStatsFromResList(resList,cols=ostats)
    if (!is.null(factor2)) {
      names(ls[['pairwise']])=paste0(names(ostats),'_',factor1,'.',factor2)
    } else {
      names(ls[['pairwise']])=paste0(names(ostats),'_',factor1)
    }
  }

  ls[['series']]=list()
  if (series) {
    resList=list()
    for (i in 1:(length(levels1)-1)) {
      res=DESeq2::results(ddsRes, contrast=c(factor1, levels1[i+1], levels1[i]))
      resList[[paste0(factor1,levels1[i+1],'_vs_',levels1[i])]]=res
    }

    ls[['series']]=.getStatsFromResList(resList,cols=ostats)
    if (!is.null(factor2)) {
      names(ls[['series']])=paste0(names(ostats),'_',factor1,'.',factor2)
    } else {
      names(ls[['series']])=paste0(names(ostats),'_',factor1)
    }
  }

  DEAllRes <- DEAllResults(ref=ls[['ref']], pairwise=ls[['pairwise']], series=ls[['series']])
  return(DEAllRes)
}




# ------------ *** DEXseq *** -----------------

# Convert PACdataset to DEXdataset
# @param PACds A PACdataset object
# @param group a group name (column in colData). If group=NULL, then use 1 or 2 factors in colData depends on the ncol in colData.
# If design has two conditions, then the formular will have two *, otherwise uses only the formular with the first condition.
# @example
# dxd=PACdataset2Dxd(PACds)
PACdataset2Dxd<-function(PACds, group=NULL) {
  #library(DEXSeq, verbose = FALSE)
  sampleData <- PACds@colData

  if (is.null(group)) {
    #TODO: two factors
    if (ncol(sampleData)==2) {
      factor1=colnames(PACds@colData)[1]
      factor2=colnames(PACds@colData)[2]
      design <- as.formula( sprintf("~ %s + %s + exon + %s:exon + %s:exon", factor1, factor2, factor2, factor1))
    } else {
      factor1=colnames(PACds@colData)[1]
      design <- as.formula( sprintf("~ %s + exon + %s:exon", factor1, factor1))
    }

    #design <- as.formula( sprintf("~ %s + exon + %s:exon", group, group))

  } else {
    factor1=group
    design <- as.formula( sprintf("~ %s + exon + %s:exon", factor1, factor1))
  }

  groupID <- PACds@anno$gene
  featureID=as.character(PACds@anno$coord);
  #print(head(cbind(groupID,featureID)));
  dxd= DEXSeq::DEXSeqDataSet( PACds@counts, sampleData, design, featureID, groupID )

  #if (disp) dxd = estimateDispersions(dxd) #No normfactor, cannot calc dispersion
  return(dxd)
}

# Perform DEXseq for all levels in one factor
# DEXSeq will not get a padj for each pairwise, but only a overall padj for all pairs.
# BUT, we can use estimatefoldchange to get all pairwise FC, and then call getDEXlist() to ROUGHLY get DElist from a pair conds.
# @param dxd DEXseqDataset
# @param factor1 One column in colData(dxd), e.g., len
# @param reflevel Base level for log2FC, like 1200, then log2fold_0140_1200
# @param alllevel TRUE then estimateExonFoldChanges for all levels
# @param norm: FALSE. If TRUE, then call estimateSizeFactors(no matter whether DE has been called before)
#              If there is no sizefactor in dds and norm=FALSE, then will not do normalization (i.e., set factor=1)
# @param returnDxr TRUE then return DEXSeqResults object, else return new dxd.
# @return DEXSeqResults object or DEXseqDataset
# @examples
# dxr1=DEXSeqOneFactor(dxd,factor1='len')
# dxr2=DEXSeqOneFactor(dxd,factor1='len', reflevel='1200')
# dxr3=DEXSeqOneFactor(dxd,factor1='len', reflevel='1200', norm=TRUE)

# dxd=DEXSeqOneFactor(dxd,factor1='len', reflevel='1200', norm=TRUE, returnDxr=FALSE)
# estimateFC for all levels
# dxd=DEXSeqOneFactor(dxd,factor1='len', reflevel='', alllevel=TRUE, norm=TRUE, returnDxr=FALSE)
DEXSeqOneFactor<-function(dxd, factor1, reflevel="", alllevel=TRUE, norm=FALSE, returnDxr=TRUE) {
  #library(DEXSeq, verbose = FALSE)
  dxd=normalizeDDS(dxd, norm=norm)

  ## one factor ~ len
  formulaFullModel = as.formula(sprintf("~ %s + exon + %s:exon", factor1,factor1))
  formulaReducedModel =as.formula(sprintf("~ %s + exon", factor1))  #~ len + exon

  if (is.null(DESeq2::dispersions(dxd))) {
    cat("Call estimateDispersions using all samples\n")
    dxd = DEXSeq::estimateDispersions(dxd)
  }

  dxd = DEXSeq::testForDEU( dxd,
                    reducedModel = formulaReducedModel,
                    fullModel = formulaFullModel )

  if (!alllevel) {
    cat('estimate FC for reflevel\n')
    dxd = DEXSeq::estimateExonFoldChanges( dxd, fitExpToVar=factor1, denominator = reflevel)
  } else {
    cat('estimate FC for all levels\n')
    lvs=levels(dxd@colData[,factor1])
    for (lv in lvs) {
      cat(lv,'estimating...\n')
      dxd = DEXSeq::estimateExonFoldChanges( dxd, fitExpToVar=factor1, denominator = lv)
    }
  }

  if (returnDxr) {
    return(DEXSeq::DEXSeqResults( dxd ))
  }
  return(dxd)
}

# Get DE list from DEXSeqResults
# @param dxr: DEXSeqResults
# @param cond1, cond2: must be in dxr; log2_cond2_cond1 (cond2/cond1)
# @param padj: padj<0.1
# @param logFC: if null then get minMaxfc of padj<XX as logFC
# @param full: FALSE; if True, then output list including all information (list[[1]]=DElist)
# @param final DE list = (padj<0.1 & (logFC=maxLogFC | logFC>providedLogFC))
#   NOTE: will get ALL columns start with log2fold.
#         Even if a dxd may be called more than one time then it would have many log2FC_ columns. But all these columns will be extracted.
# @return rownames(DE). The output is gene:coord, because dxd use gene as gene and coord as exonID.
#   NOTE: DEXseq only get one padj for all pairs, so there is no exact result for each pair cond.
#         Only if you call DEXseq for each pairwise cond one by one, but the speed is SLOW.
#         Even logFC cutoff is set, the final DEs are not always >=logFC. Because rows with FC=maxFC will be considered as DE.
#         The criteria for DE is >=logFC or =maxFC.
# @examples
# de=getDEXlist(dxr, cond1='0140', cond2='1200', padj=0.1, logFC=NULL)
# head(dxr[de,c('padj',colnames(dxr)[grep('^log2fold',colnames(dxr))])])
# de=getDEXlist(dxr, cond1='0140', cond2='1200', padj=0.1, logFC=0.5)

# dxd=DEXSeqOneFactor(dxd,factor1='len', reflevel="", norm=FALSE, returnDxr=FALSE)
# dxd = estimateExonFoldChanges( dxd, fitExpToVar="len", denominator = '0350')
# dxr=DEXSeqResults( dxd )
# de=getDEXlist(dxr, cond1='0140', cond2='1200', padj=0.1, logFC=NULL)

# pacDEsFull=getDEXlist(pacDERes, cond1='0140', cond2='1200', padj=0.1, logFC=NULL, full=TRUE)
# $DE
# "PH02Gene06999:11637629" "PH02Gene06961:10665660"
# $sample
# [1] "log2fold_1200_0140"
# $switch
# [1] FALSE
# $cutoff
# [1] "autoMin"
# $logFCThd
# [1] 0.09499061
getDEXlist<-function(dxr, cond1, cond2, padj=0.1, logFC=NULL, full=FALSE) {

  fccols=colnames(dxr)[grep('^log2fold',colnames(dxr))] #log2fold_0350_0140
  fccols=strsplit(fccols,'_')
  fccols=unique(c(unlist(lapply(fccols,'[',2)),unlist(lapply(fccols,'[',3))))
  if (!(cond1 %in% fccols)) stop(cat("cond1",cond1,"not in dxr",fccols,'\n'))
  if (!(cond2 %in% fccols)) stop(cat("cond2",cond2,"not in dxr",fccols,'\n'))

  #log2fold_0350_0140
  fccol=sprintf("log2fold_%s_%s",cond2,cond1)
  switch=FALSE
  if (fccol %in% colnames(dxr)) {
    fc=dxr[,fccol]
  } else {
    fccol=sprintf("log2fold_%s_%s",cond1,cond2)
    if (fccol %in% colnames(dxr)) {
      fc=dxr[,fccol]*(-1)
    } else {
      stop("no log2 cond1 cond2 in dxr, use estimateExonFoldChanges to get log2FC")
      #dxd = estimateExonFoldChanges( dxd, fitExpToVar="len", denominator = reflevel)
    }
    switch=TRUE
  }

  dxr=.asDf(dxr)
  if (!is.null(padj)) {
    idx=which(!is.na(dxr$padj) & dxr$padj<padj)
    if (length(idx)==0) {return(NULL)}
    dxr=dxr[idx,]
  }

  maxfc=apply(abs(dxr[,colnames(dxr)[grep('^log2fold',colnames(dxr))]]),1,max)
  minMaxfc=min(maxfc)
  fc=abs(fc[idx])
  txt='set'
  if (is.null(logFC)) {
    logFC=minMaxfc
    txt='autoMin'
  }
  idx=which(fc==maxfc | fc>=logFC)
  #If logfFC is the max FC of the row, then it is definitely DEX.
  #Otherwise it should with logFC>cutoff.
  if (length(idx)==0) {return(NULL)}
  cat(cond1,'~',cond2,'padj=',padj,'DEXlogFC(',txt,')=',logFC, 'nDEXPAC=',length(idx),'switchLog=',switch,'\n')
  if (!full) {
    return(rownames(dxr[idx,]))
  }

  dxr=list(DE=rownames(dxr[idx,]), sample=fccol, switch=switch, cutoff=txt, logFCThd=logFC)
  return(dxr)
}


# ---------------- *** APA switching *** ------------------
#  FU
#  @param s1 The PAT# from sample 1 and sample 2.
#  @param score The 3UTR length of PACs.
# s1=c(17066,14464,788,126,37)
# s2=c(48,38,5,1,1)
# score=c(0,0.5,1.5,4.0,7.0)
# r=FU(s1,s2,score)
# pvalue:
#M=(sum(s1+s2)-1)*(r^2)
#prop.trend.test(s1, s1+s2,score=score)$p.value #using R
#1-pchisq(M,1) #using chisq-test
#> s1
#[1] 50 30 20
#> s2
#[1]  30  35 200
#> score
#[1] 100 150 300
#> r
#0.5179539

#>  log2(sum(s1)/sum(s2))
#[1] -1.405992
#FU(c(50,30,20),c(30,35,200),c(100,150,300))
#FU(c(100,50),c(50,100),c(100,200)) #0.33  s1 to s2 longer
#FU(c(50,100),c(100,50),c(100,200)) #-0.33 s1 to s2 shorter
#FU(c(5,100),c(100,5),c(100,200)) #-0.9 s1 to s2 shorter
#FU(c(5,100,200),c(100,5,1),c(100,150,200)) #-0.82 s1 to s2 shorter
#FU(c(5,200),c(100,1),c(100,200)) #-0.95 s1 to s2 shorter
FU<-function(s1,s2,score) {
  p1=s1/sum(s1+s2)
  p2=s2/sum(s1+s2)
  p=rbind(p1,p2)
  u=c(1,2)
  v=score
  ubar=sum(p1)*u[1]+sum(p2)*u[2]
  vbar=sum((p1+p2)*v)

  fz=0
  for (i in 1:2) {
    for (j in 1:length(v)) {
      fz=fz+(u[i]-ubar)*(v[j]-vbar)*p[i,j]
    }
  }

  fm1=0
  for (i in 1:2) {
    fm1=fm1+(u[i]-ubar)^2*sum(p[i,])
  }

  fm2=0
  for (j in 1:length(v)) {
    fm2=fm2+(v[j]-vbar)^2*sum(p[,j])
  }

  r=fz/sqrt(fm1*fm2)
  return(r)
}


#Fu switching
# @param PACds PACdataset
# @param group one column in PACds@colData, like len
# @param cond1, cond2: 0140, 1200; level in colData of group
# @param avgPACtag 5, after get cond1 and cond2, the average PAT tag >=5
# @param avgGeneTag 10, after get cond1 and cond2, the average gene tag >=10
# @return
# gene, nPAC, geneTag1, geneTag2, avgUTRlen1, avgUTRlen2, pvalue, padj,
# change <linear: cor, logRatio=log2(geneTag@smp1/geneTag@smp2)>, PAs1, PAs2
# cr>0 (change=1) s1 to s2 longer
UTRtrend_linear<-function(PACds, group, cond1, cond2, avgPACtag=5, avgGeneTag=10) {

  ds=subsetPACds(PACds, group=group, cond1=cond1, cond2=cond2, avgPACtag=avgPACtag, avgGeneTag=avgGeneTag, avg=TRUE)
  ds=get3UTRAPAds(ds)

  gstat=geneStatForPACds(ds, avgUTRlen=TRUE, geneTag=TRUE, PAinfo=TRUE, merge=TRUE)
  colnames(gstat)=c('gene','avgUTRlen1', 'avgUTRlen2','nPAC','geneTag1', 'geneTag2','PAs1', 'PAs2')


  genes=unique(ds@anno$gene)
  p=c()
  cr=c()
  ratio=c()

  dat=.asDf(ds@counts)
  dat$UTRlen=ds@anno$toStop
  smp1=cond1
  smp2=cond2
  for (g in genes) {
    ig=which(ds@anno$gene==g)
    s1=dat[ig,smp1]
    s2=dat[ig,smp2]
    score=dat[ig,'UTRlen']
    if (sum(s1)==0 | sum(s2)==0) {
      p=c(p,NA); cr=c(cr,NA); ratio=c(ratio,NA)
      next
    }

    p=c(p, suppressWarnings(prop.trend.test(s1, s1+s2, score=score))$p.value)
    cr=c(cr,FU(s1,s2,score))
    ratio=c(ratio,log2((sum(s1)+0.1)/(sum(s2)+0.1)))
  }
  pv=p
  pv[!is.na(pv)]=p.adjust(p[!is.na(p)])

  rt=.asDf(cbind(cor=cr,logRatio=ratio,pvalue=p, padj=pv))
  rownames(rt)=genes
  rt$change=cr
  rt$change[!is.na(cr) & cr>0]=1 #smp1 to smp2 longer
  rt$change[!is.na(cr) & cr<=0]=-1 #smp1 to smp2 shorter
  rt$gene=rownames(rt)
  if (nrow(rt)!=nrow(gstat)) {
    stop("nrow(rt)!=nrow(gstat)!")
  }
  rt=merge(rt,gstat,by.x='gene',by.y='gene')
  rt=rt[,c('gene', 'nPAC', 'geneTag1', 'geneTag2', 'avgUTRlen1', 'avgUTRlen2', 'pvalue', 'padj', 'change', 'cor', 'logRatio', 'PAs1', 'PAs2')]
  rownames(rt)=rt[,'gene']
  return(rt)
}

#
# Get all DEpairs from a given gene-PA dataframe
# @param genePAs Df with [(smp1,smp2) must be the first two columns regardless what names, (gene,isDEPA 1/0,coord,PA)]
# @param mindist coord1-coord2+1>=mindist
# @param nDEPAC=1/2, then i or/and j is isDEPA
# @param pseudo added for count
# @param fisherThd: filter DE list by pvalue<fisherThd
# @param logFCThd: filter |rcij| >=logFCThd, which is Log2(ratio of PA1 between sample A~B / ratio PA2)
#        NOTE: The switching is not neccesary crossing, but only with sufficient difference of the ratio between PA1 and PA2.
#              E.g., it may be that both PA1 and PA2 are used more in sampleA, but PA2 used much more than PA2.
# @param cross FALSE.
# @param selectOne If NULL, output all PApairs; =farest choose the pair with the longest distance; =logFC choose biggest; =fisherPV choose the smallest PV.
# @return Dataframe [gene PA1 PA2 dist nDEPA nSwitchPair logFC fisherPV] 1-row if selectOne
# @examples
# DE=ddply(d, .(gene), getAllSwitchingPairsForOneGene, mindist=50, nDEPAC=1, pseudo=0.1, fisherThd=1, logFCThd=0, selectOne='logFC')
# gene               PA1               PA2 dist nDEPA nSwitchPair       logFC   fisherPV change
# 1  PH02Gene00022 PH02Gene00022.PA1 PH02Gene00022.PA2  374     2           1 -0.41404606 0.51879565     -1
# 2  PH02Gene00068 PH02Gene00068.PA1 PH02Gene00068.PA2   66     2           1  0.55019543 0.75958929      1
getAllSwitchingPairsForOneGene<-function (genePAs, mindist=50, nDEPAC=1, pseudo=1,
                                          fisherThd=0.01, logFCThd=1, cross=FALSE,
                                          selectOne=NULL) { #c('farest','logFC','fisherPV')
  genePAs$gene <- as.character(genePAs$gene)
  genePAs$PA = as.character(genePAs$PA)
  genePAs[,1:2]=genePAs[,1:2]+pseudo
  ocols=c('gene','PA1','PA2','dist','nDEPA','nSwitchPair','logFC','fisherPV','change')
  opair=matrix(nrow=0,ncol=length(ocols))
  colnames(opair)=ocols
  opair=.asDf(opair)
  opair$gene = as.character(opair$gene)
  opair$PA1 = as.character(opair$PA1)
  opair$PA2 = as.character(opair$PA2)

  if (nrow(genePAs)<=1) {
    return(opair)
  }

  for (i in 1:(nrow(genePAs)-1)) {
    for (j in (i+1):(nrow(genePAs))) {
      dist=abs(genePAs$coord[j]-genePAs$coord[i])+1
      if (dist<mindist) {
        next
      }
      nDEPA=genePAs$isDEPA[i]+genePAs$isDEPA[j] #ijDE#
      if (length(nDEPA)==0) nDEPA=0
      if (nDEPAC>0 & nDEPA<nDEPAC) {
        next
      }
      #TAPAS: PAi~PAj: rcij=log2(ejB/ejA)-log2(eiB/eiA)  |rcij|>=1 as switching
      #iPA is proximal, jPA is distal PA; rcij<0 AtoB shorter = -1; >0 longer =1
      fcj=log2(genePAs[j,2]/genePAs[j,1])
      fci=log2(genePAs[i,2]/genePAs[i,1])
      logFC=fcj-fci
      if (cross) { #cross
        if (fcj*fci>0) {
          next
        }
      }
      fisherPV=fisher.test(round(genePAs[c(i,j),1:2]))$p.value
      if (is.na(fisherPV)) fisherPV=1
      if (abs(logFC)>=logFCThd & fisherPV<fisherThd) {
        if (logFC>0) change=1 else change=-1
		opair$gene = as.character(opair$gene)
        opair=rbind(opair,c(gene=genePAs$gene[1],PA1=genePAs$PA[i],PA2=genePAs$PA[j], dist=dist, nDEPA=nDEPA, nSwitchPair=1,logFC=logFC, fisherPV=fisherPV, change=change))
      }
    }
  }


  colnames(opair)=ocols
  opair=.asDf(opair)

  if (nrow(opair)==0) {
    return(opair)
  }


  opair$nSwitchPair=nrow(opair)
  opair[,c('dist','logFC','fisherPV')]=apply(opair[,c('dist','logFC','fisherPV')],2,as.numeric)
  if (nrow(opair)>1 & !is.null(selectOne)) {
    if (selectOne=='farest') {
      i=which(opair$dist==max(opair$dist))
    } else if (selectOne=='logFC') {
      i=which(abs(opair$logFC)==max(abs(opair$logFC)))
    } else if (selectOne=='fisherPV') {
      i=which(opair$fisherPV==min(opair$fisherPV))
    }
    opair=opair[i[1], , drop=F]
  }
  colnames(opair)=ocols
  return(opair)
}


# UTR lengthening by DE, only two PAs in a gene are considered and at least one PA is DE according to DESeq2 or DEXSeq.
# First get switching PA pairs for each gene, and then get one pair from all switching pairs by selectOne param.
# @param PACds: PACdataset
# @param group: one column in PACds@colData, "len"
# @param cond1, cond2: 0140, 1200; level in colData of group
# @param avgPACtag=5
# @param avgGeneTag=10
# @param DEpadjThd=0.1: to filter DEPAC with padj<0.1
# @param DEXlogFC: to filter DEXPAC with (padj<0.1 and logFC>=DEXlogFC)
# @param aMovDEPACRes movDEPACRes
# @param mindist=50, fisherThd=0.1, logFCThd=1: to filter all switching pair
# @param selectOne=#c('farest','logFC','fisherPV'): to filter one switching pair
# @return
# gene, nPAC, geneTag1, geneTag2, avgUTRlen1, avgUTRlen2, pvalue, padj, change <DE: ...>, PAs1, PAs2
# gene nPAC geneTag1 geneTag2 avgUTRlen1 avgUTRlen2   fisherPV       logFC change               PA1               PA2 dist nDEPA
# 1  PH02Gene00022    2      124      126   218.1048  208.60317 0.51879565 -0.41404606     -1 PH02Gene00022.PA1 PH02Gene00022.PA2  374     2
# 2  PH02Gene00068    2       31       38   351.4839  354.68421 0.75958929  0.55019543      1 PH02Gene00068.PA1 PH02Gene00068.PA2   66     2
UTRtrend_DE<-function(PACds, group, cond1, cond2, avgPACtag=5, avgGeneTag=10,
                      DEpadjThd=0.1, aMovDEPACRes, # to filter DE/DEXPAC
                      mindist=50, fisherThd=0.1, logFCThd=1, selectOne='farest') { #to filter switching pair
  if (is.null(aMovDEPACRes)) {
    stop("aMovDEPACRes must be provided!")
  }

  if (!inherits(aMovDEPACRes, 'movDEPACRes')) stop("UTRtrend_DE::aMovDEPACRes must be of class movDEPACRes")

  if (!(selectOne %in% c('farest','logFC','fisherPV'))) {
    stop("selectOne must be farest|logFC|fisherPV!")
  }

  ds=subsetPACds(PACds, group=group, cond1=cond1, cond2=cond2, avgPACtag=avgPACtag, avgGeneTag=avgGeneTag, avg=TRUE)
  ds=get3UTRAPAds(ds, sortPA=TRUE)

  if (!is.null(aMovDEPACRes)) {
    depas=movSelect(aMovRes=aMovDEPACRes, condpair=c(cond1, cond2), padjThd=DEpadjThd, out='pa')
    if (length(depas)==0) {
      return(NULL)
    }
  }

  #filter gene with DEPA
  if (sum(rownames(ds@anno) %in% depas)==0) {
    cat('PACds rowname is not gene:coord, transforming...')
    paname=paste0(ds@anno$gene,':',ds@anno$coord)
    degenes=ds@anno$gene[paname %in% depas]
    depas=rownames(ds@anno)[paname %in% depas]
    cat(length(depas),'DEPACs\n')
  } else {
    degenes=ds@anno$gene[rownames(ds@anno) %in% depas]
    depas=rownames(ds@anno)[rownames(ds@anno) %in% depas]
  }

  if (length(degenes)==0) {
    return(NULL)
  }

  ds@anno=ds@anno[ds@anno$gene %in% degenes, ]
  ds@counts=ds@counts[rownames(ds@anno),]
  ds@anno$isDEPA=0
  ds@anno[depas,'isDEPA']=1

  ds=subsetPACds(ds,choosePA='apa') #filter APA genes again

  d=cbind(ds@counts,ds@anno[,c('gene','isDEPA','coord')])
  d=cbind(d,PA=rownames(d))

  #For all pairs of PA, if the distance > mindist and either one is DE, then will calculate switching for this pair.
  #From all switching pairs, get on pair as the final pair according to "selectOne".
  #nSwitchPair records the total number of switching pairs.
  #genePAs=d[d$gene=='PH02Gene00121',]
  cat('Get switching for',length(unique(d$gene)),'genes with >=1 DEPAC\n')

  #rt1=plyr::ddply(d1, .(gene), getAllSwitchingPairsForOneGene, mindist=mindist, nDEPAC=1, pseudo=1,
 #                fisherThd=fisherThd, logFCThd=logFCThd, selectOne=selectOne)

  rt=.asDf(d %>% dplyr::group_by(gene) %>% dplyr::do(getAllSwitchingPairsForOneGene(genePAs=., mindist=mindist, nDEPAC=1, pseudo=1,
                                                             fisherThd=fisherThd, logFCThd=logFCThd, selectOne=selectOne)))

  usedpa=unique(c(rt$PA1,rt$PA2))
  ds@anno=ds@anno[usedpa, ]
  ds@counts=ds@counts[usedpa,]
  gstat=geneStatForPACds(ds, avgUTRlen=TRUE, geneTag=TRUE, PAinfo=TRUE, merge=TRUE)
  colnames(gstat)=c('gene','avgUTRlen1', 'avgUTRlen2','nPAC','geneTag1', 'geneTag2','PAs1', 'PAs2')

  if (nrow(rt)!=nrow(gstat)) {
    stop("nrow(rt)!=nrow(gstat)!")
  }
  rt=merge(rt,gstat,by.x='gene',by.y='gene')
  rt=rt[,c('gene', 'nPAC', 'geneTag1', 'geneTag2', 'avgUTRlen1', 'avgUTRlen2', 'fisherPV', 'logFC', 'change', 'PA1', 'PA2', 'dist', 'nDEPA', 'nSwitchPair', 'PAs1', 'PAs2')]
  return(rt)
}


# Similar to UTRtrend_DE, but also for non-3UTR PACs
# @param PACds: PACdataset
# @param group: one column in PACds@colData, "len"
# @param cond1, cond2: 0140, 1200; level in colData of group
# @param avgPACtag=5
# @param avgGeneTag=10
# @param only3UTR=FALSE If TRUE, then first filter 3UTR PACs.
# @param nDEPAC: =0/1/2 switching pair has 1 or 2 DE PACs
# @param DEpadjThd=0.1: to filter DEPAC with padj<0.1
# @param DEXlogFC: to filter DEXPAC with (padj<0.1 and logFC>=DEXlogFC)
# @param aMovDEPACRes movDEPACRes
# @param mindist=50, fisherThd=0.1, logFCThd=1: to filter all switching pair
# @param fisherThd: filter DE list by pvalue<fisherThd
# @param logFCThd: filter |rcij| >=logFCThd, which is Log2(ratio of PA1 between sample A~B / ratio PA2)
#        NOTE: The switching is not neccesary crossing, but only with sufficient difference of the ratio between PA1 and PA2.
#              E.g., it may be that both PA1 and PA2 are used more in sampleA, but PA2 used much more than PA2.
# @param cross FALSE.
# @param selectOne c(NULL-output all pairs, 'farest','logFC','fisherPV')
# @return Dataframe [gene, nPAC, geneTag1, geneTag2, (avgUTRlen1, avgUTRlen2 - if 3UTRPACds), pvalue, padj, change <DE: ...>, PAs1, PAs2, ftr (if not 3UTRPACds)]

APAswitching_DE<-function(PACds, group, cond1, cond2, avgPACtag=5, avgGeneTag=10,
                          only3UTR=FALSE, mergeReps='pool',
                          aMovDEPACRes=NULL, DEpadjThd=0.1, nDEPAC=0, # to filter DE/DEXPAC
                          mindist=50, fisherThd=0.1, logFCThd=1, cross=FALSE, selectOne=NULL) { #to filter switching pair

  if (!is.null(selectOne)) {
    if (!(selectOne %in% c('farest','logFC','fisherPV'))) {
      stop("selectOne must be farest|logFC|fisherPV!")
    }
  }

  if (!(nDEPAC %in% 0:2)) {
    stop("nDEPAC must be 0/1/2!")
  }

  pool=FALSE; avg=FALSE
  if (mergeReps=='pool') pool=TRUE
  if (mergeReps=='avg') avg=TRUE

  ds=subsetPACds(PACds, group, cond1, cond2, avgPACtag=avgPACtag, avgGeneTag=avgGeneTag, noIntergenic=TRUE, avg=avg, pool=pool, verbose=FALSE)
  if (only3UTR) {
    ds=get3UTRAPAds(ds, sortPA=TRUE)
  }

  if (length(ds)==0) {
    return(NULL)
  }

  if (!is.null(aMovDEPACRes) & nDEPAC>0) {
    depas=movSelect(aMovRes=aMovDEPACRes, condpair=c(cond1, cond2), padjThd=DEpadjThd, out='pa')
    if (length(depas)==0) {
      return(NULL)
    }

    #filter gene with DEPA
    if (sum(rownames(ds@anno) %in% depas)==0) {
      cat('PACds rowname may not be gene:coord, try to set new PA name\n')
      paname=paste0(ds@anno$gene,':',ds@anno$coord)
      degenes=ds@anno$gene[paname %in% depas]
      depas=rownames(ds@anno)[paname %in% depas]
    } else {
      degenes=ds@anno$gene[rownames(ds@anno) %in% depas]
      depas=rownames(ds@anno)[rownames(ds@anno) %in% depas]
    }

    if (length(degenes)==0) {
      return(NULL)
    }

    ds@anno=ds@anno[ds@anno$gene %in% degenes, ]
    ds@counts=ds@counts[rownames(ds@anno),]
    ds@anno$isDEPA=0
    ds@anno[depas,'isDEPA']=1

    ds=subsetPACds(ds,choosePA='apa') #filter APA genes again
  } else if (nDEPAC==0) {
    ds@anno$isDEPA=0
  }
  d=cbind(.asDf(ds@counts), ds@anno[,c('gene','isDEPA','coord')])
  d=cbind(d, PA=rownames(d))

  #rt=plyr::ddply(d, .(gene), getAllSwitchingPairsForOneGene, mindist=mindist, nDEPAC=nDEPAC, pseudo=1,
  #               fisherThd=fisherThd, logFCThd=logFCThd, cross=cross, selectOne=selectOne)

  rt=.asDf(d %>% dplyr::group_by(gene) %>% dplyr::do(getAllSwitchingPairsForOneGene(genePAs=., mindist=mindist, nDEPAC=nDEPAC, pseudo=1,
                                                                                     fisherThd=fisherThd, logFCThd=logFCThd, cross=cross, selectOne=selectOne)))


  usedpa=unique(c(rt$PA1,rt$PA2))
  ds@anno=ds@anno[usedpa, ]
  ds@counts=ds@counts[usedpa,]
  avgUTRlen=is3UTRAPAds(ds)
  gstat=geneStatForPACds(ds, avgUTRlen=avgUTRlen, geneTag=TRUE, PAinfo=TRUE, ftr=!avgUTRlen, merge=TRUE)
  if (avgUTRlen) {
    colnames(gstat)=c('gene','avgUTRlen1', 'avgUTRlen2','nPAC','geneTag1', 'geneTag2','PAs1', 'PAs2')
  } else {
    colnames(gstat)=c('gene','nPAC','geneTag1', 'geneTag2','PAs1', 'PAs2','ftr')
  }

  rt=merge(rt,gstat,by.x='gene',by.y='gene')
  if (avgUTRlen) {
    rt=rt[,c('gene', 'nPAC', 'geneTag1', 'geneTag2', 'avgUTRlen1', 'avgUTRlen2', 'fisherPV', 'logFC', 'change', 'PA1', 'PA2', 'dist', 'nDEPA', 'nSwitchPair', 'PAs1', 'PAs2')]
  } else {
    rt=rt[,c('gene', 'nPAC', 'geneTag1', 'geneTag2','fisherPV', 'logFC', 'change', 'PA1', 'PA2', 'dist', 'nDEPA', 'nSwitchPair', 'PAs1', 'PAs2','ftr')]
  }
  return(rt)
}


# -------------- *** APA index **** ---------------

# Get ratio for each PA
# @param chooseAPA TRUE, then will filter APA PACs
# @return matrix as PACds@counts
getAPAratio<-function(PACds, chooseAPA=TRUE) {
  if (chooseAPA) {
    PACds=subsetPACds(PACds, group=NULL, cond1=NULL, cond2=NULL,
                      avgPACtag=0, avgGeneTag=0, choosePA='apa',
                      noIntergenic=TRUE, avg=FALSE, verbose=FALSE)
  }
  gs=geneStatForPACds(PACds, avgUTRlen=FALSE, geneTag=TRUE, PAinfo=FALSE, merge=TRUE)
  gs$nPAC=NULL
  smps=colnames(gs)[-1]
  colnames(gs)[-1]=paste0('gene_',colnames(gs)[-1])

  counts=.asDf(PACds@counts)
  counts$PA=rownames(PACds@counts)
  counts$gene=PACds@anno$gene
  counts=merge(counts,gs,by.x='gene',by.y='gene')
  counts[,smps]=counts[,smps]/counts[,paste0('gene_',smps)]
  rownames(counts)=counts$PA

  return(as.matrix(counts[,smps]))
}


# Get tissue-specificity by shannon
# @param PACds PACdataset
# @param shan.ratio if TRUE then use ratio instead, then only APA PACs are considered.
# @return dataframe of [H (overall TS), Q_min (min TS), Q_min_tissue, Q for each tissue]
#          Smaller H or Q means higher ts.
#          If PA=0 in a sample, then the Q=NA.
# @examples
# ## First avg tissues and remove small PACs, and then get specificity
# pp=subsetPACds(PACds, group='len',cond1=NULL, cond2=NULL, avgPACtag=5, avgGeneTag=10, noIntergenic=TRUE, avg=TRUE, verbose=TRUE)
# head(pp@counts)
# sh=specificityByShannon(pp)
# > head(sh)
# H    Q_min Q_min_cond     0140     0350     0700      1200
# PH02Gene00003.PA1 1.949082 3.534045       0700 3.726690 4.534045 3.534045  4.212117
# PH02Gene00003.PA2 1.906310 3.199092       0700 4.061588 4.521020 3.199092  4.199092
specificityByShannon <- function (PACds, shan.ratio=FALSE) {

  if (shan.ratio) {
    dat=movPAindex(PACds, method='ratio')
    cat('Using ratio for Shannon.\n')
  } else {
    dat=PACds@counts
    cat('Using count for Shannon.\n')
  }
  Ptg=dat/rowSums(dat, na.rm=TRUE)
  Ptg[Ptg==0]=NA  #min(Ptg[Ptg>0])
  Hg=rowSums(-Ptg*log2(Ptg), na.rm=TRUE)
  Qg=Hg-log2(Ptg)

  ##Q_min: the min Q across all tissues; Q_min_tissue min Q for the specified tissue
  minQ=apply(Qg, 1, min, na.rm=TRUE)

  getT<-function(par) {
    m=min(par, na.rm=TRUE)
    paste(names(par)[!is.na(par) & par==m], collapse=',')
  }

  minQt=apply(Qg, 1, getT)
  #minQt=gsub("Q_",'',minQt)

  out=data.frame(H=Hg, Q_min=minQ, Q_min_cond=minQt, as.matrix(Qg))

  Hmax=mean(Hg, na.rm=TRUE)-2*sd(Hg, na.rm=TRUE)
  Qmax=mean(minQ, na.rm=TRUE)-2*sd(minQ, na.rm=TRUE)

  cat('Tissue-specific PAC\'s H_cutoff (mean-2*sd): ',Hmax,'\n')
  cat('Tissue-specific PAC\'s Q_cutoff (mean-2*sd): ',Qmax,'\n')
  cat('Tissue-specific PAC# (H<H_cutoff): ',sum(Hg<Hmax, na.rm=TRUE),'\n')
  cat('Tissue-specific PAC# (Q<Q_cutoff): ',sum(minQ<Qmax, na.rm=TRUE),'\n')

  Hmin=mean(Hg, na.rm=TRUE)+2*sd(Hg, na.rm=TRUE)
  Qmin=mean(minQ, na.rm=TRUE)+2*sd(minQ, na.rm=TRUE)
  cat('Constitutive PAC\'s H_cutoff (mean+2*sd): ',Hmin,'\n')
  cat('Constitutive PAC\'s Q_cutoff (mean+2*sd): ',Qmin,'\n')
  cat('Constitutive PACs (H>H_cutoff): ',sum(Hg>Hmin, na.rm=TRUE),'\n')
  cat('Constitutive PACs (Q>Q_cutoff): ',sum(minQ>Qmin, na.rm=TRUE),'\n')

  return(out)
}


# geometry mean: (A*B*C)^(1/3)
# geometry.mean <- exp(mean(log(x)))
# geo_mean(c(0,0)) --> NaN
# geo_mean(c(10,10)) --> 10
# geo_mean(c(10,1)) --> 3.162278
# geo_mean(c(-10,20)) --> 10 (warning)
# EnvStats::geoMean(c(-10, 20)) --> NA
geo_mean <- function(x) {
  log_data <- log(x)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

# calculate PUI for all peaks, Shulman 2019 et al.
# PUI: deviation from the mean, if peak1 and peak2 is the same in a cluster, then PUI=0.
#      not affected by the scale of the data, which means calc_pui(y/10)=calc_pui(y*10).
# y: a matrix with rownames and colnames; row=peaks, col=clusters (or conditions)
# return: same as y but values are PUI. If a peak is not expressed in a cluster, then PUI=NA.
# ##Ex.
# y <- matrix(c(0, 0, 800, 1000, 50, 10), ncol = 2, dimnames = list(c('peak1','peak2', 'peak3'), c('c1','c2')))
# calc_pui(y)
calc_pui <- function(y, .psudo = 1) {
  peak.counts <- y
  peak.counts[,(colSums(peak.counts) == 0)] <-  NA
  peak.counts <- peak.counts + .psudo
  #utr.geo.mean = apply(peak.counts, 2, EnvStats::geoMean)
  utr.geo.mean = apply(peak.counts, 2, geo_mean)
  i.peak.counts <- y + .psudo
  out <- log2(t(t(i.peak.counts)/utr.geo.mean))
  #colnames(out) <- paste0(colnames(out),"_PUI" )
  out=.asDf(out)
  out
}


# Get geo (PUI) for each PA
# @param chooseAPA TRUE, then will filter APA PACs
# @param psudo=1; first set PAC=NA when the total read count of all PACs in a sample is 0. Then add non-NA+psudo.
# @param sort=TRUE; if false then not sort PACds (for called by movAPAindex to filter proximal Index)
# @return dataFrame with the same format as PACds@counts
## Ex.
# pA=c(0, 11); gm=geo_mean(pA+1); gmPUI=log2((pA+1)/gm)
# --> gm=3.464102; gmPUI=-1.792481  1.792481
getAPAgeo<-function(PACds, chooseAPA=TRUE, psudo=1, sort=TRUE) {
  if (chooseAPA) {
    PACds=subsetPACds(PACds, group=NULL, cond1=NULL, cond2=NULL,
                      avgPACtag=0, avgGeneTag=0, choosePA='apa',
                      noIntergenic=TRUE, avg=FALSE, verbose=FALSE)
  }

  d=cbind(gene=PACds@anno$gene, .asDf(PACds@counts))

  if (sort) d=d[order(d$gene),]

    ## geo mean for each conidtion each gene
    sp=split(d[, -1], d[, 'gene'], drop=T)
    sp=lapply(sp, calc_pui)
    #pui=rbindlist(sp, idcol='gene')
    pui=.asDf(data.table::rbindlist(sp))
    rownames(pui)=rownames(d)

    return(pui)
}

#' Plot distribution of PA index among biological samples.
#'
#' PlotCummPAindex plots distribution of PA index among biological samples by calculating the cumulative fraction.
#'
#' @param PAindex a matrix of PAindex, which is normally from movAPAindex() or movPAindex().
#' @param groupName shown as the legend title.
#' @param xlab x label.
#' @param ylab y label, default is "Cumulative fraction".
#' @return A ggplot2 plot object.
#' @examples
#' data(scPACds)
#' scPACdsCt=subsetPACds(scPACds, group='celltype', pool=TRUE)
#' scPACdsCt=get3UTRAPAds(scPACdsCt, sortPA=TRUE, choose2PA='PD')
#' ## Cummlative plot of GPI index.
#' gpi=movAPAindex(scPACdsCt, method="GPI")
#' plotCummPAindex(PAindex=gpi, groupName='cell type', xlab='GPI')
#' ## The plot of the distribution of Shannon index (tissue-specificity) for each PAC.
#' shan=movPAindex(scPACdsCt, method="shan")
#' plotCummPAindex(PAindex=shan[, -c(1:3)], groupName='cell type', xlab='Shannon index')
#' @name plotCummPAindex
#' @family comparison functions
#' @export
plotCummPAindex<-function(PAindex, groupName='Sample', xlab, ylab="Cumulative fraction") {
  tidy.ppui <- tidyr::gather(data = PAindex)
  colnames(tidy.ppui) <- c("Sample", "value")
  p <- ggplot(data = tidy.ppui, aes_string(x = "value", color = "Sample")) +
    stat_ecdf(size = 1) + theme_bw() + xlab(xlab) +
    ylab(ylab) + labs(color = groupName)
  #+ ggplot2::coord_cartesian(xlim = c(-3, 4.1))
  return(p)
}


# -------------- *** Permutation **** ---------------

# SAAP boostrap
# PAs is a 0,1 vector denoting the tags of PA1 and PA2 for a given sample.
# And then boostrap n times to get number of items with value=lbl.
# @param PAs A vector of 1...1, 2...2 denoting tags under two samples.
# @param lbl The label to be counted
# @param n The time of bootstrap
# @return A vector of numbers.
# NOTE: this function is the same as sampling from 1...totalLength, and then calculate the count <=expNum.
boostrapSa <-function (PAs, lbl=1, n=1000) {
  paNum=c()
  for (i in 1:n) {
    x=sample(PAs, replace = TRUE)
    paNum=c(paNum,sum(x==lbl)) #sampling result under theoretical distribution
  }
  return(paNum)
}

# Z-test to test whether average=mu, result is similar to T-test.
#http://www.endmemo.com/program/R/ztest.php
z.test <- function(a, mu){
  a=a[!is.na(a)]
  #a[is.infinite(a)]=max(a)
  #mu=mean(a)
  sd=sd(a, na.rm=TRUE)
  z = (mean(a) - mu)*sqrt(length(a))/sd(a)
  p <- 2 * pnorm(-abs(z),0,1)
  return(p)
}

## TODO: the pvalue calculation by Ji et al.
## Zc is the boostZ under pv=0.05.
ji.fdr <- function (a, mu) {
  #FDR=(#(|Ze|>Zc)/m)/#(|Zo|>Zc)
  sum(abs(a)>mu)/sum(abs(mu)-mu) #TODO ?? not correct, what is Zc?
}

# Given a 2x2 table with the order raw=c(PA1a=5, PA2a=20, PA1b=10, PA2b=10), get index value
#     Sa	Sb
#PA1	PA1a=5	PA1b=10 -- proximal
#PA2	PA2a=20	PA2b=10 -- distal
# @param pa2x2 The table.
# @param pseudo Add when get RED value, add count first
# @param method RUD=distal/gene, SLR=RED=short/long, RS=|ratioPA1a-ratioPA1b|
# @param avoidInf TRUE. NOTE: If the first calculation returns Inf, then mannually set pseudo=1 and calculate again (even if param pseudo=0).
# @param return An index value.
# @examples
#  getAPAIndexFor2PA(pa2x2=c(5, 0, 5, 0), method='RS', pseudo=0, avoidInf=TRUE)
#  getAPAIndexFor2PA(pa2x2=c(5, 0, 0, 5), method='RS', pseudo=0, avoidInf=TRUE)
#  getAPAIndexFor2PA(pa2x2=c(5, 0, 0, 5), method='SLR', pseudo=0, avoidInf=TRUE)
#  getAPAIndexFor2PA(pa2x2=c(5, 0, 0, 5), method='SLR', pseudo=0, avoidInf=FALSE)
#  getAPAIndexFor2PA(pa2x2=c(5, 0, 0, 5), method='RUD', pseudo=0, avoidInf=FALSE)
#  getAPAIndexFor2PA(pa2x2=c(5, 0, 0, 5), method='RUD', pseudo=0, avoidInf=TRUE)
getAPAIndexFor2PA <- function (pa2x2, method, pseudo=0, avoidInf=TRUE) {
  while(1) {
    pa2x2=pa2x2+pseudo
    if (method=='SLR') {
      idx=log2((pa2x2[1]/pa2x2[2])/(pa2x2[3]/pa2x2[4]))
    } else if (method %in% c('RUD')) {
      idx=log2((pa2x2[2]/(pa2x2[1]+pa2x2[2]))/(pa2x2[4]/(pa2x2[3]+pa2x2[4])))
    } else if (method=='RS') {
      t1=(pa2x2[1]+pa2x2[2])
      t2=(pa2x2[3]+pa2x2[4])
      #idx=(abs(pa2x2[1]/t1-pa2x2[3]/t2)+abs(pa2x2[2]/t1-pa2x2[4]/t2))/2
      idx=abs(pa2x2[1]/t1-pa2x2[3]/t2) #same as above
    }
    if (!is.infinite(idx) | !avoidInf)  return(idx)
    pseudo=1
  }
}


# Get SAAP value for a 2x2 table
# @param raw The 2x2 table
#     Sa	Sb
#PA1	5	10
#PA2	20	10
#raw=c(5, 20, 10, 10)
# @param method RUD, SLR=RED, RS
# @param n Bootstrap times. Hihger n, more likely significant (smaller pv).
# @param pseudo
# @return: c(prox1 dist1 prox2 dist2 realIndex meanPmt sdPmt pv.ztest)
# @examples
# TODO
# Sig. is highly affected by n. The pv changes a lot by different bootstrap runs.
# If there is no swithcing, pv will change also.
# saapOnce(c(5, 20, 10, 10), method='SLR', n=100)
# saapOnce(c(5, 20, 10, 10), method='RS', n=100)
# saapOnce(c(5, 20, 10, 10), method='RUD', n=100)
# saapOnce(raw=c(15, 15, 100, 100), method='RUD', n=1000)
# saapOnce(raw=c(5, 0, 5, 0), method='RS', n=1000)
saapOnce <- function(raw, method='SLR', n=100, pseudo=0) {
  raw=raw+pseudo

  #The real index value
  #If it is NA, then will not do permute, retult=NA
  #If =Inf, then pseudo=1, and then permute to avoid Inf
  rawRED=getAPAIndexFor2PA(raw, method=method, pseudo=0, avoidInf=FALSE)
  if (is.na(rawRED)) return(c(realIndex=rawRED, meanPmt=NaN, sdPmt=NaN, pv.ztest=1))
  if (is.infinite(rawRED)) {
    rawRED=getAPAIndexFor2PA(raw, method=method, pseudo=1, avoidInf=TRUE)
    raw=raw+1
  }

  tota=raw[1]+raw[2]
  totb=raw[3]+raw[4]
  a1=raw[1]/tota; a2=raw[2]/tota;
  b1=raw[3]/totb; b2=raw[4]/totb;
  avgRa=(a1+b1)/2
  avgRb=(a2+b2)/2

  expNumSa=round(avgRa*tota)
  PAsa=c(rep(1,expNumSa), rep(2,tota-expNumSa)) #This is theoretical distribution, NOT raw[1]/raw[2]
  #names(PAsa)=1:length(PAsa)
  PA1a=boostrapSa(PAs=PAsa, n=n) #Under theoretical distribution, boostrapped counts for PA1, n times with n counts.
  PA2a=tota-PA1a

  expNumSb=round(avgRa*totb)
  PAsb=c(rep(1,expNumSb), rep(2,totb-expNumSb))
  PA1b=boostrapSa(PAs=PAsb, n=n)
  PA2b=totb-PA1b

  boostMtx=cbind(PA1a,PA2a,PA1b,PA2b)
  #Get index value according to the boostrap N times results.
  boostRED=apply(boostMtx, 1, getAPAIndexFor2PA, method=method, pseudo=pseudo)

  pvz=z.test(boostRED, mu=rawRED)
  #pvt=t.test(boostRED, mu=rawRED)$p.value #similar to z-test

  #print(head(boostMtx))
  #cat('[',a1,a2,b1,b2,']',avgRa, avgRb,'\n','rawRED=',rawRED,'boostRED=',mean(boostRED),'pvz=',pvz,'pvt=',pvt,'\n')
  return(c(realIndex=rawRED, meanPmt=mean(boostRED, na.rm=TRUE), sdPmt=sd(boostRED, na.rm=TRUE), pv.ztest=pvz))
}

# Bootstrap A' and B' from mixed A and B (one step of GAAP)
# @param ounts The PA count table, the first two columns are sample columns. other columns will be outputed as input.
# @param cols Columns to be boostrap, only two items, can be 1:2 or c('condA','condB') or NULL(1:2)
# @return A new matrix of the same format as counts.
# Note: This is similar to LiuTao's doPermutation, but with simplification.
# @example
# mix1=mixABpermutation(PACds@counts, cols=1:2); colSums(mix1);  colSums(PACds@counts)
mixABpermutation <- function(counts, cols=NULL){

  if (is.null(cols)) cols=colnames(counts)[1:2]
  if (is.numeric(cols)) cols=colnames(counts)[cols]
  if (!(AinB(cols, colnames(counts)))) stop("cols not in colnames(counts)")
  if (length(cols)!=2) stop("cols should be two cols")

  #The start/stop indexes for each PA (each row)
  tot12=rowSums(counts[,cols])
  cumnum12=cumsum(tot12)
  fromTo=.asDf(cbind(start=c(0,cumnum12[-length(cumnum12)])+1,end=cumnum12))

  total1 = sum(counts[,cols[1]])
  total2 = sum(counts[,cols[2]])

  lbls1 = rep(1,total1)
  lbls2 = rep(0,total2)
  label = c(lbls1,lbls2)
  shuffleLbl = sample(label, length(label), replace = FALSE)

  #Here, only need the tottags of each PA under s1 and s2 for counting??
  #Do not need distinguish the lables for cond1 and cond2.
  #For example, PA1@cond1=16, @cond2=20, then after the 01 array of the gaap shuffle, count the number with labels=0 or 1 for the first 36 tags.
  #For each PA, get tag number with label=1 as the count for cond1; the count for cond2 is the rowsum - cond1#
  count1=apply(fromTo, 1, function(id) {
    if (id[2]<id[1]) return(0)
    return(sum(shuffleLbl[id[1]:id[2]]))}
  )
  count2=rowSums(counts[,1:2])-count1
  counts[,cols[1]]=count1
  counts[,cols[2]]=count2
  return(counts)
}

# Sampling for a lib, bootstrap each column respectively.
# If libsize<tot_tag of a sample, then sample without replace, otherwise with replace.
# @param counts PA count matrix, each column is one sample.
# @param cols Sample columns to be used, e.g., 1:2:3 or c('condA','condB') or NULL(use all columns in counts)
# @param libsize The Libsize after sampling.
# @param New matrix with the same format as counts
# NOTE: Similar to LiuTao's bootStrap with simplification.
# @example
# sub1=subSampleLibrary(PACds@counts, 1:2, libsize=1000); colSums(sub1)[1:2]
subSampleLibrary <- function(counts, cols=NULL, libsize=1e6){
  if (is.null(cols)) cols=colnames(counts)
  if (is.numeric(cols)) cols=colnames(counts)[cols]
  if (!(AinB(cols, colnames(counts)))) stop("cols not in colnames(counts)")

  for (col in cols) {
    cond1=cumsum(counts[,col])
    fromTo=.asDf(cbind(start=c(0,cond1[-length(cond1)])+1,end=cond1))

    tot1 = sum(counts[,col])
    if (tot1>libsize) {
      rand1=sample(1:tot1, libsize, replace = FALSE)
    } else {
      rand1 = sample(1:tot1, libsize, replace = TRUE)
    }
    count1=apply(fromTo, 1, function(id) {sum(rand1<=id[2] & rand1>=id[1])})
    counts[, col] = count1
  }
  return(counts)
}



gaap <- function(PACds) {
  # TODO... I feel that the amount of calculation is large, and need to compare each pair.

  #For two samples, mix and permuate
  mix1=mixABpermutation(PACds@counts, cols=1:2)

  #Then subsampled same number of tags from the two samples
  sub1=subSampleLibrary(PACds@counts, 1:2, libsize=1000)

  #Then get APAindex, DE... for this bootstrapped PACds, and get DE number

  #Do analyses for the real data

  #Get corrected DE number

}



# -------------- *** Stat output **** ---------------

# Get heatmap results from DEXSeqResults (either for single columns to get exon usage, or to get log2FC for a pair)
# - object
# - conds: c('0140','0350')
# - condPairs: matrix of two columns, cbind(c('0140','0350','1200'),c('0350','0700','0140'))
# - padjThd=0.1:
# - valueThd=NULL (padj and value filtered by getDEXlist())
# - return: heatmapResult type = list[padj=df[conds or cond1.cond2], value=df]

# Ex1. per sample exon usage
# heatDxr=getHeatMtx_DEXSeqResults(dxr, conds=c('0140','0350','1200'), padjThd=0.1, valueThd=1)
# nrow(dxr); nrow(heatDxr@padj); head(heatDxr@padj); head(heatDxr@value)
# 0140     0350      1200
# PH02Gene02460:4252966  8.236215 7.920872  7.854877
# PH02Gene02460:4253337  4.099735 4.603733  4.712841

# Ex2. pair sample log2FC
# heatDxr=getHeatMtx_DEXSeqResults(dxr, condPairs=cbind(c('0140','0350','1200'),c('0350','0700','0140')), padjThd=0.1, valueThd=1)
# 0140.0350  0350.0700  1200.0140
# PH02Gene02465:4326528   0.4194957 -0.9128486 -2.9909194
# PH02Gene02465:4331578  -0.2987850  0.6501746  2.1302765
getHeatMtx_DEXSeqResults <- function (object, conds=NULL, condPairs=NULL, padjThd=0.1, valueThd=NULL) {
  if (!inherits(object, "DEXSeqResults")) stop ("object's class must be DEXSeqResults")
  cols=colnames(object)
  oheat=list()
  object=.asDf(object)

  # Given conds, then get exon usage per sample
  # Filtering by padj<padjThd and max(cond)>=valueThd (At least one sample meets the filteration)
  if (!is.null(conds)) {
    lbls=conds
    oheat$padj=matrix(nrow=0, ncol=length(conds),dimnames = list(c(),conds))
    oheat$value=matrix(nrow=0, ncol=length(conds),dimnames = list(c(),conds))
    if (length(conds)<=1) stop("conds must at least 2 elements")
    conds=make.names(conds)
    if (sum(!(conds %in% cols))!=0) stop("conds not all in colnames(object)")
    idx=!is.na(object$padj) & object$padj<padjThd
    if (sum(idx)==0) return(oheat)
    object=object[idx,]
    if (!is.null(valueThd)) {
      minv=apply(object[,conds], 1, max, na.rm=TRUE)
      idx=which(!is.na(minv) & minv>=valueThd)
      if (length(idx)==0) return(oheat)
      object=object[idx,]
    }
    oheat$padj=matrix(rep(object$padj,length(conds)),
                      nrow=nrow(object),ncol=length(conds),
                      dimnames = list(rownames(object),lbls))
    oheat$padj=.asDf(oheat$padj)
    oheat$value=object[,conds]
    colnames(oheat$value)=lbls

    oheat=new("heatmapResults", padj=oheat$padj, value=oheat$value)
    return(oheat)
  }

  # Given condpairs, then filter columns with log2fold_
  if (!is.null(condPairs)) {

    if (is.null(rownames(condPairs))) {
      lbls=paste0(condPairs[,1], '.',condPairs[,2])
    } else {
      lbls=rownames(condPairs)
    }

    oheat$padj=matrix(1,
                      nrow=nrow(object),ncol=nrow(condPairs),
                      dimnames = list(rownames(object),lbls))
    oheat$padj=.asDf(oheat$padj)
    oheat$value=oheat$padj

    .validCondPairs(condPairs)
    logcols1=paste0('log2fold_',condPairs[,1], '_',condPairs[,2])
    logcols2=paste0('log2fold_',condPairs[,2],'_', condPairs[,1])
    #If in dxr it is log_B_A, then value*(-1) and change the column name as log_A_B (same as logcols1)
    for (i in 1:length(logcols1)) {
      if (!(logcols1[i] %in% cols)) {
        if (logcols2[i] %in% cols) {
          object[,logcols2[i]]=object[,logcols2[i]]*(-1)
          cols[which(cols==logcols2[i])]=logcols1[i]
        } else {
          stop(cat(logcols1[i],'&',logcols2[i],'not in object, use estimateExonFoldChanges to get log2FC\n'))
        }
      }

      #DE for each pair of conditions
      depac=getDEXlist(dxr=object, cond1=condPairs[i,1], cond2=condPairs[i,2], padj=padjThd, logFC=valueThd, full=FALSE)
      #If not DE, then set padj=1
      if (length(depac)>0) {
        oheat$padj[depac,lbls[i]]=object[depac,'padj']
      }
      oheat$value[,i]=object[,logcols2[i]]
    }
    #Remove rows with all padj=1
    idx=which(rowSums(oheat$padj==1)==length(lbls))
    if (length(idx)>0) {
      oheat$padj=oheat$padj[-idx,,drop=FALSE]
      oheat$value=oheat$value[-idx,,drop=FALSE]
    }
    oheat=new("heatmapResults", padj=oheat$padj, value=oheat$value)
    return(oheat)
  }

}


# Get heatmap result from DEAllResults
# if in log2FC, it is B_vs_A but in condPair it is A_vs_B, then log2FC value will be (*-1)
# @param object
# @param condPairs Matrix of two columns, "c1_vs_c2", cbind(c('0140','0350','1200'),c('0350','0700','0140'))
#        If there are rownames, then use rownames for the final column names.
# @param padjThd=0.1: filtered by object$padj
# @param valueThd: filtered by max(condPairs) abs(log2FC) >=valueThd
# @return: heatmapResult type = list[padj=df[conds or cond1.cond2], value=df]
# NOTE: if object is an interaction DEAllResults, then a condpair name will have the level of another factor,
#       e.g., len0140_vs_0350.posdn, then the final output columns is 0140.0350.posdn

# Ex1. pair sample log2FC
# heatDE=getHeatMtx_DEAllResults(resLRT_heats, condPairs=cbind(c('0140','0350','1200'),c('0350','0700','0140')), padjThd=0.1, valueThd=1)
# nrow(heatDE@padj); head(heatDE@value); head(resLRT_heats@pairwise$log2FC_len.pos)
# 0140.0350.posmd 0140.0350.posdn 0140.0350.posup 0350.0700.posmd 0350.0700.posdn
# PH02Gene00003             0.6438561            -0.5000734            0.50199196            -2.2240011            -1.1520032
# PH02Gene00004            -0.2630228             0.5849607            0.20022273            -5.3633977            -3.5235655

getHeatMtx_DEAllResults <- function (object, condPairs, padjThd=0.1, valueThd=NULL) {

  if (!inherits(object, "DEAllResults")) stop("object is not class DEAllResults")

  # object
  # > colnames(DEAllRes$ref$log2FC_len) (~len)
  # [1] "len0140_vs_1200" "len0350_vs_1200" "len0700_vs_1200"

  # colnames(resLRT_heats@series$log2FC_len.pos) ## _len.pos (~len*/+pos)
  # "len0140_vs_0700.posmd" "len0350_vs_0140.posmd"

  .validCondPairs(condPairs)

  cols1=paste0(condPairs[,1], '_vs_',condPairs[,2])
  cols2=paste0(condPairs[,2],'_vs_', condPairs[,1])

  #Return the index of each cols1 in allcols; if the order is flipped, then the index is -index.
  #.colInWhere(cols1,cols2,allcols=colnames(resLRT_heats@pairwise[[1]]))
  # return:
  #$`0140_vs_0350`
  #[1]  -4 -10 -16
  #> colnames(resLRT_heats@pairwise[[1]])[c(4,10,16)]
  #[1] "len0350_vs_0140.posmd" "len0350_vs_0140.posdn" "len0350_vs_0140.posup"
  #$`0350_vs_0700`
  #[1]  2  8 14
  .colInWhere<-function(cols1, cols2, allcols) {
    ocols=list()
    for (ci in cols1) {
      gi=grep(ci,allcols,fixed=TRUE)
      if (length(gi)>0) {
        ocols[[ci]]=gi
      } else {
        gi=grep(cols2[which(cols1==ci)],allcols,fixed=TRUE)
        if (length(gi)>0) {
          ocols[[ci]]=gi*(-1)
        }
      }
    }
    return(ocols)
  }

  ocols=c()
  where='pairwise'
  if (length(object@pairwise)>0) {
    ocols=.colInWhere(cols1,cols2,colnames(object@pairwise[[1]]))
    where='pairwise'
  }

  if (length(ocols)!=length(cols1)) {
    if (length(object@ref)>0) {
      ocols=.colInWhere(cols1,cols2,colnames(object@ref[[1]]))
      where='ref'
    }
  }

  if (length(ocols)!=length(cols1)) {
    if (length(object@series)>0) {
      ocols=.colInWhere(cols1,cols2,colnames(object@series[[1]]))
      where='series'
    }
  }

  if (length(ocols)!=length(cols1)) {
    stop(cat("object do not have",cols1,'\n'))
  }

  padj=slot(object,where)[[grep('padj',names(slot(object,where)))]]
  value=slot(object,where)[[grep('log2FC',names(slot(object,where)))]]

  #if the index of ocols is negative, then flip the cols and set log2FC*-1
  for (i in 1:length(ocols)) {
    if (sum(ocols[[i]]<0)>0) {
      ocols[[i]]=abs(ocols[[i]])
      colnames(padj)[ocols[[i]]]=gsub(cols2[i],cols1[i],colnames(padj)[ocols[[i]]],fixed=TRUE)
      colnames(value)[ocols[[i]]]=gsub(cols2[i],cols1[i],colnames(value)[ocols[[i]]],fixed=TRUE)
      value[,ocols[[i]]]=value[,ocols[[i]]]*(-1)
    }
    #len0350_vs_0140.posdn --> 0350.0140.posdn
    if (length(ocols[[i]])>1) {
      cnames=colnames(value)[ocols[[i]]]
      cnames=strsplit(cnames,'\\.')
      c2=unlist(lapply(cnames,'[',2))
      c1=gsub('_vs_','\\.',names(ocols)[i])
      colnames(value)[ocols[[i]]]=paste0(c1,'.',c2)
      colnames(padj)[ocols[[i]]]=paste0(c1,'.',c2)
    } else {
      colnames(value)[ocols[[i]]]=paste0(condPairs[i,1],'.',condPairs[i,2])
      colnames(padj)[ocols[[i]]]=paste0(condPairs[i,1],'.',condPairs[i,2])
    }
  }
  padj=padj[,abs(unlist(ocols)),drop=FALSE]
  value=value[,abs(unlist(ocols)),drop=FALSE]

  if (!is.null(rownames(condPairs))) {
    for (i in 1:nrow(condPairs)) {
      colnames(padj)=gsub(paste0(condPairs[i,1],'.',condPairs[i,2]),rownames(condPairs)[i],colnames(padj),fixed=TRUE)
    }
  }
  colnames(value)=colnames(padj)

  #filter
  value[value==Inf]=max(value[value!=Inf])
  value[value==-Inf]=min(value[value!=-Inf])
  idx=rep(TRUE, length(padj))
  if (!is.null(padjThd)) {
    idx=!is.na(padj) & padj<padjThd
  }
  if (!is.null(valueThd)) idx=idx & !is.na(value) & abs(value)>=valueThd
  idx=which(rowSums(idx)>=1)
  padj=padj[idx, ,drop=FALSE]
  value=value[idx,,drop=FALSE]

  oheat=new("heatmapResults", padj=padj, value=value)
  return(oheat)
}

#Validate condPairs
#At least two rows, and without rows of c1.c2 = c2.c1, e.g., (140,350),(350,140)
.validCondPairs<-function(condPairs) {
  #if (nrow(condPairs)<=1) stop("condPairs at least 2 rows")
  c1=paste0(condPairs[,1],'.',condPairs[,2])
  c2=paste0(condPairs[,2],'.',condPairs[,1])
  if (length(base::intersect(c1,c2))>0) stop("condPairs overlapping!")
}

# Get heatmapRes, fisherPV/log2FC (method=de/dex) or padj/cor (method=lineartrend), from swResults
# @param object
# @param condPairs Matrix of two columns, "c1_vs_c2", cbind(c('0140','0350','1200'),c('0350','0700','0140'))
#        If there are rownames, then use rownames for the final column names.
# @param padjThd=0.1: filtered by max cor (lineartrend) | logFC (APAswitching) of (condPairs)>=valueThd
# @param valueThd: filtered by max(condPairs) abs(log2FC) >=valueThd
# @param padjCol: NULL (autoset logFC/cor)
# @param valueCol: NULL (autoset padj/fisherPV), user can also provide the column names for padj and value.
# @param fixed If TRUE, then all rows in condpairs and the order of cond1.cond2 should be in object@sw!
# @return: heatmapResult type = list[padj=df[conds or cond1.cond2], value=df]
# NOTE: Because the row number of different condpairs in swRes are not the same, so first filter genes.
#       If the condname is flipped, then change the name of sw and set logFC/cor*-1.
#       Fianlly, build an NA matrix, and fill genes from respective condpairs, there will be NAs in the heatmapResult.
# @examples
#heatUtrDex=getHeatMtx_swResults(object=utrDEX, condPairs=cbind(c('0140','0350','1200'),c('0350','0700','0140')), padjThd=1, valueThd=NULL)
#heatUtr=getHeatMtx_swResults(object=utr, condPairs=cbind(c('0140','0350','1200'),c('0350','0700','0140')), padjThd=1, valueThd=NULL)
#heatSWDex=getHeatMtx_swResults(object=swDEX, condPairs=cbind(c('0140','0350','1200'),c('0350','0700','0140')), padjThd=1, valueThd=NULL)
# Ex. rename output columns:
#condPairs=matrix(cbind(c('0140','0350','1200'),c('0350','0700','0140')),nrow=3, dimnames=list(c('x','y','z'),c('c1','c2')))
#rename=getHeatMtx_swResults(object=swDEX, condPairs, padjThd=1, valueThd=NULL)
# x           y          z
# PH02Gene02437 0.67409074 1.000000000         NA
# PH02Gene02445 0.41403561 0.830806742 0.85470880
# PH02Gene02452 0.76943375 1.000000000 0.60158103
getHeatMtx_swResults <- function(object, condPairs, padjThd=0.1, valueThd=NULL, padjCol=NULL, valueCol=NULL, fixed=TRUE) {

  # Given swResults[[i]] automatelly filter the list meets the cretieria
  # Return gene, PA, or index
  # @param filters=list() given columns in swRes, such as padj/value/change, and judge to use < or >= or = for filtering.
  #   filters list can have many elements, but the used elements should have words like 'padj','pvalue','fisherPV' | 'logFC','cor','value' | change .
  # @param padj: lineartrend
  # @param fisherPV/logFC: APA switching
  # @param change: select change
  # @return gene, PA, index(row idx) (get from PA1/PA2 for APA switching; or PAs1/PAs2 for lineartrend)
  .filterOneSWRes<-function(sw, filters=list(padj=NULL, fisherPV=NULL, logFC=NULL, cor=NULL, change=NULL), out) {

    if (!(out %in% c('gene','PA','index'))) {
      stop("out must be gene or PA or index")
    }
    cols=colnames(sw)
    if (out=='gene' & !('gene' %in% cols)) stop("out=gene, but gene not in sw")

    idx=rep(TRUE,nrow(sw))

    for (n in names(filters)) {
      if (!is.null(filters[[n]])) {
        if (!(n %in% cols))  stop(cat(n,'!=NULL, but it is not in \n'))
      } else {
        next
      }

      thd=filters[[n]]
      if (n %in% c('padj','pvalue','fisherPV','pv.ztest')) idx=idx & (!is.na(sw[,n]) & sw[,n]<thd)
      if (n %in% c('logFC','cor','value','realIndex')) idx=idx & (!is.na(sw[,n]) & sw[,n]>=thd)
      if (n %in% c('change')) idx=idx & (!is.na(sw[,n]) & sw[,n]==thd)
    }

    if (out=='index') return(which(idx==TRUE))
    if (out=='gene') return(unique(sw$gene[idx]))

    if ('PA1' %in% cols & 'PA2' %in% cols) {
      PAs=unique(c(sw$PA1[idx],sw$PA2[idx]))
    } else if ('PAs1' %in% cols & 'PAs1' %in% cols) {
      PAs=unique(c(sw$PAs1[idx],sw$PAs2[idx])) #"PH02Gene00202.PA3=92;PH02Gene00202.PA4=42"
      PAs=unlist(strsplit(PAs,';'))
      PAs=gsub('=.*','',PAs)
    } else {
      stop("out=PA, but PA1/2 or PAs1/2 not in sw cols")
    }
    return(unique(PAs))
  }

  METHOD_DE=c('de','dex','deseq','dexseq','deseq2','chisq')
  METHOD_LINEAR=c('lineartrend','linear')

  if (!inherits(object,'swResults')) stop("object is not class swResults")
  .validCondPairs(condPairs)
  lbls=rownames(condPairs)
  rownames(condPairs)=paste0(condPairs[,1],'.',condPairs[,2])

  if (fixed) {
    if (!(AinB(rownames(condPairs), names(object@sw)))) stop("fixed=TRUE, but condPairs not exact as in object@sw")
  }

  object@method=tolower(object@method)
  if (is.null(padjCol))  {
    if (object@method %in% METHOD_DE) {
      padjCol='fisherPV'
    } else if (object@method %in% METHOD_LINEAR) {
      padjCol='padj'
    } else {
      stop("padjCol=NULL, but object$method not de/dex/deseq/dexseq/deseq2/lineartrend/linear")
    }
  }

  if (is.null(valueCol))  {
    if (object@method %in% METHOD_DE) {
      valueCol='logFC'
    } else if (object@method %in% METHOD_LINEAR) {
      valueCol='cor'
    } else {
      stop("valueCol=NULL, but object$method not de/dex/deseq/dexseq/deseq2/lineartrend/linear")
    }
  }


  filters=list(padj=NULL, fisherPV=NULL, logFC=NULL, change=NULL)
  filters[[padjCol]]=padjThd
  filters[[valueCol]]=valueThd

  # NOTE: Because the row number of different condpairs in swRes are not the same, so first filter genes.
  #       If the condname is flipped, then change the name of sw and set logFC/cor*-1.
  #       Fianlly, build an NA matrix, and fill genes from respective condpairs, there will be NAs in the heatmapResult.
  genes=c()
  for (i in 1:nrow(condPairs)) {
    condpair=rownames(condPairs)[i]
    if (!(condpair %in% names(object@sw))) {
      condpair=paste0(condPairs[i, 2],'.',condPairs[i, 1])
      if  (!(condpair %in% names(object@sw))) {
        stop(cat(condpair,'not in object',names(object),'\n'))
      }
    }
    #APAswitching may have many switching pairs, randomly choose 1 pair.
    dup=duplicated(object@sw[[condpair]]$gene)
    if (sum(dup)>0) {
      object@sw[[condpair]]=object@sw[[condpair]][!dup,]
    }

    g=.filterOneSWRes(object@sw[[condpair]], filters=filters, out="gene")
    if (rownames(condPairs)[i]!=condpair) {
      object@sw[[condpair]][,valueCol]=(-1)*object@sw[[condpair]][,valueCol]
      if ('change' %in% cols(object@sw[[condpair]])) {
        object@sw[[condpair]][,'change']=(-1)*object@sw[[condpair]][,'change']
      }
      names(object@sw)[which(names(object@sw)==condpair)]=rownames(condPairs)[i]
    }
    genes=c(genes,g)
  }

  genes=unique(genes)
  padj=matrix(NA,nrow=length(genes),ncol=nrow(condPairs))
  rownames(padj)=genes
  if (!is.null(lbls)) {
    colnames(padj)=lbls
  } else {
    colnames(padj)=rownames(condPairs)
  }
  padj=.asDf(padj)
  value=padj

  for (i in 1:nrow(condPairs)) {
    sw=object@sw[[rownames(condPairs)[i]]]
    mt=match(genes,sw$gene)
    if (sum(!is.na(mt))==0) next
    padj[!is.na(mt),i]=sw[,padjCol][mt[!is.na(mt)]]
    value[!is.na(mt),i]=sw[,valueCol][mt[!is.na(mt)]]
  }

  oheat=new("heatmapResults", padj=padj, value=value)
  return(oheat)
}


#Get heatmapResults from a list of betas (which is calculated from DESeqLRTAll())
# @param object: betas list from DESeqLRTAll()
# @param condPairs: two columns with the 2nd column the ref level
# @param padjThd: filter object$padj
# @param valueThd: filter object$betas
# @return heatmapResults
# @example
#condPairs=matrix(cbind(c('0350','0700','1200'),c('0140','0140','0140')),nrow=3, dimnames=list(c('x','y','z'),c('c1','c2')))
#heatBetas=getHeatMtx_betas(object=ccmtx, condPairs, padjThd=1, valueThd=NULL)
getHeatMtx_betas <- function (object, condPairs, padjThd=0.1, valueThd=NULL) {
  # > colnames(ccmtx$betas)
  # [1] "len_0350_vs_0140.posdn"      "len_0700_vs_0140.posdn"      "len_1200_vs_0140.posdn"      "pos_md_vs_dn.len0140"
  # [5] "pos_up_vs_dn.len0140"        "len0350_vs_0140.posmd_vs_dn" "len0700_vs_0140.posmd_vs_dn" "len1200_vs_0140.posmd_vs_dn"
  # [9] "len0350_vs_0140.posup_vs_dn" "len0700_vs_0140.posup_vs_dn" "len1200_vs_0140.posup_vs_dn" "padj"

  if (!inherits(object, 'list')) stop("object is not class list")
  .validCondPairs(condPairs)
  ref=unique(condPairs[,2])
  if (length(ref)!=1) stop("condPairs[,2] is not the same (ref level)")
  cnames=colnames(object$betas)
  ii=c()
  for (i in 1:nrow(condPairs)) {
    vs=sprintf("%s_vs_%s",condPairs[i,1],ref)
    idx=grep(vs,cnames)
    if (length(idx)==0) stop(cat(vs,'not in colnames of object$betas\n'))
    ii=c(ii,idx)
    orgnames=strsplit(cnames[idx],'\\.')
    c2=unlist(lapply(orgnames,'[',2))
    cnames[idx]=paste0(sprintf("%s.%s",condPairs[i,1],ref),'.',c2)
  }
  colnames(object$betas)=cnames
  colnames(object$padj)=cnames[-length(cnames)]
  value=object$betas[,ii]
  padj=object$padj[,ii]

  #filter
  value[value==Inf]=max(value[value!=Inf])
  value[value==-Inf]=min(value[value!=-Inf])
  idx=!is.na(padj) & padj<padjThd
  if (!is.null(valueThd)) idx=idx & !is.na(value) & abs(value)>=valueThd
  idx=which(rowSums(idx)>=1)
  padj=padj[idx,]
  value=value[idx,]

  if (!is.null(rownames(condPairs))) {
    for (i in 1:nrow(condPairs)) {
      colnames(padj)=gsub(paste0(condPairs[i,1],'.',condPairs[i,2]),rownames(condPairs)[i],colnames(padj),fixed=TRUE)
    }
  }
  colnames(value)=colnames(padj)

  oheat=new("heatmapResults", padj=padj, value=value)
  return(oheat)

}

# Get heatmapResults from different classes of objects.
# Wrapper function to get heatmapData from different classes
# @param object DEXSeqResults/DEAllResults/swResults/list(betas)
# @param colData One column for exonUsage-DEXSeqResults or two columns with >=2 rows, must with rownames=lbl, columns=cond or cond1/cond2
# @param padjThd/valueThd To filter sig. object
# @return A heatmapResults object [padj, value, colData]
# @examples
#colData=matrix(cbind(c('0140','0140','0140'),c('0350','0700','1200')), nrow=3, dimnames=list(c('350','700','1200')))
#getHeatmapRes(object=utr, colData=colData, padjThd=0.1, valueThd=NULL)
#getHeatmapRes(object=utrDEX, colData=colData, padjThd=0.1, valueThd=NULL)
#getHeatmapRes(object=swDEX, colData=colData, padjThd=0.1, valueThd=NULL)
#getHeatmapRes(object=dxr, colData=colData, padjThd=0.1, valueThd=NULL)
#getHeatmapRes(object=ddsNoInt_heats, colData=colData, padjThd=0.1, valueThd=NULL)
#heat=getHeatmapRes(object=resLRT_heats, colData=colData, padjThd=0.1, valueThd=NULL)
getHeatmapRes <- function(object, colData, padjThd=0.1, valueThd=NULL) {

  if (!(class(object) %in% c('DEXSeqResults','DEAllResults','swResults','list'))) {
    stop('class of object must be DEXSeqResults/DEAllResults/swResults/list')
  }

  if(is.list(object)) {
    if (!('betas' %in% names(object))) stop("object is list (change-change), but betas not in object")
  }

  if (is.null(rownames(colData))) stop("rownames of colData is null")
  #if (nrow(colData)<=1) stop("getHeatmapRes: nrow of colData must >=2")
  if (ncol(colData)==1 & (class(object) %in% c('DEAllResults','swResults'))) stop("colData must =2 cols for DEAllResults/swResults")

  if (inherits(object, 'DEXSeqResults')) { #log2fold_1200_0350
    if (ncol(colData)==1) {
      heat=getHeatMtx_DEXSeqResults(object, conds=colData, condPairs=NULL, padjThd=padjThd, valueThd=valueThd)
    } else {
      heat=getHeatMtx_DEXSeqResults(object, conds=NULL, condPairs=colData, padjThd=padjThd, valueThd=valueThd)
    }

  } else if (inherits(object, 'DEAllResults')) {
    heat=getHeatMtx_DEAllResults(object, condPairs=colData, padjThd=padjThd, valueThd=valueThd)
  } else if (inherits(object, 'swResults')) {
    heat=getHeatMtx_swResults(object, condPairs=colData, padjThd=padjThd, valueThd=valueThd)
  } else if (inherits(object, 'list')) {
    heat=getHeatMtx_betas(object, condPairs=colData, padjThd=padjThd, valueThd=valueThd)
  }
  heat@colData=.asDf(colData)
  return(heat)
}


# Make stat of heatmapResults object
# @param heat heatmapResults (padj, value, colData)
# @param ... filter params
# @param lbl For the label of nsig, e.g., all, upregulate, dn
# @return - list[nsig, tf01, <<<de01, deNum; if it is for all condpairs>>>]
#           nsig (significant DE number)
#           siglist (list of sig genes)
#           ovp (dataframe of overlapping stats, <pair, n1.all, n2.all, novp.all)
#           tf01 (01matrix)
#           de01 (01matrix, each column is one condition) - only for all condpairs
#           deNum(matrix, each column is a condition denoting the DE pair number) - only for all condpairs
# 2019-04-02 bug fixed, even heat only have 1 column, will have output, but without ovp, tf01, de01, deNum
statGeneralHeat <- function(heat, padjThd=NULL, valueThd=NULL, upThd=NULL, dnThd=NULL, lbl=NULL) {
  ol=list()
  tf=subsetHeatmap(heat, padjThd=padjThd, valueThd=valueThd, upThd=upThd, dnThd=dnThd, returnTF=TRUE)

  #nSig: DE number for each column
  nsig=colSums(tf)
  nsig=data.frame(num=nsig)
  rownames(nsig)=colnames(tf)
  if (!(is.null(lbl))) colnames(nsig)=lbl
  ol[['nsig']]=nsig

  # list for each column (row names, normally is gene)
  siglist=lapply(tf, function(par, rows) {rows[par]}, rows=rownames(tf))
  ol[['siglist']]=siglist

  # if heatmap only has one column
  if (length(siglist)<=1) {
    return(ol)
  }

  # tf1: the DE status of each gene among columns
  tf01=apply(tf,2,as.numeric) #pheatmap(tf01)
  if (nrow(tf)==1) {
    tf01=matrix(tf01, nrow=1)
    colnames(tf01)=colnames(tf)
  }

  if (length(tf01)==0) {
    cat("statGeneralHeat: no result.\n")
    return(invisible(NULL))
  }

  rownames(tf01)=rownames(tf)
  tf01=.asDf(tf01)
  ol[['tf01']]=tf01

  # overlapping stat
  ovp=matrix(nrow=0,ncol=4,dimnames = list(c(),c('pair','n1','n2','novp')))
  if (!is.null(lbl)) colnames(ovp)[2:4]=paste0(colnames(ovp)[2:4],'.',lbl)
  for (i in 1:(length(siglist)-1)) {
    for (j in (i+1):length(siglist)) {
      n1=length(siglist[[i]])
      n2=length(siglist[[j]])
      novp=length(base::intersect(siglist[[i]],siglist[[j]]))
      ovp=rbind(ovp,c(paste(names(siglist)[c(i,j)],collapse = '-'),n1,n2,novp))
    }
  }
  ovp=.asDf(ovp)
  ovp[,2:4]=lapply(ovp[,2:4],as.numeric)
  ol[['ovp']]=ovp

  # deMtx, deMtx2
  # The intensity of each gene in each condition (not cond pair) by counting the number of DE pairs involving the cond
  # Only if there are all condpairs for all conds. E.g., A-B A-C B-C; but A-B A-C can not because A two times but B/C only 1 time.
  if (ncol(heat@colData)==2 &
      nrow(heat@colData)==ncol(tf) &
      length(unique(table(unlist(heat@colData))))==1) {
    cat('All cond pairs in heat@colData, get de01 and deNum\n')
    rawConds=unique(unlist(heat@colData))
    conds=make.names(rawConds)
    heat@colData[]=apply(heat@colData,2,make.names)
    deNum=data.frame(matrix(0,nrow=nrow(tf),ncol=length(conds),dimnames = list(rownames(tf),conds)))
    for (i in 1:nrow(heat@colData)) {
      deNum[,unlist(heat@colData[i,])]=deNum[,unlist(heat@colData[i,])]+tf[,rownames(heat@colData)[i]]
    }
    de01=deNum
    de01[de01>0]=1  #0/1 matrix
    colnames(de01)=rawConds
    colnames(deNum)=rawConds
    # de01                       0140 0350 0700 1200
    # PH02Gene02460:4253337      0     1     2     1
    # PH02Gene02465:4326528      0     1     2     1
    # PH02Gene02465:4331578      0     1     2     1

    ol[['de01']]=de01
    ol[['deNum']]=deNum
  }
  return(ol)
}


#' Subset heatmapResults object
#'
#' heatmapResults is a generic class to store differential analyses results, which is part of a movRes object.
#' subsetHeatmap is to filter a heatmapResults by padj and/or value thresholds.
#'
#'
#' @param heat a heatmapResults object
#' @param padjThd a cutoff for padj
#' @param valueThd a cutoff for value, to filter by abs(..)>=valueThd
#' @param upThd a cutoff to filter value by ...>=upThd (up-regulated)
#' @param dnThd a cutoff to filter by value by ...<=dnThd (down-regulated)
#' @param returnTF If TRUE, then return a matrix with TRUE or FALSE, otherwise return a heatmapResults object.
#' @export
subsetHeatmap <- function(heat, padjThd=0.1, valueThd=NULL, upThd=NULL, dnThd=NULL, returnTF=FALSE) {
  if ((!is.null(valueThd) | !is.null(upThd)) & (!is.null(upThd) | !is.null(dnThd))  & (!is.null(valueThd) | !is.null(dnThd)) ) stop("valueThd, upThd, dnThd, only one shoud be not NULL")

  heat@value[is.na(heat@value)]=0

  if (nrow(heat@padj)==0) {
    if (returnTF) return (heat@padj)
    return(heat)
  }

  tf=heat@padj
  tf[]=TRUE
  if (!is.null(padjThd)) {
    tf=!is.na(heat@padj)
    tf=tf & heat@padj<padjThd
  }
  if (!is.null(upThd)) {
    tf=tf & !is.na(heat@value) & heat@value>=upThd
  }
  if (!is.null(dnThd)) {
    tf=tf & !is.na(heat@value) & heat@value<=dnThd
  }
  if (!is.null(valueThd)) {
    tf=tf & !is.na(heat@value) & abs(heat@value)>=valueThd
  }
  rid=rowSums(tf)>0
  if (returnTF) {
    tf=tf[rid,,drop=FALSE]
    colnames(tf)=colnames(heat@padj)
    tf=.asDf(tf)
    return(tf)
  }
  heat@padj=heat@padj[rid,,drop=FALSE]
  heat@value=heat@value[rid,,drop=FALSE]
  return(heat)
}

#' Plot heatmap.
#'
#' Given a data frame or a matrix, plot a heatmap with R package ComplexHeatmap. This functions is useful for plot results from movStat().
#'
#' plotHeatmap will automatelly determine whether to use discrete colors for a heatMtx.
#' If the number of unique elements in heatMtx (e.g., tf01 matrix from movStat()) is less than 9, then will use discrete colors.
#'
#' @param heatMtx a matrix of PAindex, which is normally from movAPAindex() or movPAindex().
#' @param refs reference condition, can beNULL/0350... Use to order other conditions with refs.
#' If refs is not NULL, then order other columns with refs, and put refs in the left most columns.
#' If refs!=NULL, cluster_cols=FALSE.
#' @param plotPre Output pdf file name prefix. If not NULL, then output plot to <plotPre>.pdf.
#' @param show_rownames show row names or not.
#' @param show_colnames show column names or not.
#' @param treeheight_row show row tree or not (=0 not show).
#' @param treeheight_col show column tree or not (=0 not show).
#' @param show_grid show grey grid or not.
#' @param main specify the title of the plot.
#' @param cluster_rows See ComplexHeatmap.
#' @param cluster_cols See ComplexHeatmap.
#' @param ... Dot arguments passed to ComplexHeatmap.
#' @return NULL
#' @examples
#' data(PACds)
#' ## get some DE results.
#' swDE=movUTRtrend(PACds, group='group', method='DE',
#' avgPACtag=10, avgGeneTag=20,
#' aMovDEPACRes=DEPAC, DEPAC.padjThd=0.01,
#' mindist=50, fisherThd=0.01, logFCThd=1, selectOne='fisherPV')
#' ## get heatmapResults.
#' heat=movRes2heatmapResults(swDE)
#' ## Filter switching genes.
#' heatUp=subsetHeatmap(heat, padjThd=0.05, valueThd=1)
## Plot heatmap using anther.embryo as reference column and sort other columns.
## From the heatmap ,we can see gene Os06g0682633 is shorter from anther to embryo (value=-8) and longer from embryo to maturePollen (value=7).
#' plotHeatmap(heatUp@value, refs=c('anther.embryo'), show_rownames=TRUE, plotPre=NULL)
#' plotHeatmap(heatUp@value, show_rownames=TRUE, plotPre=NULL, cluster_rows=TRUE)
#' @name plotHeatmap
#' @family visualization functions
#' @export
plotHeatmap <- function(heatMtx, refs=NULL, plotPre=NULL,
                        cluster_rows=TRUE, cluster_cols=FALSE,
                        show_rownames=FALSE, show_colnames=TRUE,
                        treeheight_row=0, treeheight_col=0, show_grid=FALSE, main="", ...) {

  # if (!require('ComplexHeatmap', quietly = TRUE, warn.conflicts = FALSE)) {
  #   stop("R package ComplexHeatmap is not installed!")
  # }

  heatMtx=as.matrix(heatMtx)
  cols=colnames(heatMtx)
  if (!is.null(refs)) {
    cluster_cols=FALSE
    if ( sum(refs %in% colnames(heatMtx))!=length(refs)) stop('plotHeatmap: refs not in heatMtx')
    cols=c(refs,cols[!(cols %in% refs)])
    heatMtx=heatMtx[,cols]
  }

  filename=NA
  if (!is.null(plotPre)) {
    filename=paste0(plotPre,'.pdf')
    pdf(filename)
  }

  show_row_dend=ifelse(treeheight_row==0, FALSE, TRUE)
  show_column_dend=ifelse(treeheight_col==0, FALSE, TRUE)

  #discrete
  colors=NULL
  if (max(heatMtx)<10) {
    x=unique(c(unlist(heatMtx)))
    if (length(x)<10) {
      x=x[order(x)]
      colors= RColorBrewer::brewer.pal(9, "YlOrRd")
      if (min(x)==0) {
        colors=colors[-1]
        colors=c('white', colors)
      }
      colors=colors[1:length(x)]
      names(colors)=x
    }
  }

  rect_gp =  grid::gpar(col = NA)
  if (show_grid) {
    rect_gp=grid::gpar(col = "grey80", lty = 1, lwd = 0.5)
  }

  if (!is.null(refs)) {
    if (sum(refs %in% cols)!=length(refs)) stop ('refs not all in heat')

    p=ComplexHeatmap::Heatmap(heatMtx[, refs, drop=F], cluster_rows = TRUE, cluster_columns=FALSE)

    rows=row_order(p)[[1]]
    heatMtx=heatMtx[rows,]

    ##ComplexHeatmap::
    if (is.null(colors)) {
      h=ComplexHeatmap::Heatmap(heatMtx,
                cluster_rows=FALSE, cluster_columns=FALSE,
                show_row_names=show_rownames, show_column_names=show_colnames,
                show_row_dend=show_row_dend, show_column_dend=show_column_dend, column_title=main,
                rect_gp=rect_gp, name="Value", ...)
    } else {
      h=ComplexHeatmap::Heatmap(heatMtx,
                cluster_rows=FALSE, cluster_columns=FALSE,
                show_row_names=show_rownames, show_column_names=show_colnames,
                show_row_dend=show_row_dend, show_column_dend=show_column_dend, column_title=main,
                col=colors, rect_gp=rect_gp, name="Value", ...)
    }
  } else {
    if (is.null(colors)) {
      h=ComplexHeatmap::Heatmap(heatMtx,
                cluster_rows=cluster_rows, cluster_columns=cluster_cols,
                show_row_names=show_rownames, show_column_names=show_colnames,
                show_row_dend=show_row_dend, show_column_dend=show_column_dend, column_title=main,
                rect_gp=rect_gp, name="Value", ...)
    } else {
      h=ComplexHeatmap::Heatmap(heatMtx,
                cluster_rows=cluster_rows, cluster_columns=cluster_cols,
                show_row_names=show_rownames, show_column_names=show_colnames,
                show_row_dend=show_row_dend, show_column_dend=show_column_dend, column_title=main,
                col=colors, rect_gp=rect_gp, name="Value", ...)
    }
  }
  print(h)
  if (!is.null(plotPre)) x=dev.off()
}



#' Output movStat results
#'
#' outputHeatStat outputs statistical results from movStat() to a text file and plots figures to a pdf file.
#'
#' This function will plot bar plots for heatStats$nsig, heatmaps for tf01, de01, deNum.
#' and also write statistical results to a text file.
#' @param heatStats a list from movStat().
#' @param ostatfile an output txt file name.
#' @param plotPre a file name prefix for output plots. Final output file names like <plotPre>.plots.pdf.
#' @param show_rownames whether show rownames in the heatmap, default is FALSE.
#' @return NULL
#' @examples
#' data(PACds)
#' swDE=movUTRtrend(PACds, group='group', method='DE',
#' avgPACtag=10, avgGeneTag=20,
#' aMovDEPACRes=DEPAC, DEPAC.padjThd=0.01,
#' mindist=50, fisherThd=0.01, logFCThd=1, selectOne='fisherPV')
#' ## First call the differential analysis and then call movStat to stat the results.
#' stat=movStat(object=swDE, padjThd=0.01, valueThd=1)
#' ## Output stat results into files.
#' ## The pdf file stores the plots about the number of significant DE PACs
#' ## and the overlap among different condition pairs.
#' outputHeatStat(heatStats=stat, ostatfile='3UTR_switching_DE.stat',
#'               plotPre='3UTR_switching_DE', show_rownames = TRUE)
#' @name outputHeatStat
#' @family visualization functions
#' @export
outputHeatStat <- function (heatStats, ostatfile, plotPre, show_rownames=FALSE) {
  #library(ggplot2, verbose = FALSE)

  plotfile=paste0(plotPre,'.plots.pdf')
  pdf(plotfile)
  heatStats$nsig$cond=rownames(heatStats$nsig)
  pd=reshape2::melt(heatStats$nsig, id='cond', variable.name = 'regulation', value.name = 'Number')
  #p=ggplot(pd, aes(x=cond, y=Number))  + ggplot2::geom_bar(stat="identity") + facet_grid(~regulation)
  p=ggplot(pd, aes(x=cond, y=Number))  + ggplot2::geom_bar(stat="identity") + stat_plot_theme
  print(p)

  #heatmaps
  if (!is.null(heatStats$tf01)) {
    plotHeatmap(heatStats$tf01,
                cluster_rows=TRUE, cluster_cols=FALSE,
                show_rownames=show_rownames, show_colnames=TRUE,
                treeheight_row=0, treeheight_col=0, main='Sig. or not between conditions')
  }

  if (!is.null(heatStats$de01)) {
    plotHeatmap(heatStats$de01,
                cluster_rows=TRUE, cluster_cols=FALSE,
                show_rownames=show_rownames, show_colnames=TRUE,
                treeheight_row=0, treeheight_col=0, main='Sig. or not in a condition')
  }


  if (!is.null(heatStats$deNum)) {
    plotHeatmap(heatStats$deNum,
                cluster_rows=TRUE, cluster_cols=FALSE,
                show_rownames=show_rownames, show_colnames=TRUE,
                treeheight_row=0, treeheight_col=0, main='Sig. times in a condition')
  }
  dev.off()
  cat('>>> [nsig] to ',plotfile,'\n')

  #novp: venn plot overlapping DE number among conditions
  # if (!is.null(heatStats$ovp)) {
  #   #library(VennDiagram, verbose = FALSE)
  #   plotfile=paste0(plotPre,'.ovpVenn.emf')
  #   n=nrow(heatStats$ovp)
  #   if (nrow(heatStats$ovp)>=5) n=5
  #   col = brewer.pal(n, "Set1")
  #   x=venn.diagram(heatStats$siglist[1:n], fill=col,cex=2, cat.fontface=4,  filename=plotfile)
  #   cat('>>> [ovp] to ',plotfile,'\n')
  # }

  LOG('Number of sig. items',ostatfile, append=FALSE, star=TRUE)
  LOG(heatStats$nsig, ostatfile)
  if (!is.null(heatStats$ovp)) {
    LOG('Number of overlapping items',ostatfile, append=TRUE, star=TRUE)
    LOG(heatStats$ovp, ostatfile)
  }
  cat('>>> [txt] to ',ostatfile,'\n')
}


# -------------- *** ATCG **** ---------------

#' Remove internal priming.
#'
#' removePACdsIP removes internal priming artificats from a PACdataset object.
#'
#' If there are continuous [conA] As or at least >=[sepA] As or at least one [ipGrams] in the [up, dn] window of a PAC, then it is considered as IP.

#' @param PACds a PACdataset.
#' @param bsgenome a BSGenome/FaFile to get chromosome sequences.
#' The PACds@anno$chr should all be in bsgenome$seqnames.
#' @param returnBoth If TRUE, then return a list(real, ip), otherwise return only the real PACds.
#' @param up the start position of the window to scan, see faFromPACds().
#' @param dn the end position of the window to scan, see faFromPACds().
#' If up=-10 and dn=10 then will scan the [-10, 0, +10] window (21 nt).
#' If up=-10 and dn=-1 then will scan the [-10, -1] window (10 nt, not including the PA pos).
#' @param conA the number of continuous As in a window to define IP. Default is NA.
#' @param sepA the number of total As in a window to define IP. Default is NA, which means not considering total As.
#' @param ipGrams the grams to define IP. Default is NA. If any one gram in ipGrams is located within [up, dn] region then it is IP.
#' @param chrCheck if TRUE, then all chr in PACds should be in bsgenome, otherwise will ignore those non-consistent chr rows in PACds.
#' @return Dependent on returnBoth.
#' @examples
#' data(PACds)
#' library("BSgenome.Oryza.ENSEMBL.IRGSP1")
#' bsgenome <- BSgenome.Oryza.ENSEMBL.IRGSP1
#' ## remove internal priming within -140 ~ +10 of PA, without considering total A#
#' removePACdsIP(PACds, bsgenome, returnBoth=TRUE, up=-140, dn=10, conA=6, sepA=NA)
#' ## remove internal priming within -10 ~ +10 of PA, if there is any AAAAAA or total 7 As.
#' removePACdsIP(PACds, bsgenome, returnBoth=TRUE, up=-10, dn=10, conA=6, sepA=7)
#' ## remove internal priming within -10 ~ +10 of PA, if there is any one ipGrams.
#' removePACdsIP(PACds, bsgenome, returnBoth=TRUE, up=-5, dn=10,
#'              ipGrams=c('AAAA', 'AGAA', 'AAGA', 'AAAG'))
#' @name removePACdsIP
#' @family PACdataset functions
#' @export
removePACdsIP <- function(PACds, bsgenome, returnBoth=TRUE, up=-10, dn=10, conA=NA, sepA=NA, ipGrams=NA, chrCheck=TRUE) {

  PACds@counts=asAnyMatrix(PACds@counts)
  if (!isChrConsistent(pacds=PACds, obj=bsgenome, allin=TRUE)) {
    cat('PACds chr not all in seqnames of bsgenome\n')
    if (chrCheck) {
      stop('Please check chr names!\n')
    } else {
      cat("chrCheck=FALSE, Remove rows in PACds with chr name not in bsgenome\n")
      PACds=subsetPACds(PACds, chrs=chrnames(bsgenome), verbose=TRUE)
    }
  }

  seq=faFromPACds(PACds, bsgenome, what='updn', fapre=NULL, byGrp=NULL, up=up, dn=dn)

  idx1=c()
  if (!is.na(conA)) {
    As=paste0(rep('A',conA), collapse = '')
    vAs=Biostrings::vcountPattern(As, seq)
    idx1=which(vAs!=0)
  }

  idx2=c()
  if (!is.na(ipGrams[1])) {
    mydict <- Biostrings::DNAStringSet(ipGrams)
    names(mydict)<- paste0("gram", 1:length(ipGrams), sep="")
    mypdict <- Biostrings::PDict(mydict)
    m <- Biostrings::vwhichPDict(mypdict, seq)
    vG=Biobase::listLen(m)
    idx2=which(vG!=0)
  }

  idx3=c()
  if (!is.na(sepA)) {
    vA=Biostrings::vcountPattern('A', seq)
    idx3=which(vA>=sepA)
  }

  idx=unique(c(idx1, idx2, idx3))

  cat(length(idx),'IP PACs;',length(seq)-length(idx),'real PACs\n')

  if (length(idx)==0) {
    if (returnBoth) return(list(real=PACds, ip=new("PACdataset")))
    return(PACds)
  }

  IP=PACds
  IP@counts=IP@counts[idx, ,drop=F]
  IP@anno=IP@anno[idx, ,drop=F]
  PACds@counts=PACds@counts[-idx, ,drop=F]
  PACds@anno=PACds@anno[-idx, ,drop=F]
  if (returnBoth) return(list(real=PACds, ip=IP))
  return(PACds)
}


# Get seq from a genomicRanges object. Get seq from start-end, and if strand=-, then reverse and complement the seq.
# @param gr GenomicRanges, if no names, then output seq title is chr;strand;start;end.
# @param bsgenme BSgenome
# @param fafile Output fasta file name.
# @return Dependent on fafile, if fafile!=NULL, then return a DNAStringSet object.
getSeqFromGR<-function(gr, bsgenome, fafile=NULL, width=400) {


  if (!(class(bsgenome) %in% c('BSgenome','FaFile'))) stop("bsgenome should be BSgenome/FaFile\n")

  seq <-Biostrings::getSeq(bsgenome, gr)

  if (is.null(names(gr))) {
    gr=.asDf(gr)
    names(seq)=paste0(gr$seqnames,';',gr$strand,';',gr$start,';',gr$end,';',gr$end-gr$start+1)
  }
  if (!is.null(fafile)) {
    Biostrings::writeXStringSet(seq, fafile, format='fasta', width=width)
    return()
  }
  return(seq)
}

## get chr names from bsgenome or fafile
chrnames<-function(object) {
  if (inherits(object, 'BSgenome')) {
    objchr=GenomeInfoDb::seqnames(object) # 2023/5/24 fixed bugs
  } else if (inherits(object, "FaFile")) {
    objchr=names(GenomeInfoDb::seqlengths(object))
  } else {
    stop(cat("Don't know how to get chr names from object,",class(object),'\n'))
  }
  return(objchr)
}

#' Extract sequences from a PACdataset.
#'
#' faFromPACds extracts many kinds of sequences from a PACdataset.
#'
#' This function can export the sequences surrounding PACs, the sequences of genomic regions the PACs located,
#' and the gene sequences. If export sequences of regions/genes/etc., only one region seuqence will be exported.
#' @param PACds a PACdataset.
#' @param bsgenome BSgenome or FaFile object storing chromosome seqs, or a fasta file name
#' @param what the value can be updn, pac, region, gene.
#' \itemize{
#' \item{updn}: use -up and -dn to define the upstream and downstream region around PAC's coord.
#' \item{pac}: output the range of PACs by PACds@anno$UPA_start~UPA_end.
#' \item{region}: output the genomic region sequence of PACs by PACds@ano$ftr_start, ftr_end.
#' \item{gene}: output the gene sequence by PACds@gene_start, gene_end.
#' }
#' @param byGrp to separately output sequences to different fa files.
#' The value can be NULL / ftr / c('ftr','strand') / list(ftr=c('3UTR','5UTR'), strand=c('+'),'-').
#' @param up paramter for what=updn, specifying the upstream region from the PAC.
#' @param dn paramter for what=updn, specifying the downstream region from the PAC.
#' PAC is the 0 position.
#' e.g., up=-300, dn=100 to subset 401 nt (PAC is the 0 or 301 position, upstream 1..300 [or -300..-1], PA301 [or 0], downstream 302..401 [or 1..100])
#' e.g., up=0, dn=0, will output the nucleotide at the PAC position.
#' @param fapre a prefix for output file. If fapre=NULL, then return stringSet, but this is only valid when byGrp=NULL.
#' @param chrCheck if TRUE, then all chr in PACds should be in bsgenome, otherwise will ignore those non-consistent chr rows in PACds.
#' @return File names or a stringSet. If up=-300, dn=100, then the output sequence is 401nt and PAC position is the 301st.
#' @examples
#' library("BSgenome.Athaliana.TAIR.TAIR9")
#' bsgenome <-Athaliana
#' pacds=makeExamplePACds()
#' ## Get sequences of PAC ranges.
#' faFromPACds(pacds, bsgenome, what='pac', fapre='pac')
#'
#' ## bsgenome is a fasta file
#' fapath <- 'Arab_TAIR9_chr_all.fas'
#' faFromPACds(pacds, bsgenome=fapath, what='updn', fapre=NULL,
#'             byGrp=NULL, up=-300, dn=100, chrCheck=TRUE)
#'
#' ## Get upstream 300nt and downstream 100nt sequences around PACds.
#' faFromPACds(pacds, bsgenome, what='updn', fapre='updn', up=-300, dn=100)

#' ## Get PAC sequences and output by different genomic regions, e.g.,
#' ## 3UTR PACs, intron PACs...
#' faFromPACds(pacds, bsgenome, what='updn', fapre='updn',
#'             up=-300, dn=100, byGrp='ftr')
#' faFromPACds(pacds, bsgenome, what='updn', fapre='updn',
#'             up=-300, dn=100,
#'             byGrp=list(ftr='3UTR'))
#' faFromPACds(pacds, bsgenome, what='updn', fapre='updn', up=-300, dn=100,
#'             byGrp=list(ftr='3UTR', strand=c('+','-')))

#' ## Get sequences for genes with PACs.
#' faFromPACds(pacds, bsgenome, what='region', fapre='region', byGrp='ftr')
#' @name faFromPACds
#' @family APA signal functions
#' @export
faFromPACds <- function(PACds, bsgenome, what='updn', fapre=NULL, byGrp=NULL, up=-300, dn=100, chrCheck=TRUE) {

  if (inherits(bsgenome, what='character')) {
    bsgenome=Rsamtools::FaFile(bsgenome)
    if (!file.exists(bsgenome$index)) {
      cat("bsgenome is a fasta file, indexing...\n")
      Rsamtools::indexFa(bsgenome$path)
    }
  }

  if (!isChrConsistent(pacds=PACds, obj=bsgenome, allin=TRUE)) {
    cat('PACds chr not all in seqnames of bsgenome\n')
    if (chrCheck) {
      stop('Please check chr names!\n')
    } else {
      cat("chrCheck=FALSE, Remove rows in PACds with chr name not in bsgenome\n")
      PACds=subsetPACds(PACds, chrs=chrnames(bsgenome), verbose=TRUE)
    }
  }


  anno=PACds@anno
  fapre1=fapre
  if (is.null(fapre)) fapre=what

  if (what=='updn') {
    if (!AinB(c('coord','chr','strand'), colnames(anno))) stop("what=updn, no chr/strand/coord in PACds@anno")
    anno$start=anno$coord+up
    anno$start[anno$strand=='-']=anno$coord[anno$strand=='-']-dn
    anno$end=anno$start+(dn-up)
    anno$title=paste0(rownames(anno),':',anno$chr,';',anno$strand,';',anno$coord)
  } else if (what=='pac') {
    if (!AinB(c('UPA_start','UPA_end','chr','strand'), colnames(anno))) stop("what=pac, no chr/strand/UPA_start/UPA_end in PACds@anno")
    anno$start=anno$UPA_start
    anno$end=anno$UPA_end
    anno$title=paste0(rownames(anno),':',anno$chr,';',anno$strand,';',anno$UPA_start,';',anno$UPA_end)
  } else if (what=='region') {
    if (!AinB(c('ftr_start','ftr_end','chr','strand'), colnames(anno))) stop("what=region, no chr/strand/ftr_start/ftr_end in PACds@anno")
    anno$start=anno$ftr_start
    anno$end=anno$ftr_end
    anno$title=paste0(anno$ftr,':',anno$chr,';',anno$strand,';',anno$ftr_start,';',anno$ftr_end)
  } else if (what=='gene') {
    if (!AinB(c('gene_start','gene_end','chr','strand'), colnames(anno))) stop("what=gene, no chr/strand/gene_start/gene_end in PACds@anno")
    anno$start=anno$gene_start
    anno$end=anno$gene_end
    anno$title=paste0(anno$gene,':',anno$chr,';',anno$strand,';',anno$gene_start,';',anno$gene_end)
  }
  #anno=anno[,c('chr','strand','start','end')]
  annoList=list()
  fileList=c()

  #no grup
  if (is.null(byGrp)) {
    fafile=paste0(fapre,'.fa')
    annoList=list(anno)
    fileList=c(fafile)
  }

  #If it is c('ftr','strand'), first convert to list
  if (is.character(byGrp)) {
    if (length(byGrp)>2) stop("at most two groups in byGrp")
    bys=byGrp
    byGrp=list()
    for (i in bys) {
      vals=unique(anno[,i])
      if (length(vals)>0) byGrp[[i]]=vals
    }
  }

  #read each group to annoList
  if (is.list(byGrp)) {
    if (!AinB(names(byGrp), colnames(anno))) stop("names(byGrp) not in PACds@anno")
    if (length(byGrp)>2) stop("at most two groups in byGrp")
    for (i in 1:length(byGrp[[1]])) {
      bywhat=names(byGrp)[1]
      byval=byGrp[[1]][i]
      if (length(byGrp)==1) {
        anno1=anno[anno[,bywhat]==byval, ]
        annoList[[length(annoList)+1]]=anno1
        fileList=c(fileList, paste0(fapre,'.',byval,'.fa'))
      } else {
        for (j in 1:length(byGrp[[2]])) {
          bywhat2=names(byGrp)[2]
          byval2=byGrp[[2]][j]
          anno1=anno[anno[,bywhat]==byval & anno[,bywhat2]==byval2, ]
          annoList[[length(annoList)+1]]=anno1
          fileList=c(fileList, paste0(fapre,'.',byval,byval2,'.fa'))
        }
      }
    }
  }

  if (length(fileList)>1) {
    if (is.null(fapre1)) {
      stop("more than one fa output, fapre should not be NULL\n")
    }
  }

  for (i in 1:length(annoList)) {
    anno=annoList[[i]]
    fafile=fileList[i]
    #remove dup regions
    anno=unique(anno[,c('title','start','end','chr','strand')])
    #remove end<start lines
    anno=anno[anno$end>=anno$start, ]
    anno=anno[anno$start>0 & anno$end>0, ]

    #If reach out a chr, then delete
    #otherwise, subseq will raise error: some ranges are out of bounds
    chrlen=GenomeInfoDb::seqlengths(bsgenome)
    n1=nrow(anno)
    for (chr in unique(anno$chr)) {
      idx=which(anno$chr==chr & anno$end>chrlen[chr])
      if (length(idx)>0) anno=anno[-idx, ]
    }
    n2=nrow(anno)

    gr <- GRanges(seqnames =anno$chr,
                  ranges = IRanges::IRanges(start=anno$start, end=anno$end),
                  strand = anno$strand)
    names(gr)=anno$title

    if (is.null(fapre1)) return(getSeqFromGR(gr, bsgenome, fafile=NULL))

    x=getSeqFromGR(gr, bsgenome, fafile=fafile)
    cat(nrow(anno),'>>>',fafile,'\n')
  }

  return(fileList)
}


#' Plot single nucleotide profile
#'
#' plotATCGforFAfile plots single nucleotide profile for given fasta files.
#'
#' @param faFiles fasta file names, e.g., c(xx.fa, yy.fa)
#' @param ofreq FALSE/TRUE. If true, output xx.freq, yy.freq, [pos, A C G T]
#' @param opdf TRUE/FALSE. If true, output xx.pdf, yy.pdf
#' @param refPos the position of the poly(A) site. Value can be NULL or a number.
#' E.g., refPos=300, then the relative position to the PAC (pos=0) is pos-301.
#' @param start to narrow the plot region from start~end
#' @param end to narrow the plot region from start~end
#' @param filepre prefix of output file, if mergePlots=T, then output file=<filepre>.pdf; otherwise, each file is <file>.<filepre>.pdf
#' @param mergePlots if TRUE, then combine all plots with facet into one plot. If opdf=TRUE, then filepre should be provided.
#' In the merged PDF file, the title of each plot is like <fafile> #seqNumber. Long file name will be automatically shortened.
#' @return NULL. Plot figures to PDF file or screen.
#' @examples
#' data(PACds)
#' library("BSgenome.Oryza.ENSEMBL.IRGSP1")
#' bsgenome <- BSgenome.Oryza.ENSEMBL.IRGSP1
#' faFiles=faFromPACds(PACds, bsgenome, what='updn', fapre='updn', up=-100, dn=100, byGrp='ftr')
#' faFiles=c("updn.3UTR.fa", "updn.Ext_3UTR.fa", "updn.intergenic.fa", "updn.intron.fa" )
#' ## plot single nucleotide profile for a fa file
#' plotATCGforFAfile (faFiles="updn.3UTR.fa", ofreq=TRUE, opdf=TRUE, refPos=301)
#' ## plot multiple fa files
#' fafiles=faFromPACds(PACds, bsgenome, what='updn', fapre='400nt',
#'        up=-300, dn=100, byGrp='ftr', chrCheck=FALSE)
#' plotATCGforFAfile (faFiles=fafiles, ofreq=TRUE, opdf=TRUE, refPos=301)
#' ## plot multiple fa files into one PDF
#' fafiles=faFromPACds(PACds, bsgenome, what='updn', fapre='400nt',
#'        up=-300, dn=100, byGrp='ftr', chrCheck=FALSE)
#' plotATCGforFAfile (faFiles=fafiles, ofreq=TRUE, opdf=TRUE,
#'       refPos=301, mergePlots=TRUE, filepre='allplots')
#' @name plotATCGforFAfile
#' @family APA signal functions
#' @export
plotATCGforFAfile <- function (faFiles, ofreq=FALSE, opdf=TRUE, refPos=NULL, start=NA, end=NA, filepre='', mergePlots=FALSE) {

  #library(Biostrings, verbose = FALSE)
  #library(ggplot2)

  if (filepre=='' & (opdf & mergePlots)) stop("mergePlots & opdf=TRUE, should provide output file name: filepre\n")

  dmerge=c()

  if (opdf) {
    if (mergePlots) {
      pdf(file=paste0(filepre,'.pdf'))
    }
  }


  for (fafile in faFiles) {
    seq=Biostrings::readDNAStringSet(fafile, format="fasta")
    nseq=length(seq)
    if (!is.na(start) | !is.na(end)) {
      seq=XVector::subseq(seq, start=start, end=end)
      refPos=NULL
    }
    din<-Biostrings::consensusMatrix(seq, as.prob=TRUE, baseOnly=TRUE)
    din<- .asDf(t(din))
    din$other <- NULL
    din$pos <- c(1:nrow(din))
    if (!is.null(refPos)) {
      din$pos=din$pos-refPos
    }
    opre=gsub('\\.fa$|\\.fasta$','',fafile, ignore.case = TRUE)
    omain=paste0(opre, ' #',nseq)
    if (filepre!='') {
      opre=paste0(opre,'.',filepre)
    }
    if (ofreq) {
      din=din[, c('pos','A','C','G','T')]
      ofreqfile=paste0(opre,'.freq')
      write.table(din, file=ofreqfile, col.names = TRUE, row.names = FALSE, sep="\t", quote=F)
      cat('>>>',ofreqfile,'\n')
    }

    if (mergePlots) {
      if (length(dmerge)==0) {
        dmerge=cbind(source=omain, din)
      } else {
        dmerge=rbind(dmerge, cbind(source=omain,din))
      }
    }
    din <- reshape2::melt(din,id.vars='pos',variable.name = 'base', value.name = 'freq')
    ATCGpicture<-ggplot(din,aes(x=din$pos,y=din$freq,colour=din$base)) + geom_line() +
      xlab("Position") + ylab('Base component') + theme(legend.title=element_blank()) + theme_bw() +
      ggtitle(omain)

    if (opdf & !mergePlots) {
      pdf(file=paste0(opre,'.pdf'))
      print(ATCGpicture)
      dev.off()
      cat('>>>',paste0(opre,'.pdf'),'\n')
    }
    ## only plot the last one
    if (!opdf & !mergePlots & faFiles[length(faFiles)]==fafile) print(ATCGpicture)
  }

  if (mergePlots) {
    din <- reshape2::melt(dmerge, id.vars=c('source','pos'), variable.name = 'base', value.name = 'freq')
    din$source=as.character(din$source)
    #if file name is too long, remove common chars
    #if no common chars, change the first 1-30 chars to ...
    l=min(nchar(din$source))
    if (l>=30) {
      din$source =.autoDetectCommonString(unique(din$source), sbj=din$source, beAll=TRUE)
      l=min(nchar(din$source))
      if (l>=30) {
        din$source=paste0('...', substr(din$source, 31, nchar(din$source)))
      }
    }

    ATCGpicture <- ggplot(din, aes(x=pos, y=freq, group=base)) +
          geom_line(aes(colour=base)) +
          xlab("Position") + ylab('Base component') + theme(legend.title=element_blank()) + theme_bw()
          facet_grid(. ~ source)
    ATCGpicture = ATCGpicture + facet_wrap(~ source, ncol=2)

    print(ATCGpicture)
    if (opdf) {
      dev.off()
      cat('>>>',paste0(filepre,'.pdf'),'\n')
    }
  }

}


#  genKgrams($k,$withvalue):@kgrams
#  useage: @kgrams=genKgrams(6,0);
#  withvalue: =0, like AATAAA=0
genKgrams <- function(k, withvalue=NULL) {
  s='';
  cnt=4^k;
  kgram=c()
  for (i in 0:(cnt-1)) {
    num=i;
    for (j in 1:k) {
      res=num%%4;
      num=as.integer(num/4);
      s=paste0(res, s)
    }
    kgram=c(kgram, s)
    s=''
  }
  kgram=unlist(lapply(kgram, function (x) {return ( chartr(old='0123',new='ATCG',x))}))
  if (!is.null(withvalue)) kgram=paste0(kgram,'=',withvalue)
  return(kgram)
}

#' Count kgrams
#'
#' kcount counts given kgrams from given sequences. Kgrams can be specified by character vector or by K=N.
#'
#' @param seq a DNAStringSet object.
#' @param fafile fasta file name(s). fafile or seq can only set one.
#' If fafile=c(A='fA1.fa',B='fB2.fa'), then the output depends on mergeFile=TRUE/FALSE. And the column names are set as A for count, A_perc for percent.
#' If names(fafile) is not provided, then output columns will be A1_perc, B2_perc.
#' @param mergeFile if TRUE, then when multiple file names are given by fafile, files are combined first.
#' @param from from and to defines the range of pA to search kgrams.
#' @param to from and to defines the range of pA to search kgrams.
#' @param grams specify the grams to count. Value can be c('aataa',...) or 'V1' (means AATAAA and its 1nt variants) or mouse/mm/mm10 (means PAS of mouse).
#' @param k specify k grams. E.g., k=6 means all hexamers.
#' @param kfile specify a kgram file, with each line being a kgram. The file does not contain header line.
#'        Only one parameter of grams/k/kfile can be provided.
#' @param sort TRUE/FALSE, whether to sort kgrams by their counts.
#' @param topn whether to output topn grams.
#' @param perc whether also to output percentage.
#' @return A data frame with columns [grams, count, <perc>] or [grams, file1,..., fileN, ..., file1_perc, ..., fileN_perc]
#' @examples
#' kcount(seq=seq, grams=c('AATAAA','ATTAAA'))
#' kcount(seq=seq, k=3)
#' kcount(seq=seq, k=6, from=260, to=290)
#' kcount(fafile='updn.fa', k=6, from=260, to=290)
#' kcount(fafile='updn.fa', grams='v1', from=260, to=290)
#' kcount(fafile='updn.fa', grams='mouse', from=240, to=310)
#' kcount(fafile=c('updn.fa', 'updn2.fa'), grams='mouse', from=240, to=310, mergeFile=FALSE)
#' kcount(fafile=c(A='updn.fa', B='updn2.fa'), grams='mouse', from=240, to=310, mergeFile=FALSE)
#' kcount(fafile=c('updn.fa', 'updn2.fa'), grams='mouse', from=240, to=310, mergeFile=TRUE)
#' @name kcount
#' @family APA signal functions
#' @export
kcount <- function (seq=NULL, fafile=NULL, grams=NULL, k=NULL, kfile=NULL, mergeFile=FALSE, from=NA, to=NA, sort=TRUE, topn=50, perc=TRUE) {
  #library(Biostrings, verbose = FALSE)

  if (is.null(seq) & is.null(fafile)) stop("kcount: seq and fafile =NULL")
  if ( (!is.null(grams) & !is.null(k)) | (!is.null(grams) & !is.null(kfile)) | (!is.null(k) & !is.null(kfile)) ) stop("kcount: grams or k or kfile, only one option")

  grams=getVarGrams(grams)

  if (!is.null(k)) grams=genKgrams(k)
  if (!is.null(kfile)) grams=read.table(kfile, header=FALSE)$V1
  mydict <- Biostrings::DNAStringSet(grams)
  names(mydict)<- paste0("gram", 1:length(grams), sep="")
  mypdict <- Biostrings::PDict(mydict)

  seqs=list()
  if (!is.null(fafile)) {
    if (mergeFile) {
      seqs=list(Biostrings::readDNAStringSet(fafile, format="fasta"))
    } else {
      for (fa in fafile) {
        seqs=c(seqs, list(Biostrings::readDNAStringSet(fa, format="fasta")))
      }
      if (!is.null(names(fafile))) {
        names(seqs)=names(fafile)
      } else {
        names(seqs)=.autoDetectCommonString(fafile, fafile)
      }
    }
  } else {
    seqs=list(seq)
  }
  if (length(seqs)==1) names(seqs)='count'

  cnts=c()
  for (i in 1:length(seqs)) {

    seq=seqs[[i]]

    if (!is.na(from) | !is.na(to)) {
      mlen=min(width(seq))
      end=max(c(from,to)[!is.na(c(from,to))])
      if (mlen<end) {
        id=which(width(seq)<end)
        seq=seq[-id]
        cat('to > seq len, remove',length(id), 'short seqs\n')
      }
      seq=XVector::subseq(seq, start=from, end=to)
    }


    cnt=Biostrings::vcountPDict(mypdict, seq,
                    max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
                    algorithm="auto", collapse=1)
    cnt=.asDf(cbind(grams, count=cnt)); cnt$count=as.numeric(cnt$count)
    if (i==1) {
      cnts=cnt
    } else {
      cnts=cbind(cnts, cnt[,2])
    }
  }
  colnames(cnts)=c('grams',names(seqs))

  if (perc) {
    perc=lapply(cnts[, -1, drop=F], function(par) par/sum(par))
    if (length(perc)==1) {
      names(perc)='perc'
    } else {
      names(perc)=paste0(colnames(cnts)[-1],'_perc')
    }
    cnts=cbind(cnts, perc)
  }

  if (sort) cnts=cnts[order(cnts[,2], decreasing = TRUE),]
  if (!is.null(topn)) cnts=cnts[1:min(topn, nrow(cnts)), ,drop=F]

  return(cnts)
}

#' Plot seqlogo for poly(A) signals.
#'
#' plotSeqLogo plots seqlogo for poly(A) signals by utilizing the motifStack package.
#'
#' @param pas a vector of k-grams, normally from annotateByPAS().
#' @return A plot of seqlogo.
#' @examples
#' data(PACds)
#' library("BSgenome.Oryza.ENSEMBL.IRGSP1")
#' bsgenome <- BSgenome.Oryza.ENSEMBL.IRGSP1
#' PACdsPAS=annotateByPAS(PACds, bsgenome, grams='V1',
#'                        from=-50, to=-1, label=NULL)
#' pas=PACdsPAS@anno$V1_gram[!is.na(PACdsPAS@anno$V1_gram)]
#' plotSeqLogo(pas)
#' @name plotSeqLogo
#' @family APA signal functions
#' @export
plotSeqLogo <- function(pas) {
  # if (!require('motifStack', quietly = TRUE, warn.conflicts = FALSE)) {
  #   stop("motifStack is not installed!")
  # }
  pas=Biostrings::DNAStringSet(pas)
  pfm=Biostrings::consensusMatrix(pas, baseOnly=TRUE)
  pfm=pfm[which(rownames(pfm) %in% c("A","C","G","T")),]
  rownames(pfm)=c("A","C","G","U")
  pfm.freq=apply(pfm,2, function(x) x/sum(x))
  motifStack::plotMotifLogo(pfm.freq, yaxis = TRUE, font="Arial",
                colset=c("#D00001","#2000C7","#FFB32C","#00811B"),
                xlab="",ylab="", ic.scale=TRUE)
}



# TODO: motif analyses tools
# http://www.bioconductor.org/packages/release/bioc/manuals/motifRG/man/motifRG.pdf
# http://manuals.bioinformatics.ucr.edu/home/ht-seq#TOC-String-Handling-Utilities-in-R-s-Base-Distribution
# http://bioconductor.org/packages/release/bioc/html/motifStack.html
# http://faculty.ucr.edu/~tgirke/HTML_Presentations/Manuals/Rsequences/Rsequences.pdf

# ---------------- *** mov Wrapper *** ------------------

#Get condpair df
# - obj: PACds / character like len_0700_vs_0140
# - group: sample group
# Return: df(rownames=cond1.cond2, [cond1,cond2])
# Ex.
# .getCondPairs(obj=c("len_0350_vs_0140","len_0700_vs_0140","len_1200_vs_0140"), group='len')
.getCondPairs <- function(obj, group) {

  if (inherits(obj, "PACdataset")) {
    if (!(group %in% colnames(obj@colData))) {
      stop(sprintf(".getCondPairs: group=%s not in colData(PACds)",group))
    }

    conds=as.character(unique(obj@colData[,group]))
    conds1=conds
    condPairs=t(combn(conds,2))
    colnames(condPairs)=c('cond1','cond2')
    rownames(condPairs)=paste0(condPairs[,1],'.',condPairs[,2])
  }

  #Fixed format of condpair name, normally from resultsNames(dds)
  # - rn, like len_0350_vs_0140
  if (inherits(obj, what='character')) {
    f1rn=unlist(strsplit(obj[grep(factor1,obj)],'_vs_'))
    cond2=f1rn[seq(2,length(f1rn),2)]
    cond1=gsub(paste0('^',factor1),'',f1rn[seq(1,length(f1rn),2)])
    condPairs=cbind(cond1,cond2)
    rownames(condPairs)=paste0(cond1,'.',cond2)
  }
 # for (i in 1:ncol(condPairs)) {
 #   condPairs[,i]=as.character(condPairs[,i])
 # }
  return(.asDf(condPairs))
}

#' Get DE genes
#'
#' movDEGene first summerizes all PACs in a gene to get gene expression level, and then calls DESeq2 for DE analysis.
#' It performs DE analyses for all condition pairs of the given group. If you only need a specific condition pair, call
#' subsetPACds() first to filter used conditions.
#'
#' @param PACds a PACdataset.
#' @param group a sample group name to get all conditions under the group.
#' @param method only DESeq2 is allowed for now
#' @param minSumPAT filter genes with total PAT number across all conditions >=minSumPAT
#' @return An object of movDEGeneRes, with slots[group, conds, pairwise=list[value=log2FC_<cond1.cond2>, padj=padj_<cond1.cond2>]]
#' @examples
#' data(PACds)
#' ## Detect DE genes for two conditions.
#' ## Subset two conditions first.
#' pacds=subsetPACds(PACds, group='group',cond1='anther', cond2='embryo')
#' ## Detect DE genes using DESeq2 method,
#' ## only genes with total read counts in all samples >=50 are used.
#' DEgene=movDEGene(PACds=pacds, method='DESeq2', group='group', minSumPAT=50)
#' ## Detect DE genes in all pairwise conditions.
#' DEgene=movDEGene(PACds=PACds, method='DESeq2', group='group', minSumPAT=50)
#' @name movDEGene
#' @family comparison functions
#' @export
movDEGene <- function (PACds, group, method='DESeq2', minSumPAT=10) {

  library(DESeq2, verbose = FALSE)

  method=toupper(method)
  if (!(method %in% c('DESEQ2'))) method='DESEQ2'

  dds=PACds2geneDs(PACds)

  #pre-filtering
  if (minSumPAT>1) {
    keep <- rowSums(counts(dds)) >= minSumPAT
    dds <- dds[keep,]
  }

  #all condpairs DE
  geneDE=DESeqNoInteractionAll(dds=dds, factor1=group, factor2=NULL, norm=FALSE)
  geneDERes=digDESeqNoInteractionAll(geneDE, factor1=group, factor2=NULL, ref1=NULL, ref2=NULL, pairwise=TRUE, series=FALSE, levels1=NULL, levels2=NULL)

  condPairs=.getCondPairs(PACds, group)
  geneDEHeat=getHeatmapRes(object=geneDERes, colData=condPairs, padjThd=2, valueThd=0)

  #condPairs=.getCondPairs(colnames(geneDERes@pairwise[[1]]), group)
  geneDERes=new("movDEGeneRes", group=group, method=method,  conds=condPairs, pairwise=list(padj=geneDEHeat@padj, value=geneDEHeat@value))
  return(geneDERes)
}


#' Get DE PACs
#'
#' movDEPAC calls DESeq2/DEXSeq or uses chisq.test for DE analysis.
#' It performs DE analyses for all condition pairs of the given group. If you only need a specific condition pair, call
#' subsetPACds first to filter used conditions.
#'
#' @param PACds a PACdataset.
#' @param group a sample group name to get all conditions under the group.
#' @param method can be DEXseq, DESeq2, chisq, ignoring cases.
#'               chisq is as in Shulman 2019 et al, which first performs chisq.test for all PACs in a gene to filter sig. genes and then performs chisq.test for each PAC.
#'               If there are multiple reps, then will pool reps first before chisq.
#' @param minSumPAT filter PACs with total PAT number across all conditions >=minSumPAT
#' @param chisqPadjust if T then output adjust pvalues instead of raw pvalues.
#' @return An object of movDEPACRes, including slots [group, conds,
#'                                  pairwise=list[value=log2FC_<cond1.cond2> for DESeq2/DEXSeq; 1-genePv for Chisq,
#'                                  padj=padj_<cond1.cond2>],
#'                                  PACusage(if method=DEXseq)]
#' @examples
#' data(PACds)
#' ## Detect DE PACs in all pairwise conditions using DESeq2 method,
#' ## only PACs with total read counts in all samples >=20 are used.
#' DEPAC=movDEPAC(PACds, method='DESeq2', group='group', minSumPAT=20)
#'
#' ## Detect DE PACs using DEXseq method (may be slow).
#' DEXPAC=movDEPAC(PACds, method='DEXseq', group='group', minSumPAT=20)
#'
#' ## Detect DE PACs using chisq-PACdsPAS method.
#' DEqPAC=movDEPAC(PACds, method='chisq', group='group', minSumPAT=20)

#' @name movDEPAC
#' @family comparison functions
#' @export
movDEPAC <- function (PACds, group, method='DESeq2', minSumPAT=10, chisqPadjust=TRUE) {
  #library(DESeq2, verbose = FALSE)
  #library(DEXSeq, verbose = FALSE)

  method=toupper(method)
  if (!(method %in% c('DESEQ2','DEXSEQ','CHISQ'))) stop("method must be DESEQ2/DEXSEQ/CHISQ")

  if (method=='DESEQ2') {
    dds=PACds2DESeqDds(PACds, noIntergenic = TRUE)
  } else if (method=='DEXSEQ') {
    dds=PACdataset2Dxd(PACds, group=group)
  }

  #pre-filtering
  if (method %in% c('DESEQ2','DEXSEQ') & minSumPAT>1) {
    keep <- rowSums(counts(dds)) >= minSumPAT
    dds <- dds[keep,]
  }

  if (method=='CHISQ') { # only for APA genes
    PACds=subsetPACds(PACds, group=NULL, cond1=NULL, cond2=NULL,
                      choosePA='apa',noIntergenic=TRUE, totPACtag=minSumPAT)
  }

  if (method %in% c('DESEQ2')) {
    geneDE=DESeqNoInteractionAll(dds=dds, factor1=group, factor2=NULL, norm=FALSE)
    geneDERes=digDESeqNoInteractionAll(geneDE, factor1=group, factor2=NULL, ref1=NULL, ref2=NULL, pairwise=TRUE, series=FALSE, levels1=NULL, levels2=NULL)
    condPairs=.getCondPairs(PACds, group)
    geneDEHeat=getHeatmapRes(object=geneDERes, colData=condPairs, padjThd=2, valueThd=0)

    geneDERes=new("movDEPACRes", group=group, method=method, conds=condPairs, pairwise=list(padj=geneDEHeat@padj, value=geneDEHeat@value))
    return(geneDERes)

    #condPairs=.getCondPairs(colnames(geneDERes@pairwise[[1]]), group)
  } else if (method=='DEXSEQ') {
    geneDE=DEXSeqOneFactor(dds, factor1=group, reflevel='', alllevel=TRUE, norm=FALSE, returnDxr=FALSE)
    geneDERes=DEXSeq::DEXSeqResults( geneDE )

    condPairs=.getCondPairs(PACds, group)
    geneDEHeat=getHeatmapRes(object=geneDERes, colData=condPairs, padjThd=2, valueThd=0)

    PACusage=NULL
    #TODO: check whether DEXseq's logFC is correct
    PACusage=getHeatMtx_DEXSeqResults(geneDERes, conds=unique(unlist(condPairs)), padjThd=2, valueThd=0)
    PACusage=cbind(PACusage@value,padj=PACusage@padj[,1])

    geneDERes=new("movDEPACRes", group=group, method=method,  conds=condPairs, pairwise=list(padj=geneDEHeat@padj, value=geneDEHeat@value), PACusage=PACusage)
    return(geneDERes)
  }  else if (method=='CHISQ') {

    condPairs=.getCondPairs(PACds, group)

    for (i in 1:nrow(condPairs) ) {
      c1=as.character(condPairs[i, 1])
      c2=as.character(condPairs[i, 2])
      p=subsetPACds(PACds, group=group, cond1=c1, cond2=c2, pool=TRUE, choosePA='apa')

      d=cbind(gene=p@anno$gene, .asDf(p@counts))
      d=d[order(d$gene),]

      ## chisqDEPAC for each gene
      sp=split(d[, -1], d[, 'gene'], drop=T)
      sp=lapply(sp, chisq_depac)
      pvGene=unlist(lapply(sp, '[[', 1))
      if (chisqPadjust) pvGene=p.adjust(pvGene)
      pvGene=1-pvGene ## set value as 1-pv, so that can use movStat(valueThd=0.95..) to filter DEPACs

      nPAC=Biobase::listLen(sp)-1
      pvGene=rep(pvGene, nPAC)

      sp=lapply(sp, function(par) par=par[-1])
      pvRows=unlist(sp)
      if (chisqPadjust) pvRows=p.adjust(pvRows)
      pvRows=cbind(PAC=rownames(d), pv=pvRows)
      pvGene=cbind(PAC=rownames(d), pv=pvGene)
      colnames(pvRows)[2]=paste0(c1,'.',c2)
      colnames(pvGene)[2]=paste0(c1,'.',c2)

      if (i==1) {
        padj=pvRows
        value=pvGene
      } else {
        padj=merge(padj, pvRows, by.x='PAC', by.y='PAC', all.x=TRUE, all.y=TRUE)
        value=merge(value, pvGene, by.x='PAC', by.y='PAC', all.x=TRUE, all.y=TRUE)
      }
  } #~for
    rn=padj[,'PAC']
    padj=padj[,-which(colnames(padj)=='PAC'), drop=F]
    value=value[,-which(colnames(value)=='PAC'), drop=F]
    padj=apply(padj, 2, as.numeric)
    value=apply(value, 2, as.numeric)
    padj=.asDf(padj)
    value=.asDf(value)
    rownames(padj)=rn
    rownames(value)=rn
    geneDERes=new("movDEPACRes", group=group, method=method, conds=condPairs,
                  pairwise=list(padj=padj, value=value))
    return(geneDERes)
  }
}

# get chisq.pv for all and for each row (peaks or PA)
# as in Shulman 2019 et al (test_apa + test_peaks)
# y: a matrix with rownames and colnames; row=peaks, col=clusters
# return: vector [totalPV, eachPV1, eachPV2...], totalPV is the whole chisq.pv
# Ex.
# y <- matrix(c(0, 0, 800, 1000, 50, 10), ncol = 2, dimnames = list(c('peak1','peak2', 'peak3'), c('c1','c2')))
# chisq_depac(y)
chisq_depac<-function(y) {
  pv=suppressWarnings(chisq.test(y)$p.value)
  #pv=suppressWarnings(fisher.test(y)$p.value) #to make it same as movAPAswitch
  expected <- colSums(y)/sum(y)
  chi.fit <- function(x) {
      suppressWarnings(chisq.test(p = expected, x = x)$p.value)
  }
 pvs=apply(y, 1, chi.fit)
 return(c(pv, pvs))
}

#' Analyze 3' UTR shortening/lengthening.
#'
#' movUTRtrend uses three methods for analysis of 3' UTR shortening/lengthening (or called 3' UTR switching).
#' It performs analyses for all condition pairs of the given group. If you only need a specific condition pair, call
#' subsetPACds first to filter used conditions.
#'
#' For linearTrend, switching events are not per-filtered by padj value, the change column (-1 or 1) is determined by cor value.
#' For DE/DEX, it will first filter significant switching paris by fisherThd and logFCThd, and then select one pair by selectOne.

#' @param PACds a PACdataset.
#' @param group a sample group name to get all conditions under the group.
#' @param method can be linearTrend, DE, DEX.
#' @param avgPACtag filter PACs with average tag number between cond1-cond2 >= avgPACtag.
#' @param avgGeneTag the same as avgPACtag but for gene.
#' @param DEPAC.padjThd used only for method=DE/DEX, to filter DE PACs with padj<DEPAC.padjThd.
#' @param aMovDEPACRes used only for method=DE/DEX, an object of movDEPACRes to get DE PAC results.
#' @param mindist mindist and fisherThd,logFCThd Used to filter all switching PAC pairs.
#' @param fisherThd used to filter DE PAC list from DE/DEX method by pvalue<fisherThd.
#' @param logFCThd used to filter |rcij| >=logFCThd, where |rcij| is the log2(ratio of PA1 between two samples/ratio of PA2).
#' @param selectOne used only for method=DE/DEX, can be 'farest','logFC','fisherPV' to filter one switching pair.
#' @return A movUTRTrendRes object, having the following slots (type, method, group, conds, pairwise=list(padj, value), fullList)
#' Here @fullList=[cond1.cond2, cond3.cond4, ...].
#' The full columns are: gene, nPAC, geneTag1, geneTag2, avgUTRlen1, avgUTRlen2, pvalue, padj,
#'   change <linear: cor, logRatio>, PAs1, PAs2.
#' @examples
#' data(PACds)
#' ## Detect 3'UTR lengthening/shortening events using the linear trend method.
#' ## Only PACs and genes with average read count between the two conditions >=10 and >=20 are used.
#' utr=movUTRtrend(PACds, group='group', method='linearTrend', avgPACtag=10, avgGeneTag=20)
#' ## Number of genes for analyzing, including those not significant.
#' lapply(utr@fullList, nrow)
#' head(utr@fullList[["anther.embryo"]])

#' ## Detect 3'UTR lengthening/shortening events using the DEX method.
#' ## First get DE PAC results by DEXseq.
#' DEXPAC=movDEPAC(PACds, method='DEXseq', group='group', minSumPAT=10)
#' ## Then get 3UTR switching events.
#' swDEX=movUTRtrend(PACds, group='group', method='DEX',
#'                   avgPACtag=10, avgGeneTag=20,
#'                   aMovDEPACRes=DEXPAC, DEPAC.padjThd=0.01,
#'                   mindist=50, fisherThd=0.01, logFCThd=1, selectOne='farest')
#' ## Number of genes for analyzing, including those not significant.
#' lapply(swDEX@fullList, nrow)
#' head(swDEX@fullList[["anther.maturePollen"]])

#' @name movUTRtrend
#' @family comparison functions
#' @export
movUTRtrend<-function(PACds, group, method='linearTrend',
                      avgPACtag=5, avgGeneTag=10,   #for cond1&cond2, not for all
                      aMovDEPACRes=NULL, DEPAC.padjThd=0.1,  # DE/DEX: to filter DEPAC
                      mindist=50, fisherThd=0.1, logFCThd=1, selectOne='farest') { # DE/DEX: to filter switching


  PACds@counts=asAnyMatrix(PACds@counts)
  method=tolower(method)
  if (method=='linear') method='lineartrend'
  if (!(method %in% c('lineartrend','de','dex'))) {
    stop(cat(method, 'not in linearTrend','DE','DEX\n'))
  }

  if (method!='lineartrend' & is.null(selectOne)) {
    stop("method is DE/DEX, must set selectOne as farest/logFC/fisherPV to choose one switching-PAC-pair")
  }

  if (method=='lineartrend' & !is.null(aMovDEPACRes)) stop("method=lineartrend but aMovDEPACRes!=NULL")
  if (!is.null(aMovDEPACRes)) {
    if (!inherits(aMovDEPACRes, 'movDEPACRes')) stop("aMovDEPACRes must be of class movDEPACRes")
    if (method=='de' & aMovDEPACRes@method!='DESEQ2') stop("method=de but aMovDEPACRes@method!=DESEQ2")
    if (method=='dex' & aMovDEPACRes@method!='DEXSEQ') stop("method=dex but aMovDEPACRes@method!=DEXSEQ")
  }

  condPairs=.getCondPairs(PACds, group)
  rs=list()
  for (condname in rownames(condPairs)) {
    c1=condPairs[condname, 1]
    c2=condPairs[condname, 2]
    cat(condname,'\n')
    if (method=='lineartrend') {
      rs[[condname]]=UTRtrend_linear(PACds, group=group, cond1=c1, cond2=c2, avgPACtag=avgPACtag, avgGeneTag=avgGeneTag)

    } else if (method=='de') {
      rs[[condname]]=UTRtrend_DE(PACds, group, c1, c2, avgPACtag=avgPACtag, avgGeneTag=avgGeneTag,
                                 DEpadjThd=DEPAC.padjThd, aMovDEPACRes=aMovDEPACRes, # DE/DEX: to filter DEPAC
                                 mindist=mindist, fisherThd=fisherThd, logFCThd=logFCThd, selectOne=selectOne)

    } else if (method=='dex') {
      rs[[condname]]=UTRtrend_DE(PACds, group, c1, c2, avgPACtag=avgPACtag, avgGeneTag=avgGeneTag,
                                 DEpadjThd=DEPAC.padjThd, aMovDEPACRes=aMovDEPACRes, # DE/DEX: to filter DEPAC
                                 mindist=mindist, fisherThd=fisherThd, logFCThd=logFCThd, selectOne=selectOne)
    }
  }


  rs=new("swResults",sw=rs,type='UTRtrend',method=method)
  heat=getHeatmapRes(object=rs, colData=condPairs, padjThd=NULL, valueThd=NULL)

  res=new("movUTRTrendRes", group=group, method=method, conds=condPairs, pairwise=list(padj=heat@padj, value=heat@value), fullList=rs@sw)
  return(res)
}


#' Analyze APA site switching.
#'
#' movAPAswitch analyzes APA site switching for genes with both canonical and non-canonical APA sites.
#' It performs analyses for all condition pairs of the given group. If you only need a specific condition pair, call
#' subsetPACds first to filter used conditions.
#'
#' This function is similar to movUTRtrend(method=DE/DEX), but it can be used for non-3'UTR PACs.
#' If nDEPAC>0, then first get DE PACs from aMovDEPACRes. For genes with more than two PACs, each pair of PACs with distance>=mindist is computed for switching.
#' If a pair of PACs meets the following criteria, then it is switching:
#' (1) DEPAC# (by DEPAC.padjThd) >= nDEPAC; (2) Pvalue of the Fisher's text of the two PACs between conditions < fisherThd;
#' (3) |Relative usage (RU)| of the two PACs between conditions >=logFCThd; Given PAi and PAj between samples A and B, RU=log2(ejB/ejA)-log2(eiB/eiA) (4) If cross=T, then expression levels of the two PACs are switched between conditions.
#'
#' @param PACds a PACdataset.
#' @param group a sample group name to get all conditions under the group.
#' @param avgPACtag filter PACs with average tag number between cond1-cond2 >= avgPACtag.
#' @param avgGeneTag the same as avgPACtag but for gene.
#' @param only3UTR If TRUE, then first filter 3UTR PACs.
#' @param mergeReps pool or avg replicates of the same condition, default is pool.
#' @param aMovDEPACRes an object of movDEPACRes to get DE PAC results. If nDEPAC>0, should provide.
#' @param DEPAC.padjThd to filter DE PACs with padj<DEPAC.padjThd.
#' @param nDEPAC can be NULL (0), 1 or 2, denoting that the switching pair has 0, 1 or 2 DE PACs.
#'               if nDEPAC>0 then should provide aMovDEPACRes.
#' @param mindist the min distance between two PACs to be considered for switching analyses.
#' @param fisherThd used to filter DE PAC list from DE/DEX method by pvalue<fisherThd.
#' @param logFCThd used to filter |rcij| >=logFCThd, where |rcij| is the log2(ratio of PA1 between two samples/ratio of PA2).
#' @param cross Whether the change of ratio should be crossing for PA1 and PA2.
#' @param selectOne can be 'farest','logFC','fisherPV' to output only one switching pair if a gene has multiple switching pairs.
#' @return An object of movAPASwitchRes similar to movUTRTrendRes.
#' @examples
#' # Detect 3UTR APA switching using DEXseq, which is similar to movUTRtrend(method=DEX)
#' data(PACds)
#' ## First get DE PAC results by DEXseq.
#' DEXPAC=movDEPAC(PACds, method='DEXseq', group='group', minSumPAT=10)
#' ## Then get 3UTR switching genes, usig selectOne=NULL to detect all pairs of switching PACs.
#' swDEX=movAPAswitch(PACds, group='group',aMovDEPACRes=DEXPAC,
#'                    avgPACtag=5, avgGeneTag=10,
#'                    only3UTR=TRUE,
#'                    DEPAC.padjThd=0.1, nDEPAC=1,
#'                    mindist=50, fisherThd=0.1, logFCThd=0.5, cross=FALSE, selectOne=NULL)
#'
#' lapply(swDEX@fullList,nrow)
#' head(swDEX@fullList[["anther.embryo"]])
#' head(swDEX@pairwise$padj
#'
#' ## Detect APA switching events involving non-3UTR PACs.
#' ## First get DE PAC results by DESeq2 (or DEXseq).
#' DEXPAC=movDEPAC(PACds, method='DEXseq', group='group', minSumPAT=10)
#' ## Then detect APA switching events,
#' ## using selectOne=NULL to get all pairs of switching PACs.
#' swDE=movAPAswitch(PACds, group='group',aMovDEPACRes=DEXPAC,
#'                  avgPACtag=10, avgGeneTag=20,
#'                  only3UTR=FALSE,
#'                  DEPAC.padjThd=0.1, nDEPAC=1,
#'                  mindist=50, fisherThd=0.1, logFCThd=0.5,
#'                  cross=FALSE, selectOne=NULL)
#' ## Stat the switching results.
#' stat=movStat(object=swDE, padjThd=0.1, valueThd=1)
#' stat$nsig

#' @name movAPAswitch
#' @family comparison functions
#' @export
movAPAswitch<-function(PACds, group, avgPACtag=5, avgGeneTag=10,
                       only3UTR=FALSE, mergeReps='pool',
                       aMovDEPACRes=NULL, DEPAC.padjThd=NULL, nDEPAC=NULL, # to filter DE/DEXPAC
                       mindist=50, fisherThd=0.1, logFCThd=1, cross=FALSE, selectOne=NULL) { #to filter switching pair

 if (is.null(nDEPAC)) nDEPAC=0

  if (nDEPAC>0 & is.null(aMovDEPACRes)) {
    stop("nDEPAC>0, need to provide aMovDEPACRes to get DE PAC(s)!")
  }
  if (nDEPAC>0 & is.null(DEPAC.padjThd)) {
    stop("nDEPAC>0, need to provide DEPAC.padjThd to filter DE PAC(s) from aMovDEPACRes!")
  }

  condPairs=.getCondPairs(PACds, group)
  rs=list()
  for (condname in rownames(condPairs)) {
    c1=condPairs[condname, 1]
    c2=condPairs[condname, 2]
    cat(condname,'\n')
    rs[[condname]]=APAswitching_DE(PACds=PACds, group=group, cond1=c1, cond2=c2, avgPACtag=avgPACtag, avgGeneTag=avgGeneTag,
                                   only3UTR=only3UTR, mergeReps='pool',
                                   aMovDEPACRes=aMovDEPACRes, DEpadjThd=DEPAC.padjThd, nDEPAC=nDEPAC,
                                   mindist=mindist, fisherThd=fisherThd, logFCThd=logFCThd, cross=cross, selectOne=selectOne)
  }
  if (length(rs)==0) return(NULL)

  if (is.null(aMovDEPACRes)) {
    method='de'
  }  else {
    method=aMovDEPACRes@method
  }
  rs=new("swResults",sw=rs, type='APAswitching', method=method)
  heat=getHeatmapRes(object=rs, colData=condPairs, padjThd=NULL, valueThd=NULL)

  res=new("movAPASwitchRes", group=group, method=method, conds=condPairs, pairwise=list(padj=heat@padj, value=heat@value), fullList=rs@sw)
  return(res)
}

#' Get APA index
#'
#' movAPAindex calculates APA index by different metrics for each condition.
#'
#' This funtion uses all samples to filter APA sites, and then gets index for each sample.
#' So there might be NaN for some samples without enough PAC (tag=0).
#' The final result may have NaN but will not have Inf (it has been assigned the max value for Inf element).
#' For method=WUL, NaN happens when the PAC expression in sample A for gene X is 0.
#'
#' For method=RUD or smartRUD (distal/gene), NaN happens when the PAC expression in sample A for gene X is 0.
#'
#' For method=SLR (short/long), NaN happens when both the proximal and distal PAC is 0; Inf happens when distal PAC=0 but proximal not 0.
#'
#' For method=smartRUD, the PACds should be the result of \code{\link{get3UTRAPApd}}.
#'
#' @param PACds a PACdataset.
#' @param method Can be WUL, RUD, SLR, smartRUD.
#' \itemize{
#' \item{"RUD"}: Relative Usage of Distal Poly(A) Sites (Ji 2009) (3UTRAPA with/without non-3UTR PACs, choose2PA=NULL/PD/MOST)
#' \item{"smartRUD"}: a smarter RUD index. It first calls \code{\link{get3UTRAPApd}} to obtain proximal and distal PAs in a smarter way.
#' \item{"WUL"}: Weighted UTR length (Jia 2017; Ultisky) (3UTRAPA, choose2PA=NULL/PD/MOST)
#' \item{"SLR"}: short to long ratio (Begik 2017) (3UTRAPA, choose2PA=PD/most). Similar to RED (relative expression difference) (W Li 2015; Arefeen 2018). RED is equivalent to divide the SLR index but just choose PD/most PA.
#' \item{"RS"}: Patrick et al. 2016, sum(abs_ratio/2) (3UTRAPA if choose2PA!=null; allAPA if =NULL)
#' }
#' @param choose2PA specify whether and how to choose two PACs. The value can be NULL, PD (choose only proximal and distal sites),
#' farest (choose two PACs that are with the longest distance), most (choose two PACs with the most abundance).
#' @param RUD.includeNon3UTR only used when method=RUD, to use all PA or just 3UTR PA.
#' @param clearPAT set PAC counts to 0 when it is <= clearPAT. This is important to deal with single-cell data.
#' For example, if a gene only has 1 read in a cell, then the RUD of its one pA with 1 read would be 1 (which is not very reasonable).
#' Using clearPAT=1, the gene would be set as 0 count, then its pA's RUD would be NaN.
#' Using clearPAT>=1, some genes may be removed in the final APA index list due to its both pAs are all 0s after setting clearPAT.
#' @param sRUD.oweight default is FALSE, if T and method=smartRUD, then output a list including RUD and weights matrix. Otherwise, only output a matrix of RUD.
#' This option is used to remove the effect of low expressed PACs. Otherwise, PAC with tag=1 may have high impact on RUD/SLR calculations.
#' After clearPAT, some results will be set as NaN.
#' @return A matrix/dgCMatrix with rownames denoting genes, columns are samples, and values are the APA index values. Also includes other method-specific columns.
#' @examples
#' ## smart RUD
#' data(scPACds)
#' PACds=get3UTRAPApd(scPACds)
#' ## set counts in genes of individual samples with only 1 total PATs as 0,
#' ## then the RUD of pAs in these genes would be Nan.
#' ## if not set clearPAT=1, then if one gene has only total 1 read in a sample,
#' ## then a PA with 1 read would be RUD=1!
#' rudw=movAPAindex(PACds, method="smart", sRUD.oweight=TRUE, clearPAT=1)
#' head(rudw$rud[, 1:2]);  head(rudw$weight[, 1:2])
#' rud=movAPAindex(PACds, method="smart", sRUD.oweight=F)
#' head(rud[, 1:2])
#'
#' ## WUL
#' data(PACds)
#' gs=movAPAindex(PACds, method="WUL", choose2PA=NULL)
#' gs=movAPAindex(PACds, method="SLR", choose2PA='PD', clearPAT=5)
#'
#' ## GPI index for only proximal PAC
#' geneGPI2=movAPAindex(PACds, method="GPI", choose2PA=NULL) #use all 3UTR PACs
#'
#' geneGPI=movAPAindex(PACds, method="GPI", choose2PA='PD') #use only proximal&distal PACs
#' PACds1=get3UTRAPAds(PACds, sortPA=TRUE, choose2PA='PD')
#' geneGPI1=movAPAindex(PACds1, method="GPI", choose2PA='PD')
#' identical(geneGPI1, geneGPI) ## SAME
#' ## use movPAindex to get GPI (but also including distal GPI)
#' paGeo=movPAindex(PACds1, method='geo')
#' @name movAPAindex
#' @family comparison functions
#' @export
movAPAindex <- function (PACds, method,
                         choose2PA=NULL, RUD.includeNon3UTR=FALSE,
                         clearPAT=0, sRUD.oweight=FALSE) {

  PACds@counts=asAnyMatrix(PACds@counts)

  method=toupper(method)
  if (!(method %in% c("WUL","RUD",'SLR','GPI','SMART','SM','SMARTRUD','SRUD'))) {
    stop("method must be WUL/RUD/SLR/GPI/smartRUD!")
  }

  if (clearPAT>=1) {
    cat("clearPAT>0: set elements in @counts <=", clearPAT, "to 0!\n")
    ## cannot use subset, because will remove not P/D-paired lines
    ##PACds=subsetPACds(PACds, clearPAT=clearPAT, verbose=TRUE)
    PACds@counts[PACds@counts<=clearPAT]=0
    invisible(gc())
    #idx=rowSums(pacds@counts)!=0
    #pacds=pacds[idx]
  }

  ## smartRUD from get3UTRAPApd's PACds
  .getRUDperGene<-function(pacdsPD, sRUD.oweight=FALSE) {

    if (!AinB(c('pdWhich','gene'), colnames(pacdsPD@anno))) stop("movAPAindex (smartRUD): please run get3UTRAPApd first to get 3UTR proximal and distal PACds!")

    #pacdsPD@counts=pacdsPD@counts+psudo
    did=which(pacdsPD@anno$pdWhich=='D')
    pid=which(pacdsPD@anno$pdWhich=='P')

    rud=pacdsPD@counts[did, ]/(pacdsPD@counts[pid, ]+pacdsPD@counts[did, ])
    rownames(rud)=pacdsPD@anno$gene[did]

    if (sRUD.oweight) {
      weights=round(pacdsPD@counts[pid, ]+pacdsPD@counts[did, ])
      rownames(weights)=pacdsPD@anno$gene[did]
      cat("movAPAindex (smartRUD): sRUD.oweight=TRUE, output a list(rud, weight)\n")
      return(list(rud=rud, weight=weights))
    } else {
      return(rud)
    }
  }


  if (method %in% c('SMART','SM','SMARTRUD','SRUD')) {
    mehod='SMART'
    gs=.getRUDperGene(PACds, sRUD.oweight=sRUD.oweight)
  }

  if (!is.null(choose2PA)) RUD.includeNon3UTR=FALSE

  if (method %in% c('SLR')) {
    if (is.null(choose2PA)) {
      stop("Error: method=SLR, should provide choose2PA=PD/most")
    }
  }

  if (method %in% c('WUL')) { #avg UTR len

    PACds=get3UTRAPAds(PACds, sortPA=TRUE, choose2PA=choose2PA)
    gs=geneStatForPACds(PACds, avgUTRlen=TRUE, geneTag=FALSE, PAinfo=FALSE, merge=TRUE)
    rownames(gs)=gs$gene
    gs$gene<-NULL

  } else if (method %in% c('RUD')) { #distal ratio

    if (!RUD.includeNon3UTR) {
      PACds=get3UTRAPAds(PACds, sortPA=TRUE, choose2PA=choose2PA)
    }

    #gene tagnum
    gs=geneStatForPACds(PACds, avgUTRlen=FALSE, geneTag=TRUE, PAinfo=FALSE, merge=TRUE)
    gs=.asDf(gs)
    rownames(gs)=gs$gene
    gs$gene=gs$nPAC=NULL
    gs=as.matrix(gs)

    #distal PA
    distPACds=subsetPACds(PACds, group=NULL, cond1=NULL, cond2=NULL,
                          avgPACtag=0, avgGeneTag=0, choosePA='distal',
                          noIntergenic=TRUE, avg=FALSE, verbose=FALSE)
    if (nrow(gs)<nrow(distPACds@counts)) {
      stop("Error: distal PAnum > nrow gs!")
    }
    rownames(distPACds@counts)=distPACds@anno$gene
    gs=gs[rownames(distPACds@counts),]
    gs=distPACds@counts/gs

  } else if (method %in% c('SLR')) { #short to long ratio

    PACds=get3UTRAPAds(PACds, sortPA=TRUE, choose2PA=choose2PA) #short, long
    gs=PACds@counts[seq(1,nrow(PACds@counts),2),]/PACds@counts[seq(2,nrow(PACds@counts),2),]
    rownames(gs)=unique(PACds@anno$gene)

  } else if (method %in% c('GPI')) { #geometric proximal index

    PACds=get3UTRAPAds(PACds, sortPA=TRUE, choose2PA=choose2PA) #3UTR all PAs
    gs=getAPAgeo(PACds, chooseAPA=FALSE, psudo=1, sort=FALSE)
    # get the first PA as proximal
    r=rle(PACds@anno$gene)$lengths
    id=c(1, cumsum(r)+1)
    id=id[1:(length(id)-1)]
    gs=gs[id, ]
    rownames(gs)=unique(PACds@anno$gene)
  }

  if (!inherits(gs, 'list')) {
    #Inf, Nan: set Inf for max value
    if (sum(is.infinite(as.matrix(gs)))>0) {
      gs=apply(gs, 2, function(x) {m=max(x[!is.na(x) & !is.infinite(x)]); x[is.infinite(x)]=m; return(x)})
    }
  }

  #Remove all NaN rows (due to clearPAT)
  if (inherits(gs, 'list')) { #smartRUD with weights
    nc=ncol(gs[[1]])
    rid=rowSums(is.na(gs[[1]]))==nc
    if (sum(rid)>0) {
      cat("Removed", sum(rid), "all NaN rows (genes), probably due to clearPAT>0.\n")
      gs[[1]]=gs[[1]][-rid, ]
      gs[[2]]=gs[[2]][-rid, ]
    }
  } else {
    asAnyMatrix(gs)
    nc=ncol(gs)
    rid=rowSums(is.na(gs))==nc
    if (sum(rid)>0) {
      cat("Removed", sum(rid), "all NaN rows (genes), probably due to clearPAT>0.\n")
      gs=gs[-rid, ]
    }
  }

  return(gs)
}


#' Calculate PAC index for each PAC in each condition.
#'
#' movPAindex calculates index of each PAC, while movAPAindex calculates index for APA genes.
#' method=ratio/geo is only allowed for APA genes, so these indices are only calculated for PACs in APA genes.
#'
#' @param PACds a PACDataset object.
#' @param method should be "shannon","ratio", "geo".
#' \itemize{
#' \item{"shannon"}: tissue specificity;
#' \item{"ratio"}: ratio for APA sites (single PA genes are removed first);
#' \item{"geo"}: PUI by geometry mean (=log2(pA/all_PAs_geo_mean)), as in Shulman 2019 et al (single PA genes are removed first).
#' }
#' @return A matrix of PA index.
#' \itemize{
#' \item{shannon}: return a dataframe of [H (overall TS), Q_min (min TS), Q_min_tissue, Q for each tissue], in which smaller H or Q means higher tissue-specificity.
#'  If PA=0 in one sample, then the respective Qval is NA.
#' \item{"geo/ratio"}: return a dataframe of the same row# and col# of the PACds@counts.
#' }
#' @param shan.ratio Pass to specificityByShannon.
#' @examples
#' data(PACds)
#' PAindex=movPAindex(PACds, method='ratio')
#' PAindex=movPAindex(PACds, method='geo')
#' PAindex=movPAindex(PACds, method='shannon')
#' @name movPAindex
#' @family comparison functions
#' @export
movPAindex<-function(PACds, method, shan.ratio=FALSE) {
  method=tolower(method)
  if (method=='shan') method='shannon'
  if (!method %in% c("shannon","ratio","geo")) {
    stop("method must be shannon/ratio/geo")
  }
  if (method=='shannon') {
    PAindex=specificityByShannon(PACds, shan.ratio=shan.ratio)
  } else if (method=='ratio') {
    PAindex=getAPAratio(PACds, chooseAPA=TRUE)
  } else if (method=='geo') {
    PAindex=getAPAgeo(PACds, chooseAPA=TRUE)
  } else if (method=='dex') {
    ## DEX: exon usage from DEXSeq (not realized) #TODO
    cat("please get dex index using DEXseq\n")
    PAindex=NULL
  }
  return(PAindex)
}


#' Get APA index difference between conditions
#'
#' movAPAindexDiff calculates APA index difference between conditions by different metrics.
#'
#' This function is different from APAswitching/UTR_trend. It does not filter PACds by cond1/2 but for all conditions.
#' So there might be NaN/Inf for some condition pairs when some PACs are not expressed in two specific condtions but expressed in other conds.
#' @param PACds a PACdataset.
#' @param group a sample group name to get all conditions under the group.
#' @param method can be WUL, RUD, SLR.
#' \itemize{
#' \item{"WUL"}: Weighted UTR length (Jia 2017; Ultisky) (3UTRAPA, choose2PA=NULL/PD/MOST)
#' \item{"RUD"}: Relative Usage of Distal Poly(A) Sites (Ji 2009) (3UTRAPA with/without non-3UTR PACs, choose2PA=NULL/PD/MOST)
#' \item{"SLR"}: short to long ratio (Begik 2017) (3UTRAPA, choose2PA=PD/most). Similar to RUD (relative expression difference) (W Li 2015; Arefeen 2018). RED is equivalent to divide the SLR index but just choose PD/most PA.
#' \item{"RS"}: Patrick et al. 2016, sum(abs_ratio/2) (3UTRAPA if choose2PA!=null; allAPA if =NULL)
#' }
#' @param choose2PA can be NULL, PD (proximal, distal, farest), most (most abundant two PAs).
#' @param doReps Can be `avg` or `pool` to merge replicates.
#' @param RUD.includeNon3UTR only used when method=RUD, to use all PA or just 3UTR PA.
#' @param clearPAT set PAC counts to 0 when it is < clearPAT.
#' First call movAPAindex(clearPAT), if method=RS then subsetPACds(clearPAT).
#' If doReps!=NULL, then first doReps, and then filter by clearPAT.
#' @param nPermute If >0, then will get pvalue by bootstraping.
#' nPermute>0 must be with choose2PA=PD or MOST, because pvalue can be obtained for only two PACs.
#' @return A data frame with columns being condition pairs like [cond1.cond2], and values are log2(index1/index2).
#' @examples
#' data(PACds)
#' ## RED + only for PD
#' diff=movAPAindexDiff(PACds, group='len', method="SLR", choose2PA="PD", RUD.includeNon3UTR=FALSE)
#' ## RS+only for PD
#' diff=movAPAindexDiff(PACds, group='len', method="RS", choose2PA="PD")
#' ## Get pvalue
#' geneSLRsig=movAPAindexDiff(PACds, group='len', method='SLR',
#'                           choose2PA="PD", doReps='avg', clearPAT=0, nPermute=10)
#' @name movAPAindexDiff
#' @family comparison functions
#' @export
## TODO: check pvalue with Liu Tao's results and whether the result is resonable.
movAPAindexDiff <- function (PACds, group, method, doReps='avg',
                             choose2PA=NULL, RUD.includeNon3UTR=FALSE, clearPAT=0, nPermute=0) {

  if (!(method %in% c("WUL","RUD",'SLR','RS'))) {
    stop("method must be WUL/RUD/SLR/RS!")
  }

  if (method %in% c('SLR')) {
    if (is.null(choose2PA)) {
      stop("Error: method=SLR (Short to Long Ratio), should provide choose2PA=PD(Proximal & Distal)/most")
    }
  }

  if (nPermute>0) {
    if(!(choose2PA %in% c('PD','MOST'))) stop("do permutation (nPermute>0), choose2PA must be PD/MOST!")
    if(!(method %in% c("RUD",'SLR','RS'))) stop("do permutation (nPermute>0), method must be RUD/SLR/RS!")
  }

  avg=TRUE; pool=FALSE
  if (doReps=='pool') {
    pool=TRUE; avg=FALSE;
  }

  #filter and pool
  PACds=subsetPACds(PACds, group=group, avg=avg, pool=pool, noIntergenic=TRUE, verbose=FALSE)

  #TODO: QAPA, Dapas, PD...

  condPairs=.getCondPairs(PACds, group)

  if(nPermute==0) {
    rs=matrix(nrow=0,ncol=0)

    #per sample index
    if (method %in% c('RUD')) {
      APAidx=movAPAindex(PACds, method=method,
                         choose2PA=choose2PA, RUD.includeNon3UTR=RUD.includeNon3UTR, clearPAT=clearPAT)
    } else if (method %in% c('WUL')) { #avg UTR len
      APAidx=movAPAindex(PACds, method=method, choose2PA=choose2PA, clearPAT=clearPAT)
    }  else if (method %in% c('SLR')) {
      APAidx=movAPAindex(PACds, method=method, choose2PA=choose2PA, clearPAT=clearPAT)
    } else if (method=='RS') { #ratio per sample
      if (!is.null(choose2PA)) {
        PACds=get3UTRAPAds(PACds, sortPA=TRUE, choose2PA=choose2PA, clearPAT=clearPAT)
      } else {
        #choose APA genes
        PACds=subsetPACds(PACds, choosePA='apa',
                          noIntergenic=TRUE, clearPAT=clearPAT)
      }
      PACds@counts=getAPAratio(PACds, chooseAPA=FALSE)
      PACds@anno=PACds@anno[rownames(PACds@counts),]
    }

    for (condname in rownames(condPairs)) {
      c1=condPairs[condname, 1]
      c2=condPairs[condname, 2]
      cat(condname,'\n')
      if (method!='RS') { #SLR RUD WUL
        diff=log2(APAidx[,c1]/APAidx[,c2])
      } else { #RS
        d=PACds@counts[,c(c1,c2)]
        d$gene=PACds@anno$gene
        #diff=plyr::ddply(d,.(gene), function(dd){sum(abs(dd[,1]-dd[,2]))/2})
        diff=d %>% dplyr::group_by(gene) %>% dplyr::summarise(r=sum(abs(.data[[c1]]-.data[[c2]]))/2)
        rowname=unlist(diff[,1])
        diff=unlist(diff[,2])
      }

      if(ncol(rs)==0) {
        rs=cbind(diff)
      } else {
        rs=cbind(rs,diff)
      }
      colnames(rs)[ncol(rs)]=condname
    }

    if (method!='RS') { #SLR RUD WUL
      rownames(rs)=rownames(APAidx)
    } else { #RS
      rownames(rs)=rowname
    }

    #Set Inf to row max, -Inf to row min
    if (sum(is.infinite(as.matrix(rs)))>0) {
      rs=apply(rs, 2, function(x)
      {  m=max(x[!is.na(x) & !is.infinite(x)]);
      mi=min(x[!is.na(x) & !is.infinite(x)]);
      x[x==Inf]=m; x[x==-Inf]=mi; return(x)}
      )
    }

    #remove a line if all NA
    rs=rs[rowSums(!is.na(rs))!=0, ]
    rs=.asDf(rs)

    #TODO: INF/NA/NaN/log2...puedo

    padj=rs; padj[]=NA
    rs=new("movIndexDiffRes", group=group, method=method, conds=condPairs, pairwise=list(padj=padj, value=rs))
    return(rs)

  } else { #permute

    PACds=get3UTRAPAds(PACds, sortPA=TRUE, choose2PA=choose2PA)
    rs=list()

    for (condname in rownames(condPairs)) {
      c1=condPairs[condname, 1]
      c2=condPairs[condname, 2]
      cat(condname,'\n')
      #format to vector: prox1, distal1, prox2, distal2
      pdv=PACds@counts[,c(c1,c2)]
      pdv=cbind(prox1=pdv[seq(1,nrow(pdv),2),1],dist1=pdv[seq(2,nrow(pdv),2),1],
                prox2=pdv[seq(1,nrow(pdv),2),2],dist2=pdv[seq(2,nrow(pdv),2),2])

      diff=apply(pdv, 1, saapOnce, method=method, n=nPermute)
      diff=t(diff)
      rownames(diff)=PACds@anno$gene[seq(1,nrow(PACds@anno),2)]
      colnames(diff)=c('realIndex', 'meanPmt', 'sdPmt', 'pv.ztest')
      diff=cbind(pdv, diff)
      diff=.asDf(diff)
      diff$gene=rownames(diff)
      rs[[condname]]=diff
    }

    #convert to swResults and then get heatmap
    rs=new("swResults", sw=rs, type='IndexDiff', method=method)
    heat=getHeatMtx_swResults(rs, condPairs, padjThd=2, valueThd=NULL, padjCol="pv.ztest", valueCol="realIndex", fixed=TRUE)
    res=new("movIndexDiffRes", group=group, method=method, conds=condPairs, pairwise=list(padj=heat@padj, value=heat@value), fullList=rs@sw)
    return(res)

  }#~else

}


# ---------------- *** viz *** ------------------

modelColorFTR='grey80' #gene ftr box color
modelColorTXDB="brown" #gene model txdb color
rawPAcolor='grey50' #PA point color
avgPAcolor='grey10' #avg PA point and line color
hightlightPAcolor='orange' #highlight PA point and line color

#Highlight condition, change background color only
.themeHighlightConds<-function() {
  theme(
    panel.background = element_rect(fill = "lightblue",
                                    colour = "lightblue",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                    colour = "white"),
    panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                    colour = "white")
  )
}

#' movViz generic
#'
#' movViz(object) is a generic function for plotting PAC distributions of a gene across selected samples.
#' @param ... Dot arguments passing to other functions.
#' @param object normally an object of movRes or PACdataset.
#' @family visualization functions
setGeneric("movViz", function(object, ...) standardGeneric("movViz"))

# ---------- movViz PACds ----------

#' Visualize a PACdataset
#'
#' movViz(PACdataset) plots PAC distributions of a gene across selected samples from a PACdataset object.
#' This function is actually used by many other movViz functions. PACs and conditions can be highlighted, average value will be shown.
#' PAC value can be represented by ratios instead of counts. PACs can be linked by lines to show trend.
#'
#' @param group same as `group` in subsetPACds.
#' @param highlightPAs the PAids to be highlighted in the plot.
#' @param otrack if FALSE, then return a list of plots otherwise return Tracks. Normally used for calling other movViz functions.
#' @param geneHeight default is 0.2. Set height ratio for the gene model track, the height of other data tracks are of the same height ratio: 1-geneHeight.
#' @param trackOrder a character vector to set the order of (part of) data tracks, which is useful to put some conditions at the top of the plot.
#' For example, if there are four tracks: condA, condB, condC, condD, when trackOrder=c('condC','CondB'), then the C and B tracks will be shown fist, followed by other tracks (not ordered).
#' @param simple if TRUE, then plot minimum elements: bar for avg; all tracks with the same Y limits; white background without grids.
#' @param ylimits default is c(NA,NA), set manual range for Y axis. This is applied for all tracks, and both values should be set.
#' @param showYaxis only valid when simple=TRUE. If showYaxis=FALSE, then do not show Y axis ticks, text, and line, and change X-axis color to grey.
#' parameters simple, ylimits, and showYaxis are particularly useful for controlling the plot of many conditions (e.g., single cell data) with collapseConds=FALSE.
#' @param showBox only valid when simple=TRUE; if simple=TRUE and showBox=TRUE, then show avg expression level of each PA as colored rectangle. The width of rectangle can be UPA_start~UPA_end or setted boxWidth.
#' @param boxWidth default is 20. When showBox=TRUE, boxwidth sets the width of the rectangle.
#' If boxWidth=NA and UPA_start/_end are in PACds@anno, then will use PAC range as xmin and xmax of the rectangle.
#' If boxWidth=NA and UPA_start/_end are not in PACds@anno, then will set boxWidth=20.
#' @param title add title to the plot at the top-left of the gene model track. If null, then no title.
#' @examples
#' \dontrun{
#' data(PACds)
#' load('Oryza_sativa.IRGSP-1.0.42.gff.rda')
#'
#' ## Visualize gene with PACs ##
#' gene='Os05g0451900'
#' PACds[PACds@anno$gene==gene, ]@anno
#' movViz(object=PACds, gene=gene, txdb=NULL, group='group', geneHeight=0.1)
#' movViz(object=PACds, gene=gene, txdb=gff, group='group', geneHeight=0.1)
#'
#' ## Visualize gene with DE PACs ##
#' ## Detect DE PACs using chisq-PACdsPAS method.
#' DEqPAC=movDEPAC(PACds, method='chisq', group='group', minSumPAT=20)
#' ## Visualize PACs of this gene in individual conditions.
#' ## Here the Y-axis is read count, the scale of which is different among conditions.
#' movViz(object=DEqPAC, gene=gene, txdb=gff, PACds=PACds, conds=NULL)
#' movViz(object=DEqPAC, gene=gene, txdb=NULL, PACds=PACds, conds=NULL)
#'
#' ## Only show DE PACs with padj<padjThd.
#' movViz(object=DEPAC, gene=gene, txdb=gff, PACds=PACds, collapseConds=FALSE, padjThd=0.01,
#'        showRatio=FALSE, showAllPA=TRUE)
#'
#' ## Show condition pairs in individual tracks.
#' ## If padjThd is given, then the DE PACs (padj<padjThd) will be highlighted (dashed yellow line).
#' movViz(object=DEqPAC, gene=gene, txdb=NULL, PACds=PACds,
#'        collapseConds=TRUE, padjThd=0.01, showRatio=T)
#' movViz(object=DEXPAC, gene=gene, txdb=gff, PACds=PACds,
#'        collapseConds=TRUE, padjThd=0.01, showPV=TRUE, showAllPA=TRUE)
#' movViz(object=DEPAC, gene=gene, txdb=gff, PACds=PACds,
#'      collapseConds=TRUE, padjThd=0.01, showPV=TRUE, showAllPA=TRUE, showRatio=F)
#'
#' movViz(object=DEPAC, gene=gene, txdb=gff, PACds=PACds, collapseConds=TRUE,
#'        padjThd=0.01, showPV=TRUE, showAllPA=FALSE, showRatio=F,
#'        conds=DEPAC@conds[c(1,3), ], highlightConds=DEPAC@conds[c(3), ], )
#' }
#'
#' @return A list or ggbio's Tracks. If PACds is empty, return NULL.
#' The first track is the gene model track and the height can be set by the geneHeight parameter.
#' If txdb is provided, then all transcripts in the gene model are plot.
#' If there are only intergenic PACs, then plot the intergenic; otherwise will remove intergenic PACs and only retain genomic PACs.
#' @describeIn movViz for a movRes object.
#' @export
setMethod("movViz", signature="PACdataset", def=function(object, gene, txdb=NULL,
                                                         group, conds=NULL, avg=FALSE, pool=FALSE, clearPAT=0,
                                                         showRatio=FALSE, linkPAs=FALSE,
                                                         highlightConds=NULL, highlightPAs=NULL, collapseConds=FALSE,
                                                         otrack=TRUE, geneHeight=0.2, trackOrder=NULL,
                                                         simple=FALSE, ylimits=c(NA,NA), showYaxis=TRUE,
                                                         showBox=FALSE, boxWidth=20, title=NULL) {

  if (collapseConds) simple=FALSE
  if (!simple) showBox=FALSE
  if (showBox) ylimits=c(NA,NA)

  PACds=subsetPACds(object, genes=gene, group=group, conds=conds, clearPAT=clearPAT, avg=avg, pool=pool)

  PACds@counts=as.matrix(PACds@counts)

  if (nrow(PACds@counts)==0) return(NULL)

  f=unique(PACds@anno$ftr)
  if (length(f)==1 & f[1]=='intergenic') {
    cat("Warning: shown are intergenic PAC(s)\n")
  } else if ('intergenic' %in% f) {
    PACds=subsetPACds(PACds, noIntergenic = TRUE)
  }

  #ratio
  if (showRatio) {
    PACds@counts=as.matrix(getAPAratio(PACds, chooseAPA = FALSE))
    if (nrow(PACds@counts)!=nrow(PACds@anno)) {
      stop("error getAPAratio() (showRatio=TRUE)")
    }
    PACds@anno=PACds@anno[rownames(PACds@counts), ]
    PACds@counts[is.na(PACds@counts)]=0
  }

  #gene model
  if (is.null(txdb)) {
    #If no GFF, then use ftr_start~ftr_end region.
    ext=0
    gs=min(PACds@anno$ftr_start)
    ge=max(PACds@anno$ftr_end)

    geneGR= GRanges(seqnames = PACds@anno$chr[1] ,
                    ranges =IRanges::IRanges(start=gs ,
                                    end=ge),
                    strand =as.character(PACds@anno$strand)[1],
                    transcript=PACds@anno$gene[1], model='exon')

    ##!!!! must add GenomicRanges::, otherwise return a list class, error raised!
    geneGR=GenomicRanges::split(geneGR, mcols(geneGR)$transcript)

    ## bug: remove all max.height =0.5, in autoplot
    #The whole region is box, each ftr is a line with text denoting 3UTR/intron..
    dftr=unique(PACds@anno[,c('ftr','ftr_start','ftr_end')])
    pg=ggbio::autoplot(geneGR, aes(type = model),  exon.rect.h=0.5,
                label.color = "black", color = modelColorFTR, fill = modelColorFTR) +
      ggplot2::geom_segment(aes(x=ftr_start, y=0.6, xend=ftr_end, yend=0.6, colour=ftr), data=dftr, size=5, lineend='square', alpha=0.6) +
      annotate("text", x = round((dftr$ftr_start+dftr$ftr_end)/2), y=0.65, label =dftr$ftr, colour='white') +
      viz_geneModel_theme  #Not show legend


  } else if (is.list(txdb)) {
    txdb=txdb$anno.need
    annoGene=txdb[txdb$gene_id==gene,]
    annoGene$type[annoGene$type=='intron']='gap' #only can shown utr/cds/exon/gap
    annoGene$type[annoGene$type=='CDS']='cds' #only can shown utr/cds/exon/gap
    #annoGene=annoGene[annoGene$type!='gap',]
    annoGene$type[annoGene$type %in% c('5UTR','3UTR', 'three_prime_UTR','five_prime_UTR')]='utr'
    annoGene=annoGene[annoGene$type %in% c('utr','gap','cds','exon'), ]
    geneGR= GenomicRanges::GRanges(seqnames =annoGene$seqnames ,
                    ranges =IRanges::IRanges(start=as.integer(annoGene$start) ,
                                    end=as.integer(annoGene$end)),
                    strand =as.character(annoGene$strand),
                    transcript=annoGene$Parent, model=annoGene$type)
    geneGR=GenomicRanges::split(geneGR, mcols(geneGR)$transcript) #Must convert to glist for each transcript, even only 1 transcript
    pg<-ggbio::autoplot(geneGR, aes(type = model), exon.rect.h=0.5,
                       label.color = "black", color = modelColorTXDB, fill = modelColorTXDB) +
      viz_geneModel_theme
  }

  #PA sites - triangles in annotation
  pas=PACds@anno$coord
  pg=pg + geom_point(aes(x=pas, y=0.5), shape=17, size=3, color='DarkViolet')

  if (!is.null(title)) {
    #add gene name at top-left of the gene model
    pg=pg + annotate("text", -Inf, Inf, label = title, hjust = 0, vjust = 1, color='red')
  }

  #PA expression: If multiple reps, then shown as points and show avg as vertical line
  #If only one rep, show as vertical line with a top point
  tks=list('Gene\nModel'=pg)

  #get counts
  pointDarkColors=RColorBrewer::brewer.pal(9, "Set1")
  pointLightColors=RColorBrewer::brewer.pal(9, "Pastel1")

  groups=as.character(PACds@colData[, group])
  conds=as.character(unique(groups))

  #detect outliers
  if (0) { #simple=TRUE
    dt=unlist(PACds@counts)
    outs=boxplot.stats(dt)$out

    if (length(outs)>0) {
      dt=dt[!(dt %in% outs)]
      mx=max(dt); mi=min(dt); v=mean(dt)
      up=sum(outs>=mx)
      if (up)  PACds@counts[PACds@counts>=mx]=mx
      up=sum(outs<=mi)
      if (up)  PACds@counts[PACds@counts<=mi]=mi
    }
  }

  if (collapseConds) {

    ## set ylim by ylimits or showRatio or simple
    ylims=c(NA,NA)
    if (sum(!is.na(ylimits))==2) {
      ylims=ylimits
    } else if (showRatio) {
      ylims = c(0, 1)
    }

    #manualCols the legend for point colors

    groupsCols=pointLightColors[1:length(conds)]
    names(groupsCols)=conds

    dtcond=.asDf(PACds@counts)
    #colnames(dtcond)=as.character(groups)
    dtcond=cbind(coord=PACds@anno$coord, dtcond)
    dtcond=reshape2::melt(dtcond, id='coord', variable.name = 'reps', value.name = 'nPAT')
    coldata=cbind(reps=rownames(PACds@colData), group=groups)
    dtcond=merge(dtcond, coldata, by='reps')

    #lighter points for reps
    p1 <-ggplot(data = dtcond)+ geom_point(aes(coord, nPAT, colour=group), data = dtcond)
    manualCols=groupsCols

    #hightlight PA -- plot a dashed vline for each hightlihgt PA
    hlPA=base::intersect(rownames(PACds@anno), highlightPAs)

    if (length(hlPA)>0) {
      hlcoord=PACds@anno[hlPA,'coord']
      p1=p1+geom_vline(xintercept=hlcoord, linetype="dashed", size=0.5, colour=hightlightPAcolor)
    }

    #darker point for avg
    dtcond$reps=NULL
    dtcond=dtcond %>% dplyr::group_by(group, coord) %>% dplyr::summarise(nPAT=mean(nPAT, na.rm=TRUE))
    dtcond$group=paste0(dtcond$group,'_avg')
    groupsCols[]=pointDarkColors[1:length(groupsCols)]
    names(groupsCols)=paste0(names(groupsCols),'_avg') #add _avg to use just one manualCols legend
    p1=p1 + geom_point(aes(coord, nPAT, colour=group), data = dtcond, size=3)

    #link avg points as line
    if (linkPAs)  p1=p1+geom_line(data = dtcond, aes(coord, nPAT, colour=group))

    #legend (TODO: how to align if legend positioned in the right)
    manualCols=c(manualCols, groupsCols)
    p1=p1+scale_colour_manual(name  = group, values=manualCols, limits = unique(groups))+ labs(y=NULL)

    p1=p1 + theme(legend.position='top',legend.title=element_blank(), legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0))

    p1 = p1 + scale_y_continuous(limits=ylims)

    tks=c(tks, list("PAC"=p1))

  } else {

    ## set ylim by ylimits or showRatio or simple
    ylims=c(NA,NA)
    if (sum(!is.na(ylimits))==2) {
      ylims=ylimits
    } else if (showRatio) {
      ylims = c(0, 1)
    } else if (simple) {
      ymax=max(PACds@counts)
      ymin=min(PACds@counts)
      ylims=c(ymin, ymax)
    }

    # if(showBox) {
    #   x=subsetPACds(PACds, group=group, avg = TRUE)
    #   boxScale=c(min(x@counts), max(x@counts))
    # }

    #each condition is one track
    i=TRUE
    for (cond in conds) {
      i=!i
      dtcond=.asDf(cbind(coord=PACds@anno$coord, as.matrix(PACds@counts[, which(groups==cond), drop=F])))
      addAvg=ncol(dtcond)>2
      dtcond$avg=rowMeans(as.matrix(dtcond)[,-1, drop=F], na.rm = TRUE)

      dtcond=reshape2::melt(dtcond, id='coord', variable.name = 'reps', value.name = 'nPAT')

      if (!simple) {
        #raw point
        p1 <-ggplot(data = dtcond)+ geom_point(aes(coord, nPAT), data = subset(dtcond, reps != 'avg'), colour=rawPAcolor)
        #avg line & point
        p1=p1+geom_point(data = subset(dtcond, reps == 'avg'), aes(coord, nPAT), colour=avgPAcolor, shape=20, size=3) +
          ggplot2::geom_segment(aes(x=coord,xend=coord,y=0,yend=nPAT), data = subset(dtcond, reps == 'avg'), colour=avgPAcolor,size=1)
        p1=p1 + labs(y=NULL)
      } else {
        #simple: same Y limit, not show Y ticks, not show grid and grey background, just bar for avg (not show avg point and raw points)
        boxHeight=5
        if (!showBox) { #show expression levels as black bar (only avg)

          p1=ggplot(data = dtcond)+ggplot2::geom_segment(aes(x=coord,xend=coord,y=0,yend=nPAT), data = subset(dtcond, reps == 'avg'), colour=avgPAcolor,size=1)
          p1=p1 + labs(y=NULL) + theme_classic()

        } else { #show expression levels as colored box (only avg)

          dtBox=.asDf(cbind(coord=PACds@anno$coord, as.matrix(PACds@counts[, which(groups==cond), drop=F])))
          dtBox$avg=rowMeans(as.matrix(dtBox)[,-1, drop=F], na.rm = TRUE)
          dtBox=dtBox[, c('coord','avg')]
          if (AinB(c('UPA_start','UPA_end'), colnames(PACds@anno), all=TRUE) & is.na(boxWidth)) {
            dtBox=cbind(dtBox, xmin=PACds@anno$UPA_start, xmax=PACds@anno$UPA_end)
          } else {
            if (is.na(boxWidth)) boxWidth=20
            dtBox=cbind(dtBox, xmin=dtBox$coord-floor(boxWidth/2), xmax=dtBox$coord+floor(boxWidth/2))
          }
          #print(dtBox)
          p1=ggplot() +
             ggplot2::geom_rect(data=dtBox, aes(xmin=xmin, xmax=xmax, ymin=0, ymax=5, group=coord, fill = avg))+
             scale_fill_gradientn(colours = c('white','#FFC6C6', 'red'), values=c(0,0.01,1), limits=ylims) +
             #scale_fill_gradient(low='#FFE2E2', high='red', limits=ylims) +
              guides(fill="none") #ylims
          p1=p1 + labs(y=NULL) + theme_classic()
        }

        if (!showYaxis) {
          p1= p1+ theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                        axis.line.y=element_blank(), axis.line.x =element_blank()) # element_line(color = "grey"))
        }
      }

      #highlight conds
      if (cond %in% highlightConds) {
        p1=p1 + .themeHighlightConds()
      }

      #highlight PA - change color and add vertical line
      hlPA=base::intersect(rownames(PACds@anno), highlightPAs)

      if (length(hlPA)>0) {
        hlcoord=PACds@anno[hlPA,'coord']
        dtcond1=subset(dtcond, coord %in% hlcoord & reps=='avg')
        p1=p1+geom_point(data = dtcond1, aes(coord, nPAT), colour=hightlightPAcolor, shape=20, size=3) +
          ggplot2::geom_segment(aes(x=coord,xend=coord,y=0,yend=nPAT), data = dtcond1, colour=hightlightPAcolor,size=1)
      } else if (length(highlightPAs)>0)  {
        cat(highlightPAs, "not in PACds, cannot hightlight!\n")
      }

      #link avg points
      if (linkPAs) {
        dtcond1=subset(dtcond, reps=='avg')
        p1=p1+geom_line(data = dtcond1, aes(coord, nPAT), colour=avgPAcolor)
      }

      if (sum(!is.na(ylims))==2 & !showBox) {
        #p1 = p1 + coord_cartesian(ylim=ylims, expand=FALSE)
        p1 = p1 + scale_y_continuous(limits=ylims)
      }
      tks=c(tks, list(p1))
    }
    names(tks)=c('Gene\nModel', conds)
  }#~collapse

  if (simple & !showYaxis & !showBox) {
    cat('simple=TRUE and showYaxis=FALSE: Y limits on the plot is', ylims,'\n')
  }
  if (showBox) cat('Rectangle color scale is', ylims,'\n')

  #order tracks
  if (length(tks)>2 & !is.null(trackOrder)) {
    trackOrder=cbind(order=1:length(trackOrder), cond=trackOrder)
    conds=cbind(cond=conds)
    trackOrder=merge(trackOrder, conds)
    if (nrow(trackOrder)>0) {
      cat('Ordering tracks by trackOrder...\n')
      trackOrder$order=as.numeric(trackOrder$order)
      trackOrder=trackOrder[order(trackOrder$order), , drop=F]
      conds=conds[,1]
      conds=c(trackOrder$cond, conds[!(conds %in% trackOrder$cond)])
      tks=tks[c(names(tks)[1], conds)]
    }
  }

  heights=c(geneHeight, rep((1-geneHeight)/(length(tks)-1), length(tks)-1))
  if (length(tks)>=20) names(tks)[2:length(tks)]=''

  if(otrack) return(ggbio::tracks(tks, heights=heights))
  return(tks)
}
)


#' Visualize a movRes object
#'
#' movViz(movRes) plots PAC distributions of a gene across selected samples from a movRes object.
#' This function actually calls movViz(PACdataset) for plot.
#'
#' @param object a movRes or PACdataset object.
#' @param gene a gene to plot.
#' @param txdb GFF annotation, a character of file name for .gff3 or .rda, or a list. See parseGenomeAnnotation().
#' @param conds a data frame of (cond1 cond2), specifying the conditions to plot.
#' @param avg used to subset PACds.
#' @param pool used to subset PACds.
#' @param clearPAT used to subset PACds.
#' @param showRatio If TRUE then show PAC ratio instead of PAC count.
#' @param linkPAs If TRUE, then plot lines to links PACs from the same sample.
#' @param highlightConds a two column data frame, denoting the highlight condtions.
#' Only used when collapseConds=TRUE.
#' @param collapseConds If TRUE, then collapse conds to one track, particularly useful for switching plot.
#' If FALSE, then each condition is one track.
#' @return A ggbio's Tracks.
#' @describeIn movViz for a movRes object.
#' @export
setMethod("movViz", signature="movRes", def=function(object, gene, txdb=NULL, PACds,
                                                     conds=NULL, avg=FALSE, pool=FALSE, clearPAT=0,
                                                     showRatio=FALSE, linkPAs=FALSE,
                                                     highlightConds=NULL, collapseConds=FALSE) {

  #library(ggbio, verbose = FALSE)

  pls=.movVizPreprocess(object=object, gene=gene, PACds=PACds,
                        conds=conds, avg=avg, pool=pool, clearPAT=clearPAT,
                        highlightConds=highlightConds, collapseConds=collapseConds)

  highlightConds=pls[['highlightConds']]
  if (!is.null(highlightConds)) highlightConds=paste0(highlightConds[,1], '.', highlightConds[,2])
  PACds=pls[['PACds']]
  condpairs=pls[['condpairs']]

  #If not collapse, each condition is one track, and highlight DEPA.
  if (!collapseConds) {

    tks=movViz(object=PACds, gene=gene, txdb=txdb, group=object@group, conds=unique(unlist(condpairs)),
               collapseConds=FALSE, linkPAs=linkPAs, showRatio=showRatio)

  } else {
    #Two conditions in one track.
    tks=NULL
    for (i in 1:nrow(condpairs)) {

      tks1=movViz(object=PACds, gene=gene, txdb=txdb, group=object@group, conds=unlist(condpairs[i, ]),
                  collapseConds=TRUE, linkPAs=linkPAs, otrack=FALSE, highlightPAs=NULL, showRatio=showRatio)

      if(is.null(tks1)) next

      #highligtht conds
      if(!is.null(highlightConds)) {
        if ((paste0(condpairs[i,1],'.',condpairs[i,2]) %in% highlightConds)) {
          tks1[[2]]=tks1[[2]]+.themeHighlightConds()
        }
      }

      if (is.null(tks)) {
        tks=tks1
        names(tks)=c('Gene\nModel', paste0(condpairs[i,1],'\n',condpairs[i,2]))
      } else {
        n1=names(tks)
        tks=c(tks, tks1[2])
        names(tks)=c(n1, paste0(condpairs[i,1],'\n',condpairs[i,2]))
      }
    }
    if (!is.null(tks)) tks=ggbio::tracks(tks)
  }

  if(is.null(tks)) cat("No plot, perhaps the filtered PACds is empty")
  return(tks)
}
)



# Given two condpair matrix, find matched queryCp
# If markFlip=TRUE, then return 3 columns ,the last column is T/F, with TRUE denoting the queryCp is flipped of subjectCp.
# @param queryCp: matrix/df, or a vector (cond1, cond2)
# @param subjectCp: matrix/df
# @param return: df[cond1, cond2, <flip>]
# @examples
# queryCp=matrix(c('a','b','c','d','x','y'),nrow=3, byrow = TRUE)
# queryCp=matrix(c('a1','b1','c1','d1','x1','y1'),nrow=3, byrow = TRUE)
# subjectCp=matrix(c('1','2','b','a','c','d','x','y'),nrow=4, byrow = TRUE)
# .getMatchCondPairs(queryCp, subjectCp, TRUE)
# .getMatchCondPairs(queryCp, subjectCp, F)
.getMatchCondPairs <- function(queryCp, subjectCp, markFlip=FALSE) {
  if (is.null(queryCp)) return(NULL)
  if (is.vector(queryCp)) {
    if (length(queryCp)!=2) {
      cat("Error queryCp, is vector but not length=2\n")
      return(NULL)
    }
    queryCp=matrix(queryCp, nrow=1)
  }
  if ((!is.matrix(queryCp) & !is.data.frame(queryCp))) {
    cat("Error queryCp, must be matrix or data.frame\n")
    return(NULL)
  }
  if (!(ncol(queryCp)==2)) {
    cat("Error queryCp, must be 2 cols\n")
    return(NULL)
  }
  if (!(ncol(subjectCp)==2)) stop("Error subjectCp, must be 2 cols")
  cp=paste0(subjectCp[,1],'.',subjectCp[,2])
  cp1=paste0(queryCp[,1],'.',queryCp[,2])
  cp2=paste0(queryCp[,2],'.',queryCp[,1])
  id1=which(cp %in% cp1)
  id2=which(cp %in% cp2)
  id=c(id1, id2)
  if (length(id)==0) {
    cat("No one line of queryCp in subjectCp, return NULL!\n")
    return(NULL)
  }
  cp1=subjectCp[id1, , drop=F]
  cp2=subjectCp[id2, , drop=F]
  cp=rbind(cp1,cp2)
  colnames(cp)=c('cond1', 'cond2')
  rownames(cp)=paste0(cp[,1],'.',cp[,2])
  if (markFlip) {
    cp=cbind(cp, flip=c(rep(FALSE, nrow(cp1)), rep(TRUE, nrow(cp2))))
  }

  cp =.dropColFactor(cp)
  return(.asDf(cp))
}

# drop factor to character column for factor columns in matrix or df
.dropColFactor<-function(d) {
  if (class(d) %in% c('matrix','data.frame')) {
    for (i in 1:ncol(d)) {
      if (is.factor(d[,i])) {
        d[,i]=as.character(d[,i])
      }
    }
  }
  return(d)
}

## Preprocess for all movViz functions
.movVizPreprocess <- function(object, gene, PACds,
                              conds, avg, pool, clearPAT,
                              highlightConds, collapseConds) {

  if (!is.null(highlightConds)) {
    if(!collapseConds) stop("highlightConds<>NULL but collapseConds=FALSE, cannot highlight track")
    if (!(ncol(highlightConds)==2)) stop("Error highlightConds in movViz, must be 2 cols")
    highlightConds=.getMatchCondPairs(queryCp=highlightConds, subjectCp =object@conds, markFlip=FALSE )
  }

  noIgt=FALSE
  if (class(object) %in% c('movDEGeneRes', 'movDEPACRes')) noIgt=TRUE
  PACds=subsetPACds(PACds, genes=gene, group=object@group, conds=NULL, noIntergenic = noIgt, clearPAT=clearPAT, avg=avg, pool=pool)

  #get matched conds in object@conds
  #flip if order is flipped.
  condpairs=conds
  if (is.null(conds)) {
    condpairs=object@conds
  } else {
    condpairs=.getMatchCondPairs(queryCp=condpairs, subjectCp =object@conds, markFlip=FALSE )
  }

  condpairs=.dropColFactor(condpairs)

  # highlightConds: the condpairs to highlight, are all in object@conds
  # condpairs: condpairs to plot, all in object@conds
  # PACds: filtered PACds by gene and other parameters
  return(list(highlightConds=highlightConds, condpairs=condpairs, PACds=PACds ))
}


#' Visualize a movDEGeneRes object
#'
#' movViz(movDEGeneRes) plots PAC distributions of a gene across selected samples from a movDEGeneRes object.
#'
#' @param padjThd If a PAC's padj<padjThd, then highlight the PACs. Cannot be used with highlightConds=TRUE.
#' @param showPV show gene's pvalue in the left-top of the track.
#' Only shown when collapse=TRUE, and it is irrelevent with padjThd.
#' @param PACds A PACdataset.
#' @examples
#'
#' ## Visualize DE gene
#' ## First get DE gene results.
#' DEgene=movDEGene(PACds=PACds, method='DESeq2', group='group', minSumPAT=10, norm=TRUE)
#' gene=rownames(DEgene@pairwise$padj[order(DEgene@pairwise$padj[,1]), ,drop=F])[1]
#' DEgene@pairwise$padj[gene, ]
#'
#' ## Show a DE gene without gene model.
#' movViz(object=DEgene, gene=gene, txdb=NULL, PACds=PACds)
#' ## Show a DE gene with the gene model.
#' movViz(object=DEgene, gene=gene, txdb=gff, PACds=PACds)
#'
#'## Show each condition pair in one track and show the pvalues.
#'movViz(object=DEgene, gene=gene, txdb=NULL, PACds=PACds, collapseConds=TRUE)
#'movViz(object=DEgene, gene=gene, txdb=NULL, PACds=PACds, collapseConds=TRUE, showPV=FALSE)
#'
#' ## Highlight specific conditions in blue.
#' movViz(object=DEgene, gene=gene, txdb=NULL, PACds=PACds,
#'      collapseConds=TRUE, conds=DEgene@conds[1:3, ], highlightConds=DEgene@conds[2, ])
#' ## Highlight conditions with padj<0.01.
#' movViz(object=DEgene, gene=gene, txdb=NULL, PACds=PACds, collapseConds=TRUE, padjThd=0.01)
#'
#' @describeIn movViz for a movRes object.
#' @export
setMethod("movViz", signature="movDEGeneRes", def=function(object, gene, txdb=NULL, PACds,
                                                           conds=NULL, avg=FALSE, pool=FALSE, clearPAT=0,
                                                           padjThd=NULL, highlightConds=NULL, collapseConds=FALSE,
                                                           showPV=TRUE) {

  #library(ggbio, verbose = FALSE)

  pls=.movVizPreprocess(object=object, gene=gene, PACds=PACds,
                        conds=conds, avg=avg, pool=pool, clearPAT=clearPAT,
                        highlightConds=highlightConds, collapseConds=collapseConds)

  highlightConds=pls[['highlightConds']]
  PACds=pls[['PACds']]
  condpairs=pls[['condpairs']]

  if (!is.null(highlightConds)) highlightConds=paste0(highlightConds[,1], '.', highlightConds[,2])
  #If padjThd, then set pv<hpv conds as hightlight conds
  pvs=object@pairwise$padj[gene, , drop=FALSE]
  if (!is.null(padjThd)) {
    if (!is.null(highlightConds)) stop("padjThd or highlightConds, cannot both")
    if(!collapseConds) stop("padjThd<>NULL but collapseConds=FALSE, cannot highlight track")
    hidx=which(!is.na(pvs) & pvs<padjThd)
    if (length(hidx)>0) {
      highlightConds=colnames(object@pairwise$padj)[hidx]
    }
  }

  if (!collapseConds) {
    tks=movViz(object=PACds, gene=gene, txdb=txdb, group=object@group, conds=unique(unlist(condpairs)),
               collapseConds=FALSE, linkPAs=FALSE)
  } else {
    for (i in 1:nrow(condpairs)) {
      #get condpair tracks
      tks1=movViz(object=PACds, gene=gene, txdb=txdb, group=object@group, conds=unlist(condpairs[i, ]),
                  collapseConds=TRUE, linkPAs=FALSE, otrack=FALSE)

      #highlight a track (user given or meet pvalue)
      if(!is.null(highlightConds)) {
        if ((paste0(condpairs[i,1],'.',condpairs[i,2]) %in% highlightConds)) {
          tks1[[2]]=tks1[[2]]+.themeHighlightConds()
        }
      }

      #show pvalue text
      if(showPV) {
        condname=getMovResPairwiseCondName(object,c(condpairs[i,1],condpairs[i,2]),markFlip=FALSE, check=TRUE)
        lbl=sprintf('P=%.4f\nVal=%.4f', pvs[1, condname], object@pairwise$value[gene, condname])
        tks1[[2]]=tks1[[2]]+annotate("text", -Inf, Inf,  hjust = 0, vjust = 1, label = lbl, colour='black', size=3)
      }

      if (i==1) {
        tks=tks1
      } else {
        tks=c(tks, tks1[2])
      }
    }
    names(tks)=c('Gene\nModel', paste0(condpairs[,1],'\n',condpairs[,2]))
    tks=ggbio::tracks(tks)
  }

  return(tks)
}
)

# ---------- movViz movDEPACRes ----------

#' Visualize a movDEPACRes object.
#'
#' movViz(movDEPACRes) plots PAC distributions of a gene across selected samples from a movDEPACRes object.
#' It uses additional options: [padjThd, showAllPA, showRatio].
#'
#' movViz(movDEPACRes) is Similar to movViz(DEGene), but can control PACs.
#' Use padjThd to filter DEPACs and highlight DEPACs.
#' All DEPACs's pvalue and logFC of a condpair are shown on the track.
#'
#' @param padjThd for object of movDEPACRes, if a PAC's padj<padjThd then highlight or only show (if showAllPA=FALSE) the PACs.
#' @param showAllPA If TRUE, then show all PACs including non-DE ones. Otherwise only show DEPACs with padj<padjThd.
#' @describeIn movViz for a movDEPACRes object.
#' @export
setMethod("movViz", signature="movDEPACRes", def=function(object, gene, txdb=NULL, PACds,
                                                          conds=NULL, avg=FALSE, pool=FALSE, clearPAT=0,
                                                          padjThd=NULL, highlightConds=NULL, collapseConds=FALSE,
                                                          showPV=TRUE, showAllPA=TRUE, showRatio=FALSE) {
  #library(ggbio, verbose = FALSE)
  pls=.movVizPreprocess(object=object, gene=gene, PACds=PACds,
                        conds=conds, avg=avg, pool=pool, clearPAT=clearPAT,
                        highlightConds=highlightConds, collapseConds=collapseConds)

  highlightConds=pls[['highlightConds']]
  if (!is.null(highlightConds)) highlightConds=paste0(highlightConds[,1], '.', highlightConds[,2])
  PACds=pls[['PACds']]
  condpairs=pls[['condpairs']]

  #filter DEPAC's gene PAC
  gidx=grep(gene, rownames(object@pairwise$padj))
  if (length(gidx)==0) {
    ## not gene:coord, then filter PAs from PACds
    gidx=rownames(PACds@anno[PACds@anno$gene==gene, ])
  }
  object@pairwise$padj=object@pairwise$padj[gidx, , drop=F]
  object@pairwise$value=object@pairwise$value[gidx, , drop=F]
  if (nrow(object@pairwise$padj)==0) {
    cat("Warning: no lines in object is gene=",gene,', will not highlight PAs\n')
  }

  highlightPAs=NULL

  #get DEPA for each condpair
  if (!is.null(padjThd)) {
    if (nrow(object@pairwise$padj)>0) {
    paid=rownames(object@pairwise$padj)
    highlightPAs=lapply(object@pairwise$padj, function(par) return( paid[which(!is.na(par) & par<padjThd)]))
    }
  }

  if (!collapseConds) {

    if (showAllPA) { #show all PA
      tks=movViz(object=PACds, gene=gene, txdb=txdb, group=object@group, conds=unique(unlist(condpairs)),
                 collapseConds=FALSE, linkPAs=FALSE, showRatio=showRatio, highlightPAs=unique(unlist(highlightPAs)))
    } else { #show only DEPA
      PACds=subsetPACds(PACds, PAs=unique(unlist(highlightPAs)))
      tks=movViz(object=PACds, gene=gene, txdb=txdb, group=object@group, conds=unique(unlist(condpairs)),
                 collapseConds=FALSE, linkPAs=FALSE, showRatio=showRatio, highlightPAs=unique(unlist(highlightPAs)))
    }

  } else {

    tks=NULL
    for (i in 1:nrow(condpairs)) {
      #get DEPA names of a condpair
      hlPAs=highlightPAs[[paste0(condpairs[i,1],'.',condpairs[i,2])]]

      PACds1=PACds
      if (!showAllPA) { #show DEPA
        PACds1=subsetPACds(PACds, PAs=hlPAs)
      }
      #print(PACds1@counts)

      tks1=movViz(object=PACds1, gene=gene, txdb=txdb, group=object@group, conds=unlist(condpairs[i, ]),
                  collapseConds=TRUE, linkPAs=FALSE, otrack=FALSE, highlightPAs=hlPAs, showRatio=showRatio)


      if (is.null(tks1)) next

      #highlight a user given condpairs
      if(!is.null(highlightConds)) {
        if ((paste0(condpairs[i,1],'.',condpairs[i,2]) %in% highlightConds)) {
          tks1[[2]]=tks1[[2]]+.themeHighlightConds()
        }
      }

      ## ??? debug
      #show every DEPA's pvalue and logFC. Each text line for one DEPAC.
      if(showPV & !is.null(hlPAs)) {
        condname=getMovResPairwiseCondName(object, c(condpairs[i,1],condpairs[i,2]), markFlip=FALSE, check=TRUE)
        pvs=cbind(PA=hlPAs, P=object@pairwise$padj[hlPAs, condname], logFC=object@pairwise$value[hlPAs, condname])
        pvs=pvs[order(pvs[,'PA']), 2:3, drop=F]
        lbl=apply(pvs, 1, function(par) {par=as.numeric(par); return(sprintf('P=%s Val=%.2f', format.pval(par[1]), par[2]))})
        lbl=paste(lbl, collapse='\n')
        #print(lbl)
        tks1[[2]]=tks1[[2]]+annotate("text", -Inf, Inf,  hjust = 0, vjust = 1, label = lbl, colour=hightlightPAcolor, size=3)
      }

      if (is.null(tks)) {
        tks=tks1
        names(tks)=c('Gene\nModel', paste0(condpairs[i,1],'\n',condpairs[i,2]))
      } else {
        n1=names(tks)
        tks=c(tks, tks1[2])
        names(tks)=c(n1, paste0(condpairs[i,1],'\n',condpairs[i,2]))
      }
    }
    if (!is.null(tks)) tks=ggbio::tracks(tks)
  }

  if(is.null(tks)) cat("No plot, perhaps the filtered PACds is empty, try to change options like showAllPA=TRUE\n")
  return(tks)
}
)


# ---------- movViz movUTRTrendRes ----------

#' Visualize a movUTRTrendRes object
#'
#' movViz(movUTRTrendRes) plots PAC distributions of a gene across selected samples from a movUTRTrendRes object.
#' The gene is not neccesary a significant gene, but just plot the gene.
#'
#' @describeIn movViz for a movUTRTrendRes object.
#' @export
setMethod("movViz", signature="movUTRTrendRes", def=function(object, gene, txdb=NULL, PACds,
                                                             conds=NULL, avg=FALSE, pool=FALSE, clearPAT=0,
                                                             padjThd=NULL, highlightConds=NULL, collapseConds=FALSE,
                                                             showPV=TRUE, showAllPA=TRUE, showRatio=FALSE, linkPAs=FALSE) {
  #library(ggbio, verbose = FALSE)
  pls=.movVizPreprocess(object=object, gene=gene, PACds=PACds,
                        conds=conds, avg=avg, pool=pool, clearPAT=clearPAT,
                        highlightConds=highlightConds, collapseConds=collapseConds)

  highlightConds=pls[['highlightConds']]
  if (!is.null(highlightConds)) highlightConds=paste0(highlightConds[,1], '.', highlightConds[,2])
  PACds=pls[['PACds']]
  condpairs=pls[['condpairs']]

  #get swPA for each condpair
  highlightPAs=list()
  if (!is.null(padjThd)) {
    for (i in 1:nrow(condpairs)) {
      highlightPAs[[paste0(condpairs[i,1],'.',condpairs[i,2])]]=.getPAbyGeneFromMovUTRTrendRes(object, genes=gene, condname=c(condpairs[i,1],condpairs[i,2]))
    }
  }

  if (!collapseConds) {

    if (showAllPA) {
      tks=movViz(object=PACds, gene=gene, txdb=txdb, group=object@group, conds=unique(unlist(condpairs)),
                 collapseConds=FALSE, linkPAs=FALSE, showRatio=showRatio, highlightPAs=unique(unlist(highlightPAs)))
    } else {
      PACds=subsetPACds(PACds, PAs=unique(unlist(highlightPAs)))
      tks=movViz(object=PACds, gene=gene, txdb=txdb, group=object@group, conds=unique(unlist(condpairs)),
                 collapseConds=FALSE, linkPAs=FALSE, showRatio=showRatio, highlightPAs=unique(unlist(highlightPAs)))
    }

  } else {
    tks=NULL
    for (i in 1:nrow(condpairs)) {
      hlPAs=highlightPAs[[paste0(condpairs[i,1],'.',condpairs[i,2])]]

      PACds1=PACds
      if (!showAllPA) {
        PACds1=subsetPACds(PACds, PAs=hlPAs)
      }

      tks1=movViz(object=PACds1, gene=gene, txdb=txdb, group=object@group, conds=unlist(condpairs[i, ]),
                  collapseConds=TRUE, linkPAs=linkPAs, otrack=FALSE, highlightPAs=hlPAs, showRatio=showRatio)

      if(is.null(tks1)) next

      if(!is.null(highlightConds)) {
        if ((paste0(condpairs[i,1],'.',condpairs[i,2]) %in% highlightConds)) {
          tks1[[2]]=tks1[[2]]+.themeHighlightConds()
        }
      }

      if(showPV & !is.null(hlPAs)) {
        condname=getMovResPairwiseCondName(object,c(condpairs[i,1],condpairs[i,2]),markFlip=FALSE, check=TRUE)
        lbl=sprintf('P=%.4f\nVal=%.4f', object@pairwise$padj[gene, condname], object@pairwise$value[gene, condname])
        tks1[[2]]=tks1[[2]]+annotate("text", -Inf, Inf,  hjust = 0, vjust = 1, label = lbl, colour=hightlightPAcolor, size=3)
      }

      if (is.null(tks)) {
        tks=tks1
        names(tks)=c('Gene\nModel', paste0(condpairs[i,1],'\n',condpairs[i,2]))
      } else {
        n1=names(tks)
        tks=c(tks, tks1[2])
        names(tks)=c(n1, paste0(condpairs[i,1],'\n',condpairs[i,2]))
      }
    }
    if (!is.null(tks)) tks=ggbio::tracks(tks)
  }

  if(is.null(tks)) cat("No plot, perhaps the filtered PACds is empty, try to change options like showAllPA=TRUE")
  return(tks)
}
)


#' Visualize a movAPASwitchRes object
#'
#' movViz(movAPASwitchRes) plots PAC distributions of a gene across selected samples from a movAPASwitchRes object. Same as movViz(movUTRTrendRes).
#'
#' @describeIn movViz for movAPASwitchRes and movUTRTrendRes object.
#' @export
setMethod("movViz", signature(object = "movAPASwitchRes"), function(object, gene, txdb=NULL, PACds,
                                                                    conds=NULL, avg=FALSE, pool=FALSE, clearPAT=0,
                                                                    padjThd=NULL, highlightConds=NULL, collapseConds=FALSE,
                                                                    showPV=TRUE, showAllPA=TRUE, showRatio=FALSE, linkPAs=TRUE)  {
  object=new("movUTRTrendRes", group=object@group, method=object@method, conds=object@conds, pairwise=object@pairwise, fullList=object@fullList)
  return(movViz(object=object, gene=gene, txdb=txdb, PACds=PACds,
                conds=conds, avg=avg, pool=pool, clearPAT=clearPAT,
                padjThd=padjThd, highlightConds=highlightConds, collapseConds=collapseConds,
                showPV=showPV, showAllPA=showAllPA, showRatio=showRatio, linkPAs=linkPAs))
}
)




