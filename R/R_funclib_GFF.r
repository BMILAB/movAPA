# 2019-10-29 funclibs for parsing genome annotations

#' Parse genome annotation
#'
#' parseGenomeAnnotation parses different types of genome annotations.
#'
#' Due to the complex GFF3/GTF/TxDB structure of different genome annotation files from different species,
#' this function may not be always applicable for any given file. You may need to check mannually.
#' @usage parseGenomeAnnotation(aGFF)
#' @param anAnno can be a list of anno.rna/anno.need, or .rda/.rdata/.gff3/.gtf file name, or TxDB object.
#' @return a parsed genome annotation object, which is a list of three elements (anno.rna, anno.need, anno.frame) and can be used for annotatePAC().
#' @examples
#' ## Way1: Based on an annotation file in gff3 format, You can dowonload annotation from Ensemble Plants
#' #Prepare the annotation
#' #wget -c ftp://ftp.ensemblgenomes.org/pub/plants/release-44/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.44.gff3.gz
#' #gunzip  Arabidopsis_thaliana.TAIR10.44.gff3.gz
#' gff.path <- "/path/Arabidopsis_thaliana.TAIR10.44.gff3"
#' anno <- parseGenomeAnnotation(anAnno=gff.path)
#'
#' ##way2: load from a .rda file (already processed file)
#' anno <- parseGenomeAnnotation("anno.rda")
#'
#' ##Way3: Based on a TxDb object generated from BioMart.
#' # Parse Arabidopsis Txdb
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' anno <- parseGenomeAnnotation(TxDb.Athaliana.BioMart.plantsmart28)
#' # Parse mm10 Txdb
#' BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")
#' library(TxDb.Mmusculus.UCSC.mm10.ensGene)
#' anno <- parseGenomeAnnotation(TxDb.Mmusculus.UCSC.mm10.ensGene)
#' @name parseGenomeAnnotation
#' @seealso [annotatePAC()] to annotate a PACdataset.
#' @family genome annotation functions
#' @export
parseGenomeAnnotation <- function(anAnno) {
  #library(rtracklayer)
  #library(GenomicRanges)
  #library(GenomicFeatures)

  if (class(anAnno)=='list') {
    if ( sum(names(anAnno) %in% c('anno.rna', 'anno.need', 'anno.frame'))!=3) stop("anAnno is a list, but no anno.rna/anno.need/anno.frame!")
    return(anAnno)
  }

  if (is.character(anAnno)) {
    if (grepl('\\.rda|\\.rdata', tolower(anAnno))) {
      if (!file.exists(anAnno)) {
        stop("anAnno is .rda/.rdata but file not exists!")
      }
      a = new.env()
      load(anAnno, envir = a)
      for (v in ls(a)) {
        if (class(get(v, envir = a))=='list') {
          if (!(AinB(c('anno.rna','anno.need','anno.frame'), names(get(v, envir = a))))) next
        } else {
          next
        }
        return(get(v, envir = a))
      }
      stop('No list(anno.rna, anno.need, anno.frame) in .rda file anAnno')

    } else if (grepl('\\.gff3|\\.gtf', tolower(anAnno))) {
      rt=parseGff(anAnno)

    }
    invisible(gc())
    return(rt)
  }#~chr

  if (class(anAnno)=='TxDb') {
    rt=parseTxdb(anAnno)
    invisible(gc())
    return(rt)
  }

}


#' Parse TxDb genome annotation
#'
#' parseTxdb parses genome annotation object of TxDb
#'
#' @usage parseTxdb(aGFF)
#' @param an.txdb a TxDb object
#' @return a parsed genome annotation object, which is a list of three elements (anno.rna, anno.need, anno.frame) and can be used for annotatePAC().
#' @examples
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' txdbAnno <- parseTxdb(an.txdb=TxDb.Athaliana.BioMart.plantsmart28)
#' @name parseTxdb
#' @seealso [parseGff()] to parse a Gff file.
#' @family Genome annotation functions
#' @export
parseTxdb <- function (an.txdb) {

  if(class(an.txdb)!='TxDb') stop("an.txdb not of class TxDb!")

  genes <- genes(an.txdb,columns=c("tx_type","gene_id"))
  genes <- as.data.frame(genes)

  genes <- data.frame(seqnames=as.character(genes$seqnames) ,start=as.integer(genes$start),
                      end=as.integer(genes$end),width=as.integer(genes$width),
                      strand=as.character(genes$strand),type="gene",
                      ID =as.character(genes$gene_id),biotype=as.character(genes$tx_type),
                      gene_id =as.character(genes$gene_id),Parent=NA,transcript_id=NA)
  #setdiff(colnames(genes),colnames(tari))
  rnas <- transcripts(an.txdb,columns=c("tx_name","tx_type","gene_id"))
  rnas<- as.data.frame(rnas)

  #test <- strsplit(as.character(rnas$gene_id) ,"\\s+")
  #temp3 <- paste("",lapply(test,"[[",1),sep="");
  #head(temp3)
  rnas <- data.frame(seqnames=as.character(rnas$seqnames) ,start=as.integer(rnas$start),
                     end=as.integer(rnas$end),width=as.integer(rnas$width),
                     strand=as.character(rnas$strand),type="RNA",
                     ID =as.character(rnas$tx_name),biotype=as.character(rnas$tx_type),
                     gene_id =as.character(rnas$gene_id),Parent=as.character(rnas$gene_id),
                     transcript_id=as.character(rnas$tx_name))
  # exons <- exons(an.txdb,columns=c("exon_name","tx_name","tx_type","gene_id"))
  # exons <- as.data.frame(exons)
  # head(exons)

  exons <- exonsBy(an.txdb,by=c("tx"),use.names=TRUE)
  exons <- as.data.frame(exons)
  exons <- data.frame(seqnames=as.character(exons$seqnames) ,start=as.integer(exons$start),
                      end=as.integer(exons$end),width=as.integer(exons$width),
                      strand=as.character(exons$strand),type="exon",
                      ID =as.character(exons$exon_name),biotype=NA,
                      gene_id =NA,Parent=as.character(exons$group_name),
                      transcript_id=as.character(exons$group_name))
  index <- match(exons$Parent,rnas$transcript_id)
  #which(is.na(index))
  exons$gene_id <- rnas$Parent[index]
  exons$biotype <- rnas$biotype[index]

  #==================================
  #CDS
  cdss <- cdsBy(an.txdb,by=c("tx"),use.names=TRUE)
  cdss <- as.data.frame(cdss)
  cdss <- data.frame(seqnames=as.character(cdss$seqnames) ,start=as.integer(cdss$start),
                     end=as.integer(cdss$end),width=as.integer(cdss$width),
                     strand=as.character(cdss$strand),type="CDS",
                     ID =as.character(cdss$cds_name),biotype=NA,
                     gene_id =NA,Parent=as.character(cdss$group_name),
                     transcript_id=as.character(cdss$group_name))
  index <- match(cdss$Parent,rnas$transcript_id)
  #which(is.na(index))
  cdss$gene_id <- rnas$Parent[index]
  cdss$biotype <- rnas$biotype[index]
  #head(cdss)
  #cdss <- cds(an.txdb,columns=c("cds_name","tx_name","tx_type","gene_id"))

  #==================================
  #introns
  introns <- intronsByTranscript(an.txdb,use.names=TRUE)
  introns <- as.data.frame(introns)
  introns <- data.frame(seqnames=as.character(introns$seqnames) ,start=as.integer(introns$start),
                        end=as.integer(introns$end),width=as.integer(introns$width),
                        strand=as.character(introns$strand),type="intron",
                        ID =NA,biotype=NA,
                        gene_id =NA,Parent=as.character(introns$group_name),
                        transcript_id=as.character(introns$group_name))
  index <- match(introns$Parent,rnas$transcript_id)
  #which(is.na(index))
  introns$gene_id <- rnas$Parent[index]
  introns$biotype <- rnas$biotype[index]
  #head(introns)

  #===================================================
  #five UTR
  fiveUTRs <- fiveUTRsByTranscript(an.txdb,use.names=TRUE)
  fiveUTRs <- as.data.frame(fiveUTRs)
  fiveUTRs <- data.frame(seqnames=as.character(fiveUTRs$seqnames) ,start=as.integer(fiveUTRs$start),
                         end=as.integer(fiveUTRs$end),width=as.integer(fiveUTRs$width),
                         strand=as.character(fiveUTRs$strand),type="five_prime_UTR",
                         ID =NA,biotype=NA,
                         gene_id =NA,Parent=as.character(fiveUTRs$group_name),
                         transcript_id=as.character(fiveUTRs$group_name))
  index <- match(fiveUTRs$Parent,rnas$transcript_id)
  #which(is.na(index))
  fiveUTRs$gene_id <- rnas$Parent[index]
  fiveUTRs$biotype <- rnas$biotype[index]
  #head(fiveUTRs)


  #===========================================
  #three UTR
  threeUTRs <- threeUTRsByTranscript(an.txdb,use.names=TRUE)
  threeUTRs <- as.data.frame(threeUTRs)
  threeUTRs <- data.frame(seqnames=as.character(threeUTRs$seqnames) ,start=as.integer(threeUTRs$start),
                          end=as.integer(threeUTRs$end),width=as.integer(threeUTRs$width),
                          strand=as.character(threeUTRs$strand),type="three_prime_UTR",
                          ID =NA,biotype=NA,
                          gene_id =NA,Parent=as.character(threeUTRs$group_name),
                          transcript_id=as.character(threeUTRs$group_name))
  index <- match(threeUTRs$Parent,rnas$transcript_id)
  #which(is.na(index))
  threeUTRs$gene_id <- rnas$Parent[index]
  threeUTRs$biotype <- rnas$biotype[index]

  anno.frame <- rbind(genes,rnas,exons,cdss,introns,fiveUTRs,threeUTRs)
  anno.frame$type <- factor(anno.frame$type,levels=c("gene","RNA","five_prime_UTR","exon","CDS","intron",
                                                     "three_prime_UTR"))
  #anno.frame <- anno.frame[order(anno.frame$transcript_id,anno.frame$gene_id,
  #                               anno.frame$start,anno.frame$strand,anno.frame$type),]
  anno.need <- rbind(exons,cdss,introns,fiveUTRs,threeUTRs)
  anno.rna <- rnas

  return(list(anno.need=anno.need, anno.rna=anno.rna, anno.frame=anno.frame))
}


#' Parse gff3/gtf genome annotation
#'
#' parseGff parses genome annotation file of gff3/gtf format
#'
#' Due to the complex GFF3/GFF/GTF structure of different genome annotation files from different species,
#' this function may not be always applicable for any given file. You may need to check mannually.
#' @usage parseGff(aGFF)
#' @param aGFF .gff3/.gff/.gtf file name
#' @return a parsed genome annotation object, which is a list of three elements (anno.rna, anno.need, anno.frame) and can be used for annotatePAC().
#' @examples
#' ## parse from a gff file, and save as .rda for further use.
#' gff=parseGff(aGFF='Bamboo.Hic.gff')
#' @name parseGff
#' @seealso [parseTxdb()] to parse a Txdb object.
#' @family genome annotation functions
#' @export
parseGff <- function(aGFF) {

  if (!is.character(aGFF)) stop("aGFF not a character string!")
  if (!grepl('\\.gff3|\\.gtf', tolower(aGFF))) stop('aGFF not .gff3/.gff/.gtf!')

  if (grepl('\\.gff3|\\.gff', tolower(aGFF))) {
      #------------------------------------------------------
      #Loading annotation (gff3 format)
      #-------------------------------------------------------
      gff.path=aGFF
      anno <- import.gff3(gff.path)
      anno.frame <- as.data.frame(anno,stringsAsFactors =FALSE)
      anno.frame$seqnames <- as.character(anno.frame$seqnames)
      anno.frame$strand <- as.character(anno.frame$strand)
      anno.frame$type <- as.character(anno.frame$type)
      #print("###annotation file type information")
      #print(table(anno.frame$type))
      #delete chromosome information
      anno.frame$Parent <- sub(pattern="\\S+\\:",replacement = "",anno.frame$Parent)
      anno.frame$ID <- sub(pattern="\\S+\\:",replacement = "",anno.frame$ID)
      if(length(which(anno.frame$type=="chromosome"))){
        anno.frame <- anno.frame[-which(anno.frame$type=="chromosome"),]
      }
      #instead transcript to RNA
      anno.frame$type[which(anno.frame$type == "transcript")] <-"RNA"
      #getting RNA row
      rna.id <- grep("RNA$",anno.frame$type,ignore.case = FALSE)
      anno.rna <- anno.frame[rna.id,]

    } else if (grepl('\\.gtf', tolower(aGFF))) {
      gtf.path=aGFF
      anno <- import(gtf.path,format="gtf")
      anno.frame <- as.data.frame(anno,stringsAsFactors =FALSE)
      anno.frame$seqnames <- as.character(anno.frame$seqnames)
      anno.frame$strand <- as.character(anno.frame$strand)
      anno.frame$type <- as.character(anno.frame$type)
      anno.frame$Parent <- as.character(anno.frame$transcript_id)
      anno.frame$type[which(anno.frame$type == "transcript")] <-"RNA"
      #getting RNA row
      trans.id <- grep("transcript",anno.frame$type,ignore.case = FALSE)
      if(length(trans.id)){
        anno.frame$type[which(anno.frame$type == "transcript")] <-"RNA"
        rna.id <- grep("RNA$",anno.frame$type,ignore.case = FALSE)
      }else{
        rna.id <- grep("RNA$",anno.frame$type,ignore.case = FALSE)
      }
      if(length(rna.id)==0){
        anno.frame <- add_rna(anno.frame)
        rna.id <- grep("RNA$",anno.frame$type,ignore.case = FALSE)
      }
      if(length(which(anno.frame$type=="gene"))==0){
        anno.gene <- anno.frame[rna.id,]
        anno.gene$type<- "gene"
        anno.gene$Parent <- ""
        anno.frame <- rbind(anno.frame,anno.gene)
      }
      anno.frame$ID <- anno.frame$Parent
      if(length(which(anno.frame$type=="chromosome"))){
        anno.frame <- anno.frame[-which(anno.frame$type=="chromosome"),]
      }
      anno.rna <- anno.frame[rna.id,]
    } #~gtf


    #If the comment is incomplete, losing transcript_id
    if(!length(which(colnames(anno.rna) == "transcript_id"))){
      anno.rna$transcript_id<-anno.rna$ID
    }
    #table(anno.rna$type)

    #ID=transcript:AT1G01010.1;Parent=gene:AT1G01010;biotype=protein_coding;transcript_id=AT1G01010.1
    #1	araport11	five_prime_UTR	3631	3759	.	+	.	Parent=transcript:AT1G01010.1
    #confirm that the names of transcript is consistent with exon/cds/utr
    if(length(setdiff(anno.rna$transcript_id,anno.frame$Parent))){
      stop("Not consistent between transcript id in rna and exon/cds/utr")
    }
    #anno.frame$Parent[which(anno.frame$type=="gene")]<- ""
    # anno.frame.raw <- anno.frame
    if(is.na(match("three_prime_UTR",unique(anno.frame$type)))){
      if(is.na(match("CDS",unique(anno.frame$type)))){
        warning("This annotation without CDS, we can't identify UTR region")
      }else{
        print("Extracting UTR region")
        anno.frame <- add_utr(anno.frame)
      }
    }
    #=========================================================================
    #anno.need store cds/exon/5utr/3utr information
    anno.need <- anno.frame[which(anno.frame$Parent %in% anno.rna$transcript_id),]

    need.rna.id <- grep("RNA$",anno.need$type,ignore.case = FALSE)
    if(length(need.rna.id)){
      anno.need<-anno.need[-need.rna.id,]
    }

    index <- match(anno.need$Parent,anno.rna$transcript_id)

    if(length(which(is.na(index)))){
      stop("error can't find exon/cds/5utr/3utr 's parent")
    }
    anno.need$gene_id <- anno.rna$Parent[index]

    if(is.na(match("biotype",colnames(anno.rna)))){
      anno.rna$biotype <- NA
    }
    anno.need$biotype <- anno.rna$biotype[index]
    #====================================================================
    #ann.intron stores intron information
    exon.id <- grep("exon",anno.need$type,ignore.case = FALSE)
    ann.exon <- anno.need[exon.id,]
    if(length(which(is.na(ann.exon$Parent)))){
      print("exist some exon can't find parent id ")
    }
    ann.exon <- ann.exon[order(ann.exon$Parent,ann.exon$start,ann.exon$strand),]
    ann.exon.1 <- ann.exon[seq(1,nrow(ann.exon),2),]
    ann.exon.2 <- ann.exon[seq(2,nrow(ann.exon),2),]
    ann.exon.3 <- ann.exon[seq(3,nrow(ann.exon),2),]

    keep.num1 <- min(nrow(ann.exon.1),nrow(ann.exon.2))
    ann.exon.k1<-ann.exon.1[1:keep.num1,]
    ann.exon.k2<-ann.exon.2[1:keep.num1,]
    index <- which(ann.exon.k1$Parent == ann.exon.k2$Parent)
    if(!identical(ann.exon.k1$Parent[index],ann.exon.k2$Parent[index])){
      stop("something error with extart intron region")
    }
    ann.intron1 <- ann.exon.k1[index,]
    ann.intron1$type <- "intron"
    ann.intron1$start <- ann.exon.k1$end[index]+1
    ann.intron1$end <- ann.exon.k2$start[index]-1


    keep.num2 <- min(nrow(ann.exon.2),nrow(ann.exon.3))
    ann.exon.kk2<-ann.exon.2[1:keep.num2,]
    ann.exon.k3<-ann.exon.3[1:keep.num2,]
    index <- which(ann.exon.kk2$Parent == ann.exon.k3$Parent)
    if(!identical(ann.exon.kk2$Parent[index],ann.exon.k3$Parent[index])){
      stop("something error with extart intron region")
    }
    ann.intron2 <- ann.exon.kk2[index,]
    ann.intron2$type <- "intron"
    ann.intron2$start <- ann.exon.kk2$end[index]+1
    ann.intron2$end <- ann.exon.k3$start[index]-1
    ann.intron <- rbind(ann.intron1,ann.intron2)
    ann.intron <- ann.intron[order(ann.intron$Parent,ann.intron$start,ann.intron$strand),]
    anno.need <- rbind(anno.need,ann.intron)


    #table(anno.need$type)
    rna.error <- grep("RNA$",anno.need$type,ignore.case = FALSE)
    if(length(rna.error)){
      anno.need <- anno.need[-rna.error,]
    }
    return(list(anno.need=anno.need, anno.rna=anno.rna, anno.frame=anno.frame))
}



#=========================================================
#------------------------------------------------------
#function:add_utr()
#Adding 3UTR and 5UTR region
#--------------------------------------------------------
#======================================================
add_utr <- function(anno.frame=NULL){
  anno.cds <- anno.frame[which(anno.frame$type=="CDS"),]
  anno.exon <- anno.frame[which(anno.frame$type=="exon"),]
  rna.id <- grep("RNA$",anno.frame$type,ignore.case = FALSE)
  anno.rna <- anno.frame[rna.id,]
  if(!length(which(colnames(anno.rna) == "transcript_id"))){
    anno.rna$transcript_id<-anno.rna$ID
  }
  anno.cds.frist <- anno.cds[order(anno.cds$Parent,anno.cds$start,anno.cds$strand,decreasing = FALSE),]
  anno.cds.last <- anno.cds[order(anno.cds$Parent,anno.cds$start,anno.cds$strand,decreasing = TRUE),]
  anno.cds.frist <- anno.cds.frist[!duplicated(anno.cds.frist$Parent),]
  anno.cds.last <- anno.cds.last[!duplicated(anno.cds.last$Parent),]
  index.frist <-match(anno.cds.frist$Parent,anno.rna$transcript_id)
  index.last <-match(anno.cds.last$Parent,anno.rna$transcript_id)
  if(length(which(is.na(c(index.frist,index.last))))){
    stop("Can't find cds parent based on input annotation file ")
  }
  anno.cds.frist$utr.start <- anno.rna$start[index.frist]
  anno.cds.frist$utr.end <- anno.cds.frist$start -1
  anno.cds.frist <- anno.cds.frist[which( (anno.cds.frist$utr.end- anno.cds.frist$utr.start) >=0),]

  anno.cds.last$utr.start <- anno.cds.last$end +1
  anno.cds.last$utr.end <- anno.rna$end[index.last]
  anno.cds.last <- anno.cds.last[which((anno.cds.last$utr.end- anno.cds.last$utr.start) >=0),]


  gr.first <- GRanges(seqnames =as.character(anno.cds.frist$Parent) ,
                      ranges =IRanges(start=as.integer(anno.cds.frist$utr.start) ,
                                      end=as.integer(anno.cds.frist$utr.end)),
                      strand =as.character(anno.cds.frist$strand))

  gr.last <- GRanges(seqnames =as.character(anno.cds.last$Parent) ,
                     ranges =IRanges(start=as.integer(anno.cds.last$utr.start) ,
                                     end=as.integer(anno.cds.last$utr.end)),
                     strand =as.character(anno.cds.last$strand))

  gr.exon <- GRanges(seqnames =as.character(anno.exon$Parent) ,
                     ranges =IRanges(start=as.integer(anno.exon$start) ,
                                     end=as.integer(anno.exon$end)),
                     strand =as.character(anno.exon$strand))

  ov.first <- findOverlaps(gr.first,gr.exon)
  ov.last <- findOverlaps(gr.last,gr.exon)
  ov.first <- as.data.frame(ov.first)
  ov.last <- as.data.frame(ov.last)
  colnames(ov.first)<-c("cdsID","exonID")
  colnames(ov.last) <- c("cdsID","exonID")


  ov.first$utr.start <- as.integer(anno.cds.frist$utr.start[ov.first$cdsID])
  ov.first$utr.end <- as.integer(anno.cds.frist$utr.end[ov.first$cdsID])
  ov.first$exon.start <- as.integer(anno.exon$start[ov.first$exonID])
  ov.first$exon.end <- as.integer(anno.exon$end[ov.first$exonID])
  ov.first$utr.start.r <- ov.first$exon.start
  ov.first$utr.end.r <- apply(ov.first[,c("utr.end","exon.end")],1,min)
  five.utr <- anno.exon[ov.first$exonID,]
  five.utr$start <- ov.first$utr.start.r
  five.utr$end <- ov.first$utr.end.r
  if(nrow(five.utr)){
    five.utr$type <- "five_prime_UTR"
    five.utr$type[which(five.utr$strand=="-")] <- "three_prime_UTR"
  }



  ov.last$utr.start <- as.integer(anno.cds.last$utr.start[ov.last$cdsID])
  ov.last$utr.end <- as.integer(anno.cds.last$utr.end[ov.last$cdsID])
  ov.last$exon.start <- as.integer(anno.exon$start[ov.last$exonID])
  ov.last$exon.end <- as.integer(anno.exon$end[ov.last$exonID])
  ov.last$utr.start.r <- apply(ov.last[,c("utr.start","exon.start")],1,max)
  ov.last$utr.end.r <- ov.last$exon.end
  three.utr <- anno.exon[ov.last$exonID,]
  three.utr$start <- ov.last$utr.start.r
  three.utr$end <- ov.last$utr.end.r
  if(nrow(three.utr)){
    three.utr$type <- "three_prime_UTR"
    three.utr$type[which(three.utr$strand=="-")] <- "five_prime_UTR"
  }
  utr <- rbind(three.utr,five.utr)
  utr <- utr[order(utr$Parent,utr$type,utr$start),]
  utr$width <- as.integer(utr$end-utr$start+1)

  #-------------------------------------
  #check result
  #  really.utr <- anno.frame[which(anno.frame$type %in% c("three_prime_UTR","five_prime_UTR")),]
  #  really.utr <- really.utr[order(really.utr$Parent,really.utr$type,really.utr$start),]
  # length(unique(really.utr$Parent))
  # length(unique(utr$Parent))
  # identical(utr$start,really.utr$start)
  # identical(utr$end,really.utr$end)
  # identical(utr$strand,really.utr$strand)
  #  write.table(really.utr,file="really_utr.txt",col.names = TRUE,row.names = FALSE,sep="\t",
  #              quote=FALSE)
  #  write.table(utr,file="build_utr.txt",col.names = TRUE,row.names = FALSE,sep="\t",
  #              quote=FALSE)
  anno.frame <-rbind(anno.frame,utr)
  return(anno.frame)
}

#=========================================================
#------------------------------------------------------
#function:add_rna()
#Adding RNA region
#--------------------------------------------------------
#======================================================
add_rna <- function(anno.frame=NULL){
  anno.exon <- anno.frame[which(anno.frame$type=="exon"),]
  anno.exon.order <- anno.exon[order(anno.exon$gene_id,anno.exon$transcript_id,
                                     anno.exon$strand,anno.exon$start,decreasing = FALSE),]
  anno.exon.rev <- anno.exon.order[nrow(anno.exon.order):1,]

  anno.exon.order.unique <- anno.exon.order[!duplicated(anno.exon.order$transcript_id),]
  anno.exon.rev.order <- anno.exon.rev[!duplicated(anno.exon.rev$transcript_id),]
  anno.rna <- anno.exon.order.unique
  index <- match(anno.rna$transcript_id,anno.exon.rev.order$transcript_id)
  anno.rna$end <- anno.exon.rev.order$end[index]
  anno.rna$Parent <- anno.rna$gene_id
  anno.rna$type <- "mRNA"
  anno.frame <-rbind(anno.frame,anno.rna)
  return(anno.frame)
}
