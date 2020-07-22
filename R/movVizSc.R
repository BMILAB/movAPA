#' Visualize PACs of a gene in single cells
#'
#' movVizSC visualizes PACs of a gene in single cells.
#'
#' This function is useful for visualizing distributions of expression levels of PACs among single cells and cell types.

#' @usage movVizSC(scPACds, gene, cellGroupName, cellGroupColors=NULL, txdb = NULL, PAwidth = 50, showRatio = FALSE)
#' @param movVizSC a PACdataset storing single cell PACs. The PACds@anno should have columns chr/strand/coord.
#' If there is no colData in PACds, then the sample label will be set as groupN.
#' @param gene a gene to plot.
#' @param cellGroupName group name of cell type in column names of scPACds@colData.
#' @param cellGroupColors a named vector denoting colors for cell groups, e.g., c('SC="yellow",ES="red",RS="blue").
#'                        its names must be the same as the labels of cell groups.
#' @param txdb txdb annotation, a character of file name for .gff3 or .rda, or a list. See parseGenomeAnnotation().
#' @param PAwidth expand the width of a PAC by PAwidth for the plot.
#' @param showRatio default is FALSE. If TRUE then show PAC ratio instead of PAC count.
#' @return a plot with the top panel showing expression levels of PACs in single cells, the middle panel showing expression levels or ratios of PACs among cell types, the bottome panel showing the gene model and PAC.
#' @examples
#' data("scPACds")
#' load('txdbmm10.rda')
#' #gene <- unique(scPACds@anno$gene)[2]
#' gene='ENSMUSG00000019969'
#' movVizSC(scPACds, gene, cellGroupName="celltype", cellGroupColors=NULL, txdb=txdbmm10)
#' movVizSC(scPACds, gene, cellGroupName="celltype", cellGroupColors = c(SC="yellow",ES="red",RS="blue"), txdb=txdbmm10, showRatio = T)
#' movVizSC(scPACds, gene, cellGroupName="celltype", PAwidth=100, txdb=NULL, showRatio = T)
#' @name movVizSC
#' @seealso [movViz()] to visualize a PACdataset.
#' @export
movVizSC <- function(scPACds, gene, cellGroupName, cellGroupColors=NULL, txdb = NULL, PAwidth = 50, showRatio = FALSE){

  if (!require('millefy', quietly = TRUE, warn.conflicts = FALSE)) {
    stop("R package millefy is not installed!")
  }

  print(sprintf("Begin Plot: %s", Sys.time()))

  if (!(cellGroupName %in% colnames(scPACds@colData))) stop(cat(cellGroupName,"not in colnames(scPACds@colData)!"))
  cellGroupName <- factor(scPACds@colData[, cellGroupName])
  cts <- levels(cellGroupName)

  # heatmap color(PA count)
  breaks <- 100
  pal <- colorRampPalette(c('white','blue'))
  pal <- pal(breaks)

  # cell type colors, like c(SC='blue', RS='red')
  if (is.null(cellGroupColors)) {
    if (length(cts)<=9) {
      cols=brewer.pal(9, "Set1")
      cellGroupColors=cols[1:length(cts)]
      names(cellGroupColors)=cts
    } else {
      stop("There are >9 cell groups, please reduce the cell group number before the plot!")
    }

  } else {
    if ( length(intersect(names(cellGroupColors), cts))!=length(cts) ) {
      stop("Names of cellGroupColors not the same as labels of cell groups!")
    }
    cellGroupColors=cellGroupColors[cts]
  }

  bal <- list()
  for (j in names(cellGroupColors)) {
    # j <- 1
    bal[[j]] <- colorRampPalette(c('white', cellGroupColors[j]))
    bal[[j]] <- bal[[names(cellGroupColors[j])]](breaks)
  }

  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))
  heights <- c(12,2)
  # gene_id turn to gene_name
  # same.gene <- select(org.Mm.eg.db, keys = gene, columns = c("SYMBOL", "GENENAME"),
  #                     keytype = "ENSEMBL")
  PAcounts <- scPACds@counts
  ex.samegene <- t(PAcounts[which(scPACds@anno$gene %in% gene),])

  ex.samegene <- na.omit(ex.samegene)

  ex.samegene.old <- ex.samegene

  group_colors <- cellGroupColors

  e <- try({
    isAppropriateColorLabels(group_colors, cellGroupName)
  }, silent=TRUE)
  if(class(e) == "try-error"){
    cat("Group colors must be a named vector.\n")
  }

  cellGroupColors1 <- makeColorLabels(group_colors, cellGroupName)

  ex.samegene[is.na(ex.samegene)] <- 0 #OK?
  min_value <- min(ex.samegene)
  max_value <- max(ex.samegene)

  sel_danger <- rowSums(ex.samegene) == 0
  sel_safe <- !(sel_danger)
  sel_safe[which(sel_danger)[1]] <- TRUE
  dc1_value <- runDiffusionMap(ex.samegene[sel_safe,])
  dc1 <- numeric(nrow(ex.samegene))
  dc1[sel_safe] <- dc1_value
  dc1[!sel_safe] <- dc1[which(sel_danger)[1]]

  neworder <- orderWithinGroups(dc1, cellGroupName)
  ex.samegene <- as.matrix(ex.samegene[neworder,])
  cellGroupColors2 <- cellGroupColors1[neworder]
  # neworder <- order(dc1, decreasing = T)
  # ex.samegene <- as.matrix(ex.samegene[neworder,])

  # set canvas
  grid.newpage()
  pushViewport(plotViewport(c(1,0,1,0)))
  # design plot
  pushViewport(viewport(layout = grid.layout(3, 1, heights = heights)))

  # first plot
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  plotScHeatmap(ex.samegene, breaks, pal, min_value, max_value, txdb, cellGroupName, cellGroupColors2, PAwidth, gene, scPACds)
  popViewport()

  # second plot
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  if(!showRatio){
    plotSumPACount(ex.samegene.old, txdb, cellGroupName, cellGroupColors, PAwidth, gene, scPACds)
  }else{
    plotRatioPACount(ex.samegene.old, txdb, bal, cellGroupName, cellGroupColors, PAwidth, gene, scPACds)
  }

  popViewport()

  # third plot
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  plotGeneAnnotation(txdb, scPACds, gene, PAwidth)
  popViewport()

  # fourth plot
  # pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
  # plotAxisTrack(scPACds, gene, PAwidth)
  # popViewport()


  popViewport()
  popViewport()


}
plotScHeatmap <- function(ex.samegene, breaks, pal, min_value, max_value, txdb, cellGroupName, cellGroupColors2, PAwidth, gene, scPACds){
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))

  pushViewport(plotViewport(c(0.5,0,0.5,0)))
  pushViewport(viewport(layout = grid.layout(2, length(col_ratio),
                                             heights = unit(c(2,12), c("lines", "null") ),
                                             widths=col_ratio)))

  # plot PA heatmap
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
  plotHeatmapSC(ex.samegene ,breaks, pal, max_value, PAwidth, txdb, scPACds, gene)
  popViewport()


  # add a label about expressions of PAs
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
  plotColorKey(pal[1:breaks], min_value, max_value)
  popViewport()

  # add PA info label
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  plotColorLabel(rev(cellGroupColors2))
  popViewport()

  # add cell type name
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  plotGroupLabel(cellGroupName)
  popViewport()

  popViewport()
  popViewport()

  # return(list(ex.samegene, cellGroupColors2))

}

isAppropriateColorLabels <- function(group_colors, groups){
  is.vector(group_colors) & !any(is.na(names(group_colors))) & all(groups %in% names(group_colors))
}

makeColorLabels <- function(group_colors, groups){
  group_colors[levels(groups)[groups]]
}
runDiffusionMap<- function(mat){
  e <- try({
    dif <- DiffusionMap(as.ExpressionSet(as.data.frame(mat)));
    cat(sprintf("Eigenvalue of DC1: %f\n", eigenvalues(dif)[1]));
    res <- dif@eigenvectors[, "DC1"];
  }, silent=TRUE)

  if(class(e) == "try-error"){
    cat("There was a problem when running diffusion map. Trying PCA instead...\n")
    e <- try({
      pca =  prcomp(log10(mat+1),
                    center = TRUE,
                    scale. = TRUE);
      cat(sprintf("The standard deviations of PC1: %f\n", pca$sdev[1]));
      res <- pca$x[, "PC1"];
    }, silent=TRUE)

    if(class(e) == "try-error"){
      res <- rep(1, nrow(mat))
    }else{
      return(res)
    }
  }else{
    return(res)
  }
}


orderWithinGroups <- function(values, groups){
  unlist(as.vector(mapply(function(x,y){y[order(x, decreasing = T)]}, split(values, groups), split(seq_along(values), groups))))
}

plotHeatmapSC <- function(mat, breaks, pal, fixed_max, PAwidth, txdb, scPACds, gene){
  nr <- nrow(mat)
  nc <- ncol(mat)

  x_chr <- scPACds@anno$chr[which(scPACds@anno$gene %in% gene)]
  x_start <- scPACds@anno$coord[which(scPACds@anno$gene %in% gene)]
  x_end <- scPACds@anno$coord[which(scPACds@anno$gene %in% gene)]+1

  if(is.null(txdb)){
    ftr <- scPACds@anno$ftr[which(scPACds@anno$gene %in% gene)]
    ftr_start <- scPACds@anno$ftr_start[which(scPACds@anno$gene %in% gene)]
    ftr_end <- scPACds@anno$ftr_end[which(scPACds@anno$gene %in% gene)]

    minstart <- min(ftr_start)
    maxend <- max(ftr_end)

    minstart <- minstart-((maxend-minstart)*0.1)
    maxend <- maxend+((maxend-minstart)*0.1)

    for (i in 1:nc) {
      # x = unit((gene.start[i]-minstart)/(maxend - minstart + 1), "npc")
      # x = unit(((i/nc)-0.14), "npc")
      grid.raster(matrix(pal[as.numeric(cut2(mat[,i], breaks, fixed_max = fixed_max))], nrow = nr),
                  x = unit(((x_start[i]-minstart)/(maxend - minstart + 1)), "npc"), default.units = "native",
                  y = unit(0.5, "npc"), interpolate = F, width = unit(PAwidth/800, "npc"), height = 1)
    }
  }else{
    dt_gtf <- gtfToDtExon(txdb)
    transcript_ids <- c()
    for (n in 1:length(x_start)) {
      transcript.id <- dt_gtf[seqnames == x_chr[n] & start <= x_end[n] & end >= x_start[n], unique(transcript_id)]
      transcript_ids <- c(transcript_ids, transcript.id)
    }
    transcript_ids <- unique(transcript_ids)

    dt_gtf %>% setkey(transcript_id)
    dt_gtf_local <- dt_gtf[transcript_ids]

    minstart <- min(dt_gtf_local$start)
    maxend <- max(dt_gtf_local$end)

    minstart <- minstart-((maxend-minstart)*0.1)
    maxend <- maxend+((maxend-minstart)*0.1)

    for (i in 1:nc) {
      # x = unit((gene.start[i]-minstart)/(maxend - minstart + 1), "npc")
      # x = unit(((i/nc)-0.14), "npc")
      grid.raster(matrix(pal[as.numeric(cut2(mat[,i], breaks, fixed_max = fixed_max))], nrow = nr),
                  x = unit(((x_start[i]-minstart)/(maxend - minstart + 1))+0.025, "npc"), default.units = "native",
                  y = unit(0.5, "npc"), interpolate = F, width = unit(PAwidth/800, "npc"), height = 1)
    }
  }

}
cut2 <- function (x, breaks, labels = NULL, include.lowest = FALSE, right = TRUE,
                  dig.lab = 3L, ordered_result = FALSE, fixed_max, ...)
{
  if (!is.numeric(x))
    stop("'x' must be numeric")
  if (length(breaks) == 1L) {
    if (is.na(breaks) || breaks < 2L)
      stop("invalid number of intervals")
    nb <- as.integer(breaks + 1)
    dx <- diff(rx <- range(x, na.rm = TRUE))
    rx[2L] <- fixed_max
    if (dx == 0) {
      dx <- abs(rx[1L])
      breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000,
                        length.out = nb)
    }
    else {
      breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
      breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] +
                               dx/1000)
    }
  }
  else nb <- length(breaks <- sort.int(as.double(breaks)))
  if (anyDuplicated(breaks))
    stop("'breaks' are not unique")
  codes.only <- FALSE
  if (is.null(labels)) {
    for (dig in dig.lab:max(12L, dig.lab)) {
      ch.br <- formatC(0 + breaks, digits = dig, width = 1L)
      if (ok <- all(ch.br[-1L] != ch.br[-nb]))
        break
    }
    labels <- if (ok)
      paste0(if (right)
        "("
        else "[", ch.br[-nb], ",", ch.br[-1L], if (right)
          "]"
        else ")")
    else paste("Range", seq_len(nb - 1L), sep = "_")
    if (ok && include.lowest) {
      if (right)
        substr(labels[1L], 1L, 1L) <- "["
      else substring(labels[nb - 1L], nchar(labels[nb -
                                                     1L], "c")) <- "]"
    }
  }
  else if (is.logical(labels) && !labels)
    codes.only <- TRUE
  else if (length(labels) != nb - 1L)
    stop("lengths of 'breaks' and 'labels' differ")
  code <- .bincode(x, breaks, right, include.lowest)
  if (codes.only)
    code
  else factor(code, seq_along(labels), labels, ordered = ordered_result)
}
plotColorKey <- function(colors, min_color, max_color){
  pushViewport(viewport(layout = grid.layout(
    1,2,
    widths = unit(c(4,1), c("null","null"))
  )))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))


  grid.raster(matrix(colors,nrow=1), y = 0.25, interpolate = F, width = 1, height = 0.5, vjust = 0)
  grid.rect(x=0, y = 0.25, width=1, height = 0.5, gp = gpar(fill = NA, col="black"), vjust = 0, hjust = 0)

  # grid.text(min_color, x=-0.1, y=0.24, hjust = 0, vjust = 1, gp = gpar(cex = 0.4))
  # if(max_color > 0 & min_color < 0){
  #   grid.text(0, x=(0-min_color)/(max_color-min_color), y=0.24, hjust = 0, vjust = 1, gp = gpar(cex = 0.4))
  # }
  grid.text(min_color, x=0, y=0.24, hjust = 0, vjust = 1, gp = gpar(cex = 0.8))

  grid.text(max_color, x=1, y=0.24, hjust = 1, vjust = 1, gp = gpar(cex = 0.8))

  popViewport()
  popViewport()
}
plotColorLabel <- function(cellGroupColors){

  pushViewport(viewport(yscale = c(1, length(cellGroupColors)+1)))

  for(i in seq_along(cellGroupColors)){
    grid.rect(x=unit(0, "npc"), y=i, width = unit(1, "npc"), height = 1,
              default.units = "native", vjust = 0, hjust = 0,
              gp = gpar(fill = cellGroupColors[i], col = NA))
  }

  popViewport()
}
plotGroupLabel <- function(groups){

  lv <- levels(factor(groups))

  # pushViewport(viewport(yscale = c(1, length(lv)+1)))
  pushViewport(viewport(layout = grid.layout(
    length(lv),1,
    heights = unit(sapply(lv, function(x){table(factor(groups))[x]}), "null"),
    widths = unit(1, "null")
  )))

  for(i in seq_along(lv)){
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))

    grid.text(lv[i],
              x = 0.9,
              y= 0.5,
              hjust = 1,
              gp = gpar(cex = 0.7))

    popViewport()
  }

  popViewport()
}
plotSumPACount <- function(ex.samegene, txdb, cellGroupName, cellGroupColors, PAwidth, gene, scPACds){
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))
  pushViewport(viewport(layout = grid.layout(1, length(col_ratio), widths=col_ratio)))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  plotSum(ex.samegene, cellGroupColors, cellGroupName, PAwidth, txdb, scPACds, gene)
  popViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
  com <- levels(as.factor(cellGroupName))
  x.min <- c(NULL)
  x.max <- c(NULL)
  for (m in 1:length(com)) {
    if(ncol(ex.samegene) == 1){
      mt <- sum(ex.samegene[which(cellGroupName == com[m]),])
    }else{
      mt <- as.matrix(colSums(ex.samegene[which(cellGroupName == com[m]),]))
    }

    x.min <- c(x.min, min(mt))
    x.max <- c(x.max, max(mt))
  }
  plotYvalue(length(com), x.max, x.min)
  popViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  plotColorLabel(rev(unique(cellGroupColors)))
  popViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  plotGroupLabel(levels(cellGroupName))
  popViewport()

  popViewport()

}
plotSum <- function(mat, cellGroupColors, cellGroupName, PAwidth, txdb, scPACds, gene){
  # nr <- nrow(mat)
  nc <- ncol(mat)
  ct <- levels(cellGroupName)
  # max_value <- max(colSums(mat))
  x_chr <- scPACds@anno$chr[which(scPACds@anno$gene %in% gene)]
  x_start <- scPACds@anno$coord[which(scPACds@anno$gene %in% gene)]
  x_end <- scPACds@anno$coord[which(scPACds@anno$gene %in% gene)]+1

  if(is.null(txdb)){
    ftr <- scPACds@anno$ftr[which(scPACds@anno$gene %in% gene)]
    ftr_start <- scPACds@anno$ftr_start[which(scPACds@anno$gene %in% gene)]
    ftr_end <- scPACds@anno$ftr_end[which(scPACds@anno$gene %in% gene)]

    minstart <- min(ftr_start)
    maxend <- max(ftr_end)
  }else{
    dt_gtf <- gtfToDtExon(txdb)
    transcript_ids <- c()
    for (n in 1:length(x_start)) {
      transcript.id <- dt_gtf[seqnames == x_chr[n] & start <= x_end[n] & end >= x_start[n], unique(transcript_id)]
      transcript_ids <- c(transcript_ids, transcript.id)
    }
    transcript_ids <- unique(transcript_ids)

    dt_gtf %>% setkey(transcript_id)
    dt_gtf_local <- dt_gtf[transcript_ids]

    minstart <- min(dt_gtf_local$start)
    maxend <- max(dt_gtf_local$end)
  }

  # pushViewport(viewport(xscale = c(minstart-((maxend-minstart)*0.1), maxend+((maxend-minstart)*0.1)), yscale = c(0,1), clip = "off"))

  minstart <- minstart-((maxend-minstart)*0.1)
  maxend <- maxend+((maxend-minstart)*0.1)

  pushViewport(viewport(layout = grid.layout(length(ct),1)))

  breaks <- 100
  # bottom <- 0.5
  # h <- 1

  # pushViewport(viewport(xscale = c(minstart, maxend), yscale = c(0,length(transcript_ids)+1), clip = "off"))
  for (i in 1:length(ct)) {
    if(nc == 1){
      ctmat <- sum(mat[which(cellGroupName == ct[i]),])
    }else{
      ctmat <- as.matrix(colSums(mat[which(cellGroupName == ct[i]),]))
    }
    # ctmat <- as.matrix(mat[which(cellGroupName == ct[i]),])
    max_value <- max(ctmat)

    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1, clip = "on"))
    # pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1, clip = "on"))


    for (j in 1:nc) {
      # [as.numeric(cut2(ctmat[,j], breaks, fixed_max = fixed_max))]
      # y = bottom+h*(i-1)
      if(is.null(txdb)){
        grid.rect(x = unit(((x_start[j]-minstart)/(maxend - minstart + 1)), "npc"),
                y=unit(0, "npc"), width = unit(PAwidth/800, "npc"), height = 2*(ctmat[j]/max_value),
                gp = gpar(fill = cellGroupColors[ct[i]], col = NA))
        # gp = gpar(fill = NA, col = "grey")
      }else{
        grid.rect(x = unit(((x_start[j]-minstart)/(maxend - minstart + 1))+0.025, "npc"),
                  y=unit(0, "npc"), width = unit(PAwidth/800, "npc"), height = 2*(ctmat[j]/max_value),
                  gp = gpar(fill = cellGroupColors[ct[i]], col = NA))
      }
      # grid.raster(matrix(cellGroupColors[which(names(cellGroupColors) == ct[i])], nrow =  max_value+1),
      #             x = unit((x_start[j]-minstart)/(maxend - minstart + 1), "npc"),
      #             y = unit(0, "npc"), interpolate = F,
      #             width = unit(PAwidth/800, "npc"), height = ctmat[j])
    }
    popViewport()
  }
  popViewport()
  # popViewport()

}
plotYvalue <- function(nr, max_color, min_color){

  pushViewport(viewport(layout = grid.layout(nr,1)))
  for(i in 1:nr){
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(max_color[i],
              x = 0.05, y= 0.95,
              hjust = 0,
              vjust = 1,
              gp = gpar(cex = 0.5)
    )
    grid.text(min_color[i],
              x = 0.05, y= 0.05,
              hjust = 0,
              vjust = 0,
              gp = gpar(cex = 0.5)
    )
    popViewport()
  }
  popViewport()
}
plotGeneAnnotation <- function(txdb, scPACds, gene, PAwidth){
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))

  pushViewport(viewport(layout = grid.layout(5, length(col_ratio), widths=col_ratio)))

  chr <- scPACds@anno$chr[which(scPACds@anno$gene %in% gene)]
  start <- scPACds@anno$coord[which(scPACds@anno$gene %in% gene)]
  end <- scPACds@anno$coord[which(scPACds@anno$gene %in% gene)]+1

  # gene annotation
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  if(is.null(txdb)){
    plotGeneModels2(scPACds, gene, chr, start, end)
  }else{
    plotGeneModels(txdb, PAwidth, chr, start, end)
  }
  popViewport()

  # annotation label
  # pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  # plotGroupLabel("GENECODE")
  # popViewport()

  # axis location
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
  plotAxisTrack(scPACds, gene)
  popViewport()

  popViewport()
}
plotGeneModels <- function(txdb, PAwidth, x_chr, x_start, x_end, show_transcript_id = TRUE){

    dt_gtf <- gtfToDtExon(txdb)
    transcript_ids <- c()
    for (n in 1:length(x_start)) {
      transcript.id <- dt_gtf[seqnames == x_chr[n] & start <= x_end[n] & end >= x_start[n], unique(transcript_id)]
      transcript_ids <- c(transcript_ids, transcript.id)
    }
    transcript_ids <- unique(transcript_ids)

    dt_gtf %>% setkey(transcript_id)
    dt_gtf_local <- dt_gtf[transcript_ids]

    minstart <- min(dt_gtf_local$start)
    maxend <- max(dt_gtf_local$end)

    minstart <- minstart-((maxend-minstart)*0.1)
    maxend <- maxend+((maxend-minstart)*0.1)

    if(length(transcript_ids)==0){
      cat("No gene models in this region\n")
    }else{

      pushViewport(viewport(xscale = c(minstart, maxend), yscale = c(0,length(transcript_ids)+1), clip = "off"))


      bottom <- 1
      h <- 1
      hh <- 0.3
      len_arrow = unit(hh/2, "lines")
      for(i in 1:length(transcript_ids)){

        # i <- 1
        starts <- dt_gtf_local$start[which(dt_gtf_local$transcript_id == transcript_ids[i])]
        ends <- dt_gtf_local$end[which(dt_gtf_local$transcript_id == transcript_ids[i])]


        grid.lines(x=c(starts, ends), y=rep(bottom+h*(i-1),2),
                   default.units = "native", gp = gpar(col = "#554C4D"))
        for(k in 1:length(starts)){

          for(j in 1:length(x_start)){
            # print(starts[j]) #### debug
            # j <- 1
            # k <- 1
            if(x_start[j] %in% seq(starts[k], ends[k])){
              len <- PAwidth/1000
              grid.rect(
                x = unit(((x_start[j]-minstart)/(maxend - minstart + 1)), "npc"),
                y = bottom+h*(i-1),
                width = unit(len, "npc"),
                # width = unit(1, "npc"),
                height = hh,
                gp = gpar(fill = "#ff9900"),
                default.units = "native",
                hjust = 0,

              )
            }


          }
        }

        if("strand" %in% colnames(dt_gtf_local)){
          strand <- dt_gtf_local$strand[which(dt_gtf_local$transcript_id == transcript_ids[i])]
          for (l in 1:length(strand)) {
            if(strand[l] == "+"){
            gene_text <- sprintf("> %s ", transcript_ids[i])
          }else if (strand[l] == "-"){
            # gene_text <- sprintf("< %s (%s)", transcript_ids[i], transcript_names[i])
            gene_text <- sprintf("< %s ", transcript_ids[i])
          }else{
            gene_text <- sprintf("%s (%s)", transcript_ids[i])
          }
          }

        }else{
          gene_text <- sprintf("%s (%s)", transcript_ids[i])
        }

        # xpos <- (min(x_end, max(ends)) + max(x_start, min(starts)) - 2* x_start)/2/(x_end - x_start + 1)
        if(max(ends) < maxend | min(starts) > minstart){
          xpos <- (min(max(ends), max(x_end)) + max(min(starts), min(x_start)) - 2* min(starts))/(2*2)/(max(ends) - min(starts) + 1)

        }else{
          xpos <- (min(max(ends), max(x_end)) + max(min(starts), min(x_start)) - 2* min(starts))/2/(max(ends) - min(starts) + 1)
        }
        if(show_transcript_id){
          grid.text(gene_text,
                    x = unit(xpos, "npc"), y=unit(rep(bottom+h*(i-1)+hh,2), "native"),
                    hjust = 0.5, vjust = 0,
                    gp = gpar(fontsize = 20, cex = 0.4)
          )
        }

      }
      popViewport()
    }

}

gtfToDtExon <- function(txdb = txdbmm10){
  dt_gtf <- txdbmm10$anno.need %>% subset(type == "exon")
  # dt_gtf <- txdb$anno.need
  # dt_gtf$end[which(dt_gtf$strand == "+")] <- dt_gtf$end[which(dt_gtf$strand == "+")]+1000
  # dt_gtf$start[which(dt_gtf$strand == "-")] <- dt_gtf$start[which(dt_gtf$strand == "-")]+1000

  as.data.table(dt_gtf[,c("seqnames", "start", "end", "strand", "type", "gene_id", "transcript_id")])
}

plotGeneModels2 <- function(scPACds, gene, x_chr, x_start, x_end){

  ftr <- scPACds@anno$ftr[which(scPACds@anno$gene %in% gene)]
  ftr_start <- scPACds@anno$ftr_start[which(scPACds@anno$gene %in% gene)]
  ftr_end <- scPACds@anno$ftr_end[which(scPACds@anno$gene %in% gene)]

  minstart <- min(ftr_start)
  maxend <- max(ftr_end)

  pushViewport(viewport(xscale = c(minstart-((maxend-minstart)*0.1), maxend+((maxend-minstart)*0.1)), yscale = c(0,1), clip = "off"))

  minstart <- minstart-((maxend-minstart)*0.1)
  maxend <- maxend+((maxend-minstart)*0.1)

  for(i in 1:length(ftr)){

      # i <- 1
      starts <- ftr_start[i]
      ends <- ftr_end[i]


      for(k in 1:length(starts)){

        for(j in 1:length(x_start)){
          # print(starts[j]) #### debug
          # j <- 1
          # k <- 1
          if(x_start[j] %in% seq(starts[k], ends[k])){
            len <- 0.05
            # grid.rect(
            #   x = unit(((x_start[j]-minstart)/(maxend - minstart + 1)), "npc"),
            #   y = unit(0.5,"npc"),
            #   width = unit(len, "npc"),
            #   height = 0.5,
            #   gp = gpar(fill = "#ff9900"),
            #   default.units = "native",
            #   hjust = 0,
            # )

            grid.lines(
              x = unit(((x_start[j]-minstart)/(maxend - minstart + 1)), "npc"),
              y = unit(c(0.3,0),"npc"),
              gp = gpar(fill = "#ff9900", col = "#ff9900"),
              arrow = arrow(length = unit(0.1, "inches"),
                            ends = "last", type = "closed")
            )
          }


        }
      }
      grid.lines(x=c(starts, ends), y=unit(0, "npc"),
                 default.units = "native", gp = gpar(col = "#554C4D"))

    }
    popViewport()
  }

plotAxisTrack <- function(scPACds, gene){
  # col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))

  ftr <- scPACds@anno$ftr[which(scPACds@anno$gene %in% gene)]
  ftr_start <- scPACds@anno$ftr_start[which(scPACds@anno$gene %in% gene)]
  ftr_end <- scPACds@anno$ftr_end[which(scPACds@anno$gene %in% gene)]
  PAstrand <- scPACds@anno$strand[which(scPACds@anno$gene %in% gene)]

  pushViewport(viewport(layout = grid.layout(length(ftr))))

  pushViewport(viewport(xscale = c(min(ftr_start), max(ftr_end)), clip = "off"))
  grid.xaxis(gp = gpar(cex = 0.6))
  if(length(unique(ftr)) != 1){
    for(i in seq_along(ftr)){
    pushViewport(viewport(layout = grid.layout(i)))
    grid.text(sprintf("%s %s:%d-%d", ftr[i], PAstrand, ftr_start[i], ftr_end[i]),
              x = 0.5,
              y= 0.5+0.5*(i-1),
              hjust = 0.5,
              gp = gpar(cex = 0.7))
    popViewport()

    }
  }else{
    grid.text(sprintf("%s %s:%d-%d", ftr, PAstrand, ftr_start, ftr_end),
              x = 0.5,
              y= 0.5,
              hjust = 0.5,
              gp = gpar(cex = 0.7))
  }



  # popViewport()
  # popViewport()

  popViewport()
}

plotRatioPACount <- function(ex.samegene, txdb, bal, cellGroupName, cellGroupColors, PAwidth, gene, scPACds){
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))
  pushViewport(viewport(layout = grid.layout(1, length(col_ratio), widths=col_ratio)))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  plotRatio(ex.samegene, bal, cellGroupName, cellGroupColors, PAwidth, txdb, scPACds, gene)
  popViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
  com <- levels(as.factor(cellGroupName))
  plotYvalue(length(com), rep(1,length(com)), rep(0,length(com)))
  popViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  plotColorLabel(rev(cellGroupColors))
  popViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  plotGroupLabel(levels(cellGroupName))
  popViewport()

  popViewport()

}
plotRatio <- function(mat, bal, cellGroupName, cellGroupColors, PAwidth, txdb, scPACds, gene){
  nc <- ncol(mat)
  ct <- levels(cellGroupName)
  # max_value <- max(colSums(mat))

  x_chr <- scPACds@anno$chr[which(scPACds@anno$gene %in% gene)]
  x_start <- scPACds@anno$coord[which(scPACds@anno$gene %in% gene)]
  x_end <- scPACds@anno$coord[which(scPACds@anno$gene %in% gene)]+1

  if(is.null(txdb)){
    ftr <- scPACds@anno$ftr[which(scPACds@anno$gene %in% gene)]
    ftr_start <- scPACds@anno$ftr_start[which(scPACds@anno$gene %in% gene)]
    ftr_end <- scPACds@anno$ftr_end[which(scPACds@anno$gene %in% gene)]

    minstart <- min(ftr_start)
    maxend <- max(ftr_end)
  }else{
    dt_gtf <- gtfToDtExon(txdb)
    transcript_ids <- c()
    for (n in 1:length(x_start)) {
      transcript.id <- dt_gtf[seqnames == x_chr[n] & start <= x_end[n] & end >= x_start[n], unique(transcript_id)]
      transcript_ids <- c(transcript_ids, transcript.id)
    }
    transcript_ids <- unique(transcript_ids)

    dt_gtf %>% setkey(transcript_id)
    dt_gtf_local <- dt_gtf[transcript_ids]

    minstart <- min(dt_gtf_local$start)
    maxend <- max(dt_gtf_local$end)
  }

  minstart <- minstart-((maxend-minstart)*0.1)
  maxend <- maxend+((maxend-minstart)*0.1)

  pushViewport(viewport(layout = grid.layout(length(ct),1)))

  breaks <- 100

  if(nc != 1){
    t.ratio <- as.data.frame(NULL)
    for (i in 1:length(ct)) {
      ctmat <- as.matrix(colSums(mat[which(cellGroupName == ct[i]),]))
      pair.ratio <- ctmat/sum(ctmat)
      pair.ratio <- as.matrix(pair.ratio)
      pair.ratio <- na.omit(pair.ratio)

      # ctmat <- as.matrix(mat[which(cellGroupName == ct[i]),])
      max_value <- 1

      pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1, clip = "on"))

      for (j in 1:nc) {
        # [as.numeric(cut2(ctmat[,j], breaks, fixed_max = fixed_max))]
        # y = bottom+h*(i-1)
        if(is.null(txdb)){
          grid.rect(x = unit(((x_start[j]-minstart)/(maxend - minstart + 1)), "npc"),
                    y=unit(0, "npc"), width = unit(PAwidth/800, "npc"), height = 2*(pair.ratio[j]/max_value),
                    gp = gpar(fill = cellGroupColors[ct[i]], col = NA))
          # gp = gpar(fill = NA, col = "grey")
        }else{
          grid.rect(x = unit(((x_start[j]-minstart)/(maxend - minstart + 1))+0.025, "npc"),
                    y=unit(0, "npc"), width = unit(PAwidth/800, "npc"), height = 2*(pair.ratio[j]/max_value),
                    gp = gpar(fill = cellGroupColors[ct[i]], col = NA))
        }
      }
      popViewport()
    }


  #   for (i in 1:length(ct)) {
  #     cell.name <- as.character(scPACds@colData$group[which(scPACds@colData$celltype == ct[i])])
  #     ctmat <- as.matrix(mat[which(rownames(mat) %in% cell.name),])
  #     pair.ratio <- apply(ctmat, 1, function(m){
  #       t.ratio <- rbind(t.ratio,m/sum(m));
  #       names(t.ratio) <- colnames(ctmat);
  #       return(t.ratio);
  #     })
  #     pair.ratio <- do.call("rbind",pair.ratio)
  #     pair.ratio <- as.matrix(pair.ratio)
  #     pair.ratio <- na.omit(pair.ratio)
  #
  #   # ctmat <- as.matrix(mat[which(cellGroupName == ct[i]),])
  #     max_value <- max(pair.ratio)
  #   # pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1, yscale = c(0, max_value), clip = "on"))
  #   pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1, clip = "on"))
  #
  #   for (j in 1:nc) {
  #     # [as.numeric(cut2(ctmat[,j], breaks, fixed_max = fixed_max))]
  #     # y = bottom+h*(i-1)
  #
  #     # grid.rect(x = unit(((x_start[j]-minstart)/(maxend - minstart + 1)), "npc"),
  #     #           y=unit(0.5, "npc"), width = unit(PAwidth/800, "npc"), height = ctmat[j],
  #     #           gp = gpar(fill = unlist(bal[ct[i]]), col = NA)
  #               # gp = gpar(fill = NA, col = "grey")
  #     # )
  #     if(sum(pair.ratio[,j]) == 0){
  #       next;
  #     }else{
  #       if (is.null(txdb)) {
  #         grid.raster(matrix(unlist(bal[ct[i]])[as.numeric(cut2(pair.ratio[,j], breaks, fixed_max = max_value))], nrow =  nrow(pair.ratio)),
  #                 x = unit((x_start[j]-minstart)/(maxend - minstart + 1), "npc"),
  #                 y = unit(0.5, "npc"), interpolate = F, default.units = "native",
  #                 width = unit(PAwidth/800, "npc"), height = 1)
  #       }else{
  #         grid.raster(matrix(unlist(bal[ct[i]])[as.numeric(cut2(pair.ratio[,j], breaks, fixed_max = max_value))], nrow =  nrow(pair.ratio)),
  #                     x = unit(((x_start[j]-minstart)/(maxend - minstart + 1))+0.025, "npc"),
  #                     y = unit(0.5, "npc"), interpolate = F, default.units = "native",
  #                     width = unit(PAwidth/800, "npc"), height = 1)
  #       }
  #
  #     }
  #
  #   }
  #   popViewport()
  # }
  }else{
    print("Warning: showRatio=TRUE is only applicable for APA gene.")
  }


  popViewport()

}
