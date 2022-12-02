#!/usr/bin/env Rscript

library(circlize)
library(tidyverse)

# Visualize simulated structural variants

# seg.cn must be a data frame 
# plot allele-specific/total CN on selected regions
plot_allele_cn <- function(seg.cn, info.type = "AB", ref_start = 0, ref_end = 0) {
  seg.cn$CNt = seg.cn$A + seg.cn$B
  
  chr.order <- unique(seg.cn$chromosome)
  seg.list <- split(x = seg.cn, f = seg.cn$chromosome)
  seg.list <- seg.list[order(order(chr.order))]
  
  seg.max <- lapply(X = seg.list, FUN = function(x) x[nrow(x), "end.pos"])
  seg.pos <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
  #seg.max <- cumsum(as.numeric(do.call(rbind, seg.max)))
  seg.max <- cumsum(do.call(rbind, (seg.max)))
  chr.offset <- 0
  for (i in 1:length(seg.pos)){
    seg.pos[[i]] <- seg.pos[[i]] + chr.offset
    colnames(seg.pos[[i]]) <- c("abs.start", "abs.end")
    chr.offset <- seg.max[i]
  }
  seg.max <- sapply(X = seg.pos, FUN = function(x) x[nrow(x), "abs.end" ])
  abs.list <- mapply(cbind, seg.list, seg.pos, SIMPLIFY = FALSE)
  abs.segments <- do.call(rbind, abs.list)
  
  if(ref_start > 0){
    xlab = "Position (10kb)"
    gap = 1e4
    by = 5e5
  }else{
    xlab = "Position (100kb)"
    gap = 1e5
    by = 5e6
  }
  
  if (info.type == "AB") {
    na_As <- is.na(abs.segments$A)
    max_A <- max(abs.segments$A, na.rm = TRUE)
    abs.segments$A[na_As] <- abs.segments$CNt[na_As]
    plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
         y = c(-0.1, (max_A + 0.1)), type = "n",
         ylab = "Copy number", xlab = xlab,
         xaxt = "n",  yaxt = "n", xaxs = "i")
    axis(labels = 0:max_A, at = 0:max_A, side = 2, line = 0, las = 1)
    segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
             y0 = (abs.segments$B - 0.1), y1 = (abs.segments$B - 0.1),
             col = "blue", lwd = 5, lend = 1)
    segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
             y0 = (abs.segments$A + 0.1), y1 = (abs.segments$A + 0.1),
             col = "red", lwd = 5, lend = 1)
  } else {
    min_CNt <- min(abs.segments$CNt, na.rm = TRUE)
    max_CNt <- max(abs.segments$CNt, na.rm = TRUE)
    plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
         y = c(min_CNt, max_CNt), type = "n",
         ylab = "Copy number", xlab = xlab,
         xaxt = "n", yaxt = "n", xaxs = "i")
    axis(labels = min_CNt:max_CNt,
         at = min_CNt:max_CNt,
         side = 2, line = 0, las = 1)
    segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
             y0 = abs.segments$CNt, y1 = abs.segments$CNt, col = "red",
             lwd = 5, lend = 1)
  }
  abline(v = c(0, seg.max), lty = 3)
  if(ref_start > 0){
    abline(v = c(ref_start, ref_end), lty = 2)
  }
  axis(labels = as.character(round(seq(abs.list[[1]]$start.pos[1] / gap,
                                       abs.list[[1]]$end.pos[nrow(abs.list[[1]])] / gap, by = 50), 0)),
       at = seq(abs.list[[1]]$abs.start[1],
                abs.list[[1]]$abs.end[nrow(abs.list[[1]])], by = by),
       outer = FALSE, cex = par("cex.axis") * par("cex"), side = 1,
       line = 1)
  
  for (i in 1:length(abs.list)){
    max.pos <- nrow(abs.list[[i]])
    mtext(chr.order[i], side = 3, line = 0,
          at = sum(abs.list[[i]]$abs.start[1],
                   abs.list[[i]]$abs.end[max.pos]) / 2)
  }
}


# seg.cn must be a data frame
# plot allele-specific/total CN across the whole genome
plot_genome_cn <- function(seg.cn, info.type = "T", plot_pos = T) {
  chr.order <- sort(unique(seg.cn$chromosome))
  seg.list <- split(x = seg.cn, f = seg.cn$chromosome)
  seg.list <- seg.list[order(order(chr.order))]
  
  seg.max <- lapply(X = seg.list, FUN = function(x) x[nrow(x), "end.pos"])
  seg.pos <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
  # seg.max <- cumsum(as.numeric(do.call(rbind, seg.max)))
  seg.max <- cumsum(do.call(rbind, (seg.max)))
  chr.offset <- 0
  for (i in 1:length(seg.pos)){
    seg.pos[[i]] <- seg.pos[[i]] + chr.offset
    colnames(seg.pos[[i]]) <- c("abs.start", "abs.end")
    chr.offset <- seg.max[i]
  }
  seg.max <- sapply(X = seg.pos, FUN = function(x) x[nrow(x), "abs.end" ])
  abs.list <- mapply(cbind, seg.list, seg.pos, SIMPLIFY = FALSE)
  abs.segments <- do.call(rbind, abs.list)
  if(plot_pos){
    xlab = "Position (Mb)"
  }else{
    xlab = ""
  }
  
  if (info.type == "AB") {
    na_As <- is.na(abs.segments$A)
    max_A <- max(abs.segments$A, na.rm = TRUE)
    abs.segments$A[na_As] <- abs.segments$CNt[na_As]
    plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
         y = c(-0.1, (max_A + 0.1)), type = "n",
         ylab = "Copy number", xlab = xlab,
         xaxt = "n",  yaxt = "n", xaxs = "i")
    axis(labels = 0:max_A, at = 0:max_A, side = 2, line = 0, las = 1)
    segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
             y0 = (abs.segments$B - 0.1), y1 = (abs.segments$B - 0.1),
             col = "blue", lwd = 5, lend = 1)
    segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
             y0 = (abs.segments$A + 0.1), y1 = (abs.segments$A + 0.1),
             col = "red", lwd = 5, lend = 1)
  } else {
    min_CNt <- min(abs.segments$CNt, na.rm = TRUE)
    max_CNt <- max(abs.segments$CNt, na.rm = TRUE)
    plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
         y = c(min_CNt, max_CNt), type = "n",
         ylab = "Copy number", xlab = xlab,
         xaxt = "n", yaxt = "n", xaxs = "i")
    axis(labels = min_CNt:max_CNt,
         at = min_CNt:max_CNt,
         side = 2, line = 0, las = 1)
    segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
             y0 = abs.segments$CNt, y1 = abs.segments$CNt, col = "red",
             lwd = 5, lend = 1)
  }
  abline(v = c(0, seg.max), lty = 3)
  for (i in 1:length(abs.list)){
    max.pos <- nrow(abs.list[[i]])
    mtext(chr.order[i], side = 3, line = 0,
          at = sum(abs.list[[i]]$abs.start[1],
                   abs.list[[i]]$abs.end[max.pos]) / 2)
  }
  
  if(plot_pos){
    axis(labels = as.character(round(seq(abs.list[[1]]$start.pos[1] / 1e6,
                                         abs.list[[1]]$end.pos[nrow(abs.list[[1]])] / 1e6, by = 50), 0)),
         at = seq(abs.list[[1]]$abs.start[1],
                  abs.list[[1]]$abs.end[nrow(abs.list[[1]])], by = 5e7),
         outer = FALSE, cex = par("cex.axis") * par("cex"), side = 1,
         line = 1)
  }
}


# Assume the same format as ShatterSeek input
get_SV <- function(fname){
  svs = read_tsv(fname, show_col_types = F)
  if(nrow(svs) == 0) return(NA)
  # svs %>% rowwise() %>%  mutate(svclass= strsplit(extra, "=")[[1]][2]) %>% as.data.frame() -> svs

  # sv.first and sv.second are dataframes with columns: "chr","start","end"
  sv.first <- svs %>% select(chr=chrom1, start=start1) %>% as.data.frame()
  sv.first$chr = paste0("chr", sv.first$chr)
  sv.first$end = sv.first$start

  sv.second <- svs %>% select(chr=chrom2, start=end2) %>% as.data.frame()
  sv.second$chr = paste0("chr", sv.second$chr)
  sv.second$end = sv.second$start

  return(list(svs = svs, sv.first=sv.first, sv.second=sv.second))
}


# cutoff not used, showing all regions
get_CN <- function(fname, cutoff = 5*10e4){
  cns = read_tsv(fname, show_col_types = F)

  # cns %>% rowwise() %>% mutate(cnA = as.numeric(strsplit(strsplit(strsplit(extra, "=")[[1]][2],":")[[1]][3],",")[[1]][1]), cnB = as.numeric(strsplit(strsplit(strsplit(extra, "=")[[1]][2],":")[[1]][4],"}")[[1]][1])) -> cns

  cns$chromosome = paste0("chr", cns$chromosome)

  cns %>% select(chr = chromosome, start, end, value1 = cnA, value2 = cnB) -> cnsel

  # cnsel -> cnsel %>% filter(value1 != 1 | value2 != 1)  %>% filter(end-start > cutoff)

  return(cnsel)
}


plot_SV <- function(fname){
  res = get_SV(fname)
  svs = res$svs
  sv.first = res$sv.first
  sv.second = res$sv.second
  
  # same color as shatterseek:	DEL_color='darkorange1',DUP_color='blue1', t2tINV_color="forestgreen",h2hINV_color="black",
  # pa <- c("darkorange1","blue1", "black", "forestgreen", "purple")
  pa <- c("blue1","red", "black", "forestgreen", "purple")
  names(pa) <- c("DEL", "DUP", "H2HINV", "T2TINV", "BND")
  col <- pa[svs$svclass]
  
  outfile <- str_replace(fname, ".tsv", ".png")
  name <- ""
  png(outfile, height=600, width=600)
  

  circos.initializeWithIdeogram(species="hg38")
  circos.genomicLink(sv.first, sv.second, col=col)
  circos.clear()
  
  title(name)
  dev.off()
}


plot_CN_SV <- function(fsv, fcnv, cutoff = 5*10e4){
  # get CN data
  cnsel = get_CN(fcnv, cutoff)
  # cnsel = cnsel %>% mutate(value1 = 1, value2 = 1)   # for normal data
  
  # CN colors for each haplotype
  MAXCN = 5
  # cols = c("#6283A9","#f0f0f0", "#B9574E", "#3b0107")
  cols = c("#bdd7e7","#f0f0f0", '#fdcc8a','#fc8d59','#e34a33','#b30000')
  col_fun = colorRamp2(breaks = seq(0:MAXCN)-1, colors = cols)

  # output file
  outfile <- str_replace(fsv, ".tsv", ".png")
  name <- ""
  jpeg(outfile, height=800, width=800)
  
  circos.initializeWithIdeogram(species="hg38", chromosome.index = paste0("chr", seq(1:22)))

  circos.genomicHeatmap(cnsel, col = col_fun, side = "inside", border = "white")
  
  # get SV data
  res = get_SV(fsv)
  if(!is.na(res)){
    svs = res$svs
    sv.first = res$sv.first
    sv.second = res$sv.second
    
    # SV colors
    # pa <- c("#377eb8","#d95f02","#33a02c")
    # names(pa) <- c("DEL","INS","INV")
    # pa <- c("#1F78B4","#FF7F00","#33A02C","#E31A1C")
    # names(pa) <- c("DEL","DUP","H2HINV","T2TINV")
    pa <- c("blue1","red", "black", "forestgreen", "purple")
    svtype <- c("DEL", "DUP", "H2HINV", "T2TINV", "BND")
    names(pa) <- c("DEL", "DUP", "H2HINV", "T2TINV", "BND")    
    col_fun_sv = colorRamp2(svtype, pa)
    col <- pa[svs$svclass]
    
    # lgd_heatmap = Legend(at = seq(0:3)-1, col_fun = col_fun, legend_gp = gpar(col = 0:3), title_position = "topleft", title = "Copy number")
    # lgd_links = Legend(at = c("DEL","DUP","H2HINV","T2TINV"), type = "lines", itle_position = "topleft", title = "Links")
    # lgd_list_vertical = packLegend(lgd_heatmap, lgd_links)
    # lgd_list_vertical
    
    circos.genomicLink(sv.first, sv.second, col=col)    
  }
  
  circos.clear()

  #draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

  title(name)
  dev.off()
}



plot_with_zoom <- function(){
  extend_chromosomes = function(bed, chromosome, prefix = "zoom_") {
    zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
    zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
    rbind(bed, zoom_bed)
  }
  
  
  normal_chr_index = 1:22
  zoomed_chr_index = 23
  prefix = "zoom_"
  zoomed_chr = c("chr22")
  zoomed_index = paste0(prefix, c("chr22"))
  
  cytoband = read.cytoband()
  cytoband_df = cytoband$df %>% filter(V1 %in% chromosome)
  chromosome = cytoband$chromosome[normal_chr_index]
  
  xrange = c(cytoband$chr.len[normal_chr_index], cytoband$chr.len[zoomed_chr])
  
  # normalize in normal chromsomes and zoomed chromosomes separately
  sector.width = c(xrange[normal_chr_index] / sum(xrange[normal_chr_index]), 
                   xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index])) 
  
  circos.par(start.degree = 90)
  circos.initializeWithIdeogram(extend_chromosomes(cytoband_df, zoomed_chr), 
                                sector.width = sector.width)
  
  bed = generateRandomBed(500)
  circos.genomicTrack(extend_chromosomes(cnsel, zoomed_chr),
                      panel.fun = function(region, value, ...) {
                        circos.genomicPoints(region, value, pch = 16, cex = 0.3)
                      })
  
  # Add a link from original chromosome to the zoomed chromosome 
  circos.link(zoomed_chr, get.cell.meta.data("cell.xlim", sector.index = zoomed_chr),
              zoomed_index, get.cell.meta.data("cell.xlim", sector.index = zoomed_index),
              col = "#00000020", border = NA)
  circos.clear()
  
  
  
}
