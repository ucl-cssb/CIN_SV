#!/usr/bin/env Rscript

library(circlize)
library(tidyverse)
library(igraph)
library(ggtree)
library(ape)
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

# Visualize simulated structural variants

############### basic settings #############
chr_info = readRDS("~/Gdrive/git/CIN_SV/data/chr_info.rds")
chr_info_hg38 = readRDS("~/Gdrive/git/CIN_SV/data/chr_info_hg38.rds")

# colors in CN heatmap
# cn_colors1 = c("#6283A9","#bdd7e7","#f0f0f0","#FCAE91", "#B9574E", "#76000D", "#8B0000", "#000000")
# Max CN to show in heatmap
# MAX_CN = 6
# # For absolute CN 0 to 6 (obtained from https://colorbrewer2.org, 5 classes, sequential data)
# cn_colors1 = c("#6283A9","#bdd7e7","#f0f0f0", '#fdcc8a','#fc8d59','#e34a33','#b30000')
MAX_CN = 12
# For absolute CN 0 to 8 (obtained from https://colorbrewer2.org, 8 classes, sequential data)
cn_colors1 = c("#6283A9","#bdd7e7","#f0f0f0", '#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026', '#000000')

# MAX_CN = 7
# cn_colors1 = c("#6283A9","#bdd7e7","#f0f0f0", '#fdd49e', '#fdbb84','#fc8d59','#ef6548','#d7301f','#990000')

# For relative CN -4 to 4
cn_colors2 = c('#08519c','#3182bd', '#6baed6', '#9ecae1', "#f0f0f0", '#fdcc8a','#fc8d59','#e34a33','#b30000')


# CN colors for each haplotype in circos plot
# MAX_CN = 5
# cols = c("#6283A9","#f0f0f0", "#B9574E", "#3b0107")
# cols = c("#bdd7e7","#f0f0f0", '#fdcc8a','#fc8d59','#e34a33','#b30000')
col_fun = colorRamp2(breaks = seq(0:MAX_CN)-1, colors = cn_colors1)

# # CN code in Nature paper from Funnell et al.
# cn_colours <- structure(
#   c(
#     "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
#     "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
#   ),
#   names = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
# )

# map SV colors
# pa <- c("#377eb8","#d95f02","#33a02c")
# names(pa) <- c("DEL","INS","INV")
# pa <- c("#1F78B4","#FF7F00","#33A02C","#E31A1C")
# names(pa) <- c("DEL","DUP","H2HINV","T2TINV")
# same color as shatterseek:	DEL_color='darkorange1',DUP_color='blue1', t2tINV_color="forestgreen",h2hINV_color="black",
# pa <- c("darkorange1","blue1", "black", "forestgreen", "purple")
sv_color <- c("blue1","red", "black", "forestgreen", "purple")
svtype = c("DEL", "DUP", "H2HINV", "T2TINV", "BND")
names(sv_color) <- svtype

# hcn_colors <- c(
#   `Balanced` = "#f7f7f7",
#   `A-Hom` = "#a6dba0",
#   `B-Hom` = "#c2a5cf",
#   `A-Gained` = "#1b7837",
#   `B-Gained` = "#762a83"
# )

hcn_colors <- c("#1b7837", "#a6dba0", "#762a83", "#c2a5cf", "#d5d5d4")
hcn_types <- c("A (gained)", "A (hom)", "B (gained)", "B (hom)", "Balanced")
# color functions for haplotype graph of lines
# '#edf8fb',
# colfunc1<-c('#ccece6','#99d8c9','#66c2a4','#2ca25f','#006d2c')
# colfunc1<-colorRampPalette(c("darkolivegreen1","green4"));
# colfunc1<-colorRampPalette(c("#ccece6","#006d2c"));
#colfunc1<-colorRampPalette(c("#d9f0d3","#00441b"));
#colfunc1<-colorRampPalette(c("#fee0b6","#fee0b6"));
colfunc1<-colorRampPalette(c("lightgrey","lightgrey"));

# colfunc2<-colorRampPalette(c("deepskyblue","dodgerblue4"));
# '#f1eef6',
# colfunc2<-c('#d0d1e6','#a6bddb','#74a9cf','#2b8cbe','#045a8d')
#colfunc2<-colorRampPalette(c('#e7d4e8','#40004b'))
#colfunc2<-colorRampPalette(c('#df9536','#df9536'))
colfunc2<-colorRampPalette(c('darkgrey','darkgrey'))

colfunc3<-colorRampPalette(c("white","gray"));


theme_notxt = theme(legend.position = "none",
               strip.text.x = element_blank(),
               strip.text.y = element_blank(),
               # strip.text.y.left = element_text(size=6, angle = 0),
               strip.text.y.left = element_blank(),
               strip.background = element_blank(),
               axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_blank(),
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank(),
               panel.spacing.y = unit(0, "lines"),
               panel.spacing.x = unit(0, "lines"),
               panel.border = element_rect(color = "#d8d8d8", fill = NA, size = .2),
               panel.background = element_rect(fill = "#f0f0f0")
)
 
              
# not show sample names
theme0 = theme(legend.position = "none",
               strip.text.x = element_text(size=5, angle = 0),
               strip.text.y = element_blank(),
               # strip.text.y.left = element_text(size=6, angle = 0),
               strip.text.y.left = element_blank(),
               strip.background = element_blank(),
               axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_blank(),
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank(),
               panel.spacing.y = unit(0, "lines"),
               panel.spacing.x = unit(0, "lines"),
               panel.border = element_rect(color = "#d8d8d8", fill = NA, size = .2),
               panel.background = element_rect(fill = "#f0f0f0")
)


theme1 = theme(legend.position = "none",
               #strip.text.x = element_blank(),
               strip.text.x = element_text(size=6, angle = 0),
               #strip.text.y = element_blank(),
               strip.text.y.left = element_text(size=6, angle = 0),
               strip.background = element_blank(),
               axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_blank(),
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank(),
               panel.spacing.y = unit(0, "lines"),
               panel.spacing.x = unit(0, "lines"),
               panel.border = element_rect(color = "#d8d8d8", fill = NA, size = .2),
               panel.background = element_rect(fill = "#f0f0f0")
)


# only not show samples on the left
theme2 = theme(legend.position = "none",
               # strip.text.x = element_blank(),
               # strip.text.y = element_blank(),
               # strip.text.y.left = element_text(size=6, angle = 0),
               strip.text.y.left = element_blank(),
               strip.background = element_blank(),
               axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_blank(),
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank(),
               panel.spacing.y = unit(0, "lines"),
               panel.spacing.x = unit(0, "lines"),
               panel.border = element_rect(color = "#d8d8d8", fill = NA, size = .2),
               panel.background = element_rect(fill = "#f0f0f0")
)


############### functions #############
# Plot CNPs of multiple samples as a heatmap for easy comparison
# Format of d_seg: sample, chrom, start, end, cn
plot.cn.heatmap <- function(d_seg, main, type="absolute", theme = theme1, cn_colors = cn_colors1, allele_specific = F, tsize = 12){
  d_seg$cn = round(d_seg$cn)
  d_seg$chrom = factor(d_seg$chrom, levels=paste("",c(1:22, "X"),sep=""))
  d_seg$pos = (d_seg$end + d_seg$start) / 2
  d_seg$width = (d_seg$pos - d_seg$start) * 2
  
  if(type=="absolute"){
    print("Plot absolute copy number")
    d_seg %>% dplyr::mutate(cn=if_else(cn > MAX_CN, MAX_CN, cn)) -> d_seg
    cn_vals = c("0", "1", "2", "3", "4", "5", "6")
  }else{
    print("Plot relative copy number")
    d_seg %>% dplyr::mutate(cn=if_else(cn > max_rcn, max_rcn, cn)) %>% dplyr::mutate(cn=if_else(cn < min_rcn, min_rcn, cn)) -> d_seg
    cn_vals = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4")
  }
  
  d_seg$cn = as.factor(d_seg$cn)
  # unique(d_seg$cn)
  
  p <- ggplot(data=d_seg, aes(x=pos, y=sample, width=width)) +
    geom_tile(aes(fill=cn)) +
    facet_grid(sample ~ chrom, scales="free", space = "free", switch='both') +
    scale_color_manual(values = cn_colors, limits=cn_vals) +
    scale_fill_manual(values = cn_colors, limits=cn_vals) +
    scale_y_discrete(breaks = unique(d_seg$sample), expand = c(0,0))+
    scale_x_discrete(expand = c(0,0)) + theme
  p = p + ggtitle(main) + theme(plot.title = element_text(size=tsize))
  
  return(p)
}


# Plot haplotype CN state of multiple samples as a heatmap for easy comparison
# Format of d_seg: sample, chrom, start, end, cn, type
plot.hcn.heatmap <- function(d_seg, main, hcn_colors, hcn_types, theme = theme1){
  d_seg$chrom = factor(d_seg$chrom, levels=paste("",c(1:22, "X"),sep=""))
  d_seg$pos = (d_seg$end + d_seg$start) / 2
  d_seg$width = (d_seg$pos - d_seg$start) * 2
  
  print("Plot haplotype-specific copy number")
  # hcn_vals = c("Balanced",   "A (hom)",    "B (hom)",    "A (gained)", "B (gained)")
  d_seg$type = as.factor(d_seg$type)
 
  p <- ggplot(data=d_seg, aes(x=pos, y=sample, width=width)) +
    geom_tile(aes(fill=type)) +
    facet_grid(sample ~ chrom, scales="free", space = "free", switch='both') +
    scale_color_manual(values = hcn_colors, limits = hcn_types) +
    scale_fill_manual(values = hcn_colors, limits = hcn_types) +
    scale_y_discrete(breaks = unique(d_seg$sample), expand = c(0,0))+
    scale_x_discrete(expand = c(0,0)) + theme
  p = p + ggtitle(main)
  
  return(p)
}



# seg.cn must be a data frame 
# plot allele-specific/total CN on selected regions
plot_allele_cn <- function(seg.cn, info.type = "AB", ref_start = 0, ref_end = 0, xgap = 100, cn_gap = 1) {
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
  
  if(xgap == 10){
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
    max_B <- max(abs.segments$B, na.rm = TRUE)
    max_CN = max_A 
    if(max_B > max_A){
      max_CN = max_B
    }    
    abs.segments$A[na_As] <- abs.segments$CNt[na_As]
    plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
         y = c(-0.1, (max_CN + 0.1)), type = "n",
         ylab = "Copy number", xlab = xlab,
         xaxt = "n",  yaxt = "n", xaxs = "i")
    axis(labels = seq(0, max_CN, by = cn_gap), at = seq(0, max_CN, by = cn_gap), side = 2, line = 0, las = 1)
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
    axis(labels = seq(min_CNt, max_CNt, by = cn_gap),
         at = seq(min_CNt, max_CNt, by = cn_gap),
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
plot_genome_cn <- function(seg.cn, info.type = "T", plot_pos = T, xgap = 1000) {
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
    # xlab = "Position (Mb)"
    if(xgap == 1000){
      xlab = "Position (Mb)"
      gap = 1e6
      by = 5e7
    }else{
      xlab = "Position (100kb)"
      gap = 1e5
      by = 5e6
    }    
  }else{
    xlab = ""
  }
  
  if (info.type == "AB") {
    na_As <- is.na(abs.segments$A)
    max_A <- max(abs.segments$A, na.rm = TRUE)
    max_B <- max(abs.segments$B, na.rm = TRUE)
    max_CN = max_A 
    if(max_B > max_A){
      max_CN = max_B
    }
    abs.segments$A[na_As] <- abs.segments$CNt[na_As]
    plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
         y = c(-0.1, (max_CN + 0.1)), type = "n",
         ylab = "Copy number", xlab = xlab,
         xaxt = "n",  yaxt = "n", xaxs = "i")
    axis(labels = 0:max_CN, at = 0:max_CN, side = 2, line = 0, las = 1)
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
    # axis(labels = as.character(round(seq(abs.list[[1]]$start.pos[1] / 1e6,
    #                                      abs.list[[1]]$end.pos[nrow(abs.list[[1]])] / 1e6, by = 50), 0)),
    #      at = seq(abs.list[[1]]$abs.start[1],
    #               abs.list[[1]]$abs.end[nrow(abs.list[[1]])], by = 5e7),
    #      outer = FALSE, cex = par("cex.axis") * par("cex"), side = 1,
    #      line = 1)
    axis(labels = as.character(round(seq(abs.list[[1]]$start.pos[1] / gap,
                                         abs.list[[1]]$end.pos[nrow(abs.list[[1]])] / gap, by = 50), 0)),
         at = seq(abs.list[[1]]$abs.start[1],
                  abs.list[[1]]$abs.end[nrow(abs.list[[1]])], by = by),
         outer = FALSE, cex = par("cex.axis") * par("cex"), side = 1,
         line = 1)    
  }
}


make.tree <- function(d, labels = NA) {
  nedge <- nrow(d)
  nleaf <- (nedge + 2)/2
  nnode <- nleaf - 1
  
  mytree <- list()
  mytree$edge <- as.matrix(d[, c(1, 2)])
  mytree$Nnode <- as.integer(nnode)
  
  # find leaves 
  if(is.na(labels)){
    mytree$tip.label <- paste(1:nleaf)
  }
  else{
    mytree$tip.label <- labels
  }
  
  class(mytree) <- "phylo"
  checkValidPhylo(mytree)
  
  plot(mytree)
  
  return(mytree)
}


plot_lineage_tree <- function(infile, sel_nodes){
  lineages <- read_tsv(infile, col_names = T,
                       col_types = cols())
  colnames(lineages) <- c("node", "parent")
  
  # dead_nodes=lineages[which(lineages$flag==-1),]$node
  # curr_nodes=lineages[which(lineages$flag==0),]$node
  
  edges=as.matrix(lineages[,c(2,1)])
  edges=apply(edges, 2, as.character)
  tree=graph_from_edgelist(edges)
  
  tnodes = as.numeric(V(tree)$name)
  ncolor = rep("white", length(V(tree)))
  # ncolor[match(dead_nodes, tnodes)]="lightblue"
  ncolor[match(sel_nodes, tnodes)]="lightgrey"
  

  # main="tumor cell growth"
  main = ""
  fname=str_replace(infile,".tsv","_node.png")
  png(fname)
  ptree = plot.igraph(tree, layout=layout_as_tree, vertex.color=ncolor, main=main)
  ptree
  dev.off()
}


# Assume the same format as ShatterSeek input
get_SV <- function(fname){
  svs = read_tsv(fname, show_col_types = F) 
  if(nrow(svs) == 0) return(NA)
  # all inter-chromosomal translocations are now labelled as BND in Delly
  if("extra" %in% colnames(svs)){
    svs %>% rowwise() %>%  mutate(svclass= strsplit(extra, "=")[[1]][2])  %>% dplyr::rename(chrom1 = chr1, chrom2 = chr2, start1 = coord1, end2 = coord2) %>% mutate(svclass = case_when(chrom1 != chrom2 ~ "BND", chrom1 == chrom2 & strand1 == "-" & strand2 == "+" ~ "DUP", chrom1 == chrom2 & strand1 == "+" & strand2 == "-" ~ "DEL", chrom1 == chrom2 & strand1 == "+" & strand2 == "+" ~ "H2HINV", chrom1 == chrom2 & strand1 == "-" & strand2 == "-" ~ "T2TINV")) %>% as.data.frame() -> svs
  }

  # sv.first and sv.second are dataframes with columns: "chr","start","end"
  sv.first <- svs %>% dplyr::select(chr=chrom1, start=start1) %>% as.data.frame()
  sv.first$chr = paste0("chr", sv.first$chr)
  sv.first$end = sv.first$start

  sv.second <- svs %>% dplyr::select(chr=chrom2, start=end2) %>% as.data.frame()
  sv.second$chr = paste0("chr", sv.second$chr)
  sv.second$end = sv.second$start

  return(list(svs = svs, sv.first=sv.first, sv.second=sv.second))
}


# cutoff not used, showing all regions, cutoff = 5*10e4
get_CN <- function(fname){
  cns = read_tsv(fname, show_col_types = F) 
  
  # for RCK output
  if("extra" %in% colnames(cns)){
    cns %>% rowwise() %>% mutate(cnA = as.numeric(strsplit(strsplit(strsplit(extra, "=")[[1]][2],":")[[1]][3],",")[[1]][1]), cnB = as.numeric(strsplit(strsplit(strsplit(extra, "=")[[1]][2],":")[[1]][4],"}")[[1]][1])) %>% dplyr::rename(chromosome = chr) -> cns
  }

  cns$chromosome = paste0("chr", cns$chromosome)

  cns %>% dplyr::select(chr = chromosome, start, end, value1 = cnA, value2 = cnB) -> cnsel

  # cnsel -> cnsel %>% filter(value1 != 1 | value2 != 1)  %>% filter(end-start > cutoff)

  return(cnsel)
}


plot_SV <- function(fname, ref = "hg38"){
  res = get_SV(fname)
  svs = res$svs
  sv.first = res$sv.first
  sv.second = res$sv.second
  
  col <- sv_color[svs$svclass]
  
  outfile <- str_replace(fname, ".tsv", ".png")
  name <- ""
  png(outfile, height=600, width=600)
  
  circos.initializeWithIdeogram(species=ref)
  circos.genomicLink(sv.first, sv.second, col=col)
  circos.clear()
  
  title(name)
  dev.off()
}


plot_CN_normal <- function(cnv, ref = "hg38", size = 5, col = col_fun){
  cnsel = get_CN(fcnv)
  cnsel = cnsel %>% mutate(value1 = 1, value2 = 1)   # for normal data
  
  outfile <- file.path(dirname(fcnv), "sv_normal.pdf")
  name <- ""
  #png(outfile, height=800, width=800)
  pdf(outfile, height=size, width=size)
  circos.initializeWithIdeogram(species=ref, chromosome.index = paste0("chr", seq(1:22)), plotType = c("ideogram", "labels"))
  circos.genomicHeatmap(cnsel, col = col_fun, side = "inside", border = "white")
  circos.clear()
  
  title(name)
  dev.off() 
}


plot_CN_SV <- function(fsv, fcnv, ref = "hg38", size = 5, col = col_fun, pa = sv_color){
  # get CN data
  cnsel = get_CN(fcnv)
  # cnsel = cnsel %>% mutate(value1 = 1, value2 = 1)   # for normal data
  
  # get SV data
  # print(fsv)
  res = get_SV(fsv)
  
  # output file
  # outfile <- str_replace(fsv, ".tsv", ".png")
  outfile <- str_replace(fsv, ".tsv", ".pdf")
  name <- ""
  #png(outfile, height=800, width=800)
  pdf(outfile, height=size, width=size)
  circos.initializeWithIdeogram(species=ref, chromosome.index = paste0("chr", seq(1:22)), plotType = c("ideogram", "labels"))
  circos.genomicHeatmap(cnsel, col = col_fun, side = "inside", border = "white")
  
  # if(!is.na(res)){
  if(length(res) > 0){    
    svs = res$svs 
    sv.first = res$sv.first
    sv.second = res$sv.second

    df_svt = data.frame(svc = 1:length(pa), svclass = svtype)
    svs$id = 1:nrow(svs)
    svs_num = merge(svs, df_svt, by = c("svclass"), all.x = T)
    svs_num = svs_num[order(svs_num$id), ]
    
    col <- pa[svs_num$svc]
    # print(svs)
    # print(col)
    
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


# not work so far
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

# based on scripts in https://github.com/dmcblab/InfoGenomeR
t_col <- function(color, percent = 1, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (percent)*255,
               names = name)
  return(t.col)
}

# based on scripts in https://github.com/dmcblab/InfoGenomeR
# h indicate haplotype, must be at the sample haplotype
curve_plot <- function(x1,x2,t1,t2,color, h, lwd = sv_lwd){
  if( !is.na(h) && h > 0){
    a = - abs (10/(x1-x2)^2);
  }else{
    a = + abs (10/(x1-x2)^2);
  }
  cat("curve plot, a:", a, "\n")
  b=(t1-t2)/(x1-x2) - a*(x1+x2);
  c=t1-a*x1*x1 - b*x1;
  curve(a*x*x+b*x+c, x1, x2, add=TRUE, col=color, lwd=lwd)
}


# based on scripts in https://github.com/dmcblab/InfoGenomeR
# revise color to be consistent with others
plot_chr <- function(SV, t, fplot){
  t<-cbind("Sample.1",t)
  
  maxt=0;
  for(i in 1:22){
    if(maxt<nrow(t[t[,2]==i,])){
      maxt=nrow(t[t[,2]==i,])
    }
  }	
  maxt=maxt;	
  
  
  pdf(fplot,width=20, height=15)
  
  xmax=max(t[,4])
  plot(c(-10000000,xmax+(maxt+10)*10000000),c(0,120), xlab="", ylab="", col="white",cex.axis=1.5, cex.main=1.5, cex.lab=1.5, cex.sub=1.5,axes=0)
  j=0
  chr=1
  d=data.frame(chrom=0,key=0,value=0, expanded=0)
  dindex=1
  
  xlim=xmax+(maxt+10)*10000000;
  
  colgrad=5
  text(xmax+maxt*10000000 - 105000000, 122, "B",cex=1.5)
  text(xmax+maxt*10000000 - 125000000, 122, "A",cex=1.5)
  
  for(i in 1:colgrad){
    text(xmax+maxt*10000000 - 140000000, 120 - 3.5*i+1.75, i,cex=1.5)
    rect(xmax+maxt*10000000 - 110000000, 120 - 3.5*i , xmax+maxt*10000000-100000000, 120 - 3.5*(i-1), col=colfunc1(colgrad)[i], lwd=1, border=NA);
    rect(xmax+maxt*10000000 - 130000000, 120 - 3.5*i , xmax+maxt*10000000-120000000, 120 - 3.5*(i-1), col=colfunc2(colgrad)[i], lwd=1, border=NA);
  }
  for(i in 1:22){
    if(i==23){  
      it="X"
    }else{
      it=i
    }
    text(-10000000, i*5, it, cex= 1.5);
  }
  
  
  for (i in 1:nrow(t)){
    print(t[i,2])
    ##	if(t[i,2]!="X"){
    chrindex=as.integer(as.character(t[i,2]))*5;
    ##	}else{
    ##		chrindex=22;
    ##	}
    
    ##        if((t[i,2])!="X"){
    if(chrindex != chr ){
      j=1;
    } else { j=j+1;
    }
    #                seg=round(t[i,18]);
    segcol="gray"
    #	if(seg>14){
    #		segcol=colfunc1(14)[seg];
    #	}else{
    #		segcol="red2";
    #	}
    
    point1index=t[i,3]+(j-1)*5000000;
    point2index=t[i,4]+(j-1)*5000000;
    
    if(point2index-point1index<5000000){
      point1index=point1index-1500000;
      point2index=point2index+1500000;
    }
    
    #                if(is.na(t[i,29])==T||(i>1&&t[i,1]==t[i-1,1]&&is.na(t[i-1,29])==T)||(i<nrow(t)&&t[i,1]==t[i+1,1]&&is.na(t[i+1,29])==T) || t$balanced[i] == 1){
    #		if((is.na(t[i,29])==T || (t$balanced[i]  == 1) &&
    #		!(i>1 && i<nrow(t) && t[i-1,1] == t[i,1] && t[i+1,1] == t[i,1] && t[i-1,"phased"]== 1 && t[i+1,"phased"]==1)) 
    #		){
    
    if(0){
      
      ##		if(is.na(t[i,26])==T||t[i,28]==1||(i>1&&t[i,1]==t[i-1,1]&&is.na(t[i-1,26])==T)||(i<nrow(t)&&t[i,1]==t[i+1,1]&&is.na(t[i+1,26])==T)){
      #	                points(point1index,chrindex, pch=21, bg=segcol,cex=0.7, lwd=0.5)
      #			##text(point1index,chrindex-0.2,t[i,3],cex=0.5)
      #	                points(point2index,chrindex, pch=21, bg=segcol,cex=0.7, lwd=0.5)
      ## text(point2index,chrindex-0.4,t[i,4], cex=0.5)
      #	                segments(point1index,chrindex,point2index,chrindex, lwd=0.5)
      rect(point1index, chrindex -0.25 , point2index, chrindex+0.25, col=segcol, lwd=1);
      #                        text(point1index,chrindex,t[i,3],cex=0.25)
      #                        text(point2index,chrindex,t[i,4],cex=0.25)
      
      if(i>1 && t[i-1,2] == t[i,2]){
        if(!pexpanded){
          segments(px, py, point1index , chrindex , col="black", lwd=1);
        }else{
          if(t[i-1,"allele_1"]!=0)
            segments(px1, py1, point1index , chrindex , col="black", lwd=1);
          if(t[i-1,"allele_2"]!=0)
            segments(px2, py2, point1index , chrindex , col="black", lwd=1);
        }
      }
      pexpanded=0;
      px=point2index;
      py=chrindex;
      
      expanded=0;
    }else{			
      
      if(i>1 && t[i-1,2] == t[i,2]){
        if(!pexpanded){
          if(t[i,"allele_1"]!=0)
            segments(px, py, point1index , chrindex+0.5 , col="black", lwd=1);
          if(t[i,"allele_2"]!=0)
            segments(px, py, point1index , chrindex -0.5, col="black", lwd=1);
          
        }else{
          if(t[i-1,"allele_1"]!=0 && t[i,"allele_1"]!=0)
            segments(px1, py1, point1index , chrindex+0.5 , col="black", lwd=1);
          if(t[i-1,"allele_2"]!=0 && t[i,"allele_2"]!=0)
            segments(px2, py2, point1index , chrindex-0.5 , col="black", lwd=1);
        }
      }
      
      
      if(t[i,5]!=0){
        ACN_i=0.5;
        seg=t[i,5];
        if(seg<6){
          segcol=colfunc1(6)[seg]
        }else{
          segcol="green4";
        }
        if(i>1 && t[i-1,2]==t[i,2] && t[i-1,5] == t[i,5]){
          rect(p1_point1index, chrindex +ACN_i -0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1, border=NA);
        }else{
          rect(point1index, chrindex +ACN_i -0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1, border=NA);
        }
        p1_point1index=point1index
        #				text(point1index,chrindex+ACN_i,t[i,3],cex=0.25)
        #				text(point2index,chrindex+ACN_i,t[i,4],cex=0.25)
        
        px1=point2index;
        py1=chrindex+ACN_i;
        #                   points(point1index,chrindex+ACN_i, pch=21, bg=segcol,cex=0.7, lwd=0.5)
        #        	                points(point2index,chrindex+ACN_i, pch=21, bg=segcol,cex=0.7, lwd=0.5)
        #               	        segments(point1index,chrindex+ACN_i,point2index,chrindex+ACN_i, lwd=0.5)
      }
      if(t[i,6]!=0){
        ACN_i=-0.5;
        seg=t[i,6];
        if(seg<6){
          segcol=colfunc2(6)[seg]
        }else{
          segcol="dodgerblue4";
        }
        if(i>1 && t[i-1,2]==t[i,2] && t[i-1,6] == t[i,6]){
          rect(p2_point1index, chrindex +ACN_i -0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1,border=NA);
        }else{
          rect(point1index, chrindex +ACN_i -0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1, border=NA);
        }
        p2_point1index=point1index
        #                               text(point1index,chrindex+ACN_i,t[i,3],cex=0.25)
        #                                text(point2index,chrindex+ACN_i,t[i,4],cex=0.25)
        
        px2=point2index;
        py2=chrindex+ACN_i;
        #                                points(point1index,chrindex+ACN_i, pch=21, bg=segcol,cex=0.7, lwd=0.5)
        #                               points(point2index,chrindex+ACN_i, pch=21, bg=segcol,cex=0.7, lwd=0.5)
        #                              segments(point1index,chrindex+ACN_i,point2index,chrindex+ACN_i, lwd=0.5)
      }
      pexpanded=1;
      expanded=1;
    }
    
    d[dindex,1]=chrindex/5;
    d[dindex,2]=t[i,3];
    d[dindex,3]=point1index;
    d[dindex,4]=expanded;
    dindex=dindex+1;
    d[dindex,1]=chrindex/5;
    d[dindex,2]=t[i,4];
    d[dindex,3]=point2index;
    d[dindex,4]=expanded;
    dindex=dindex+1;
    
    chr=chrindex;
  }
  
  
  for(i in 1:nrow(SV)){
    if(SV[i,2]!="X"){
      SVi2=SV[i,2];
    }else{
      SVi2=22;
    }
    if(SV[i,4]!="X"){
      SVi4=SV[i,4];
    }else{
      SVi4=22;
    }
    
    d1=d[d$chrom==SVi2,];
    chr1=as.integer(as.character(SVi2));
    x1=d1[d1$key==SV[i,3],3];
    x1_expanded=d1[d1$key==SV[i,3],4];
    #		x1_expanded=is.na(SV[i,20]);
    #		print(SV[i,3]);
    d1=d[d$chrom==SVi4,];
    chr2=as.integer(as.character(SVi4));	
    x2=d1[d1$key==SV[i,5],3];
    x2_expanded=d1[d1$key==SV[i,5],4];
    #		x2_expanded=is.na(SV[i,21]);
    
    a=15/(x1-x2)^2
    
    if(x1_expanded==0 && x2_expanded==0){
      if(SV[i,1]=="<DUP>")
        # col="#E8501F"
        col = pa[2]
      if(SV[i,1]=="<DEL>")
        # col="#76BF72"
        col = pa[1]
      if(SV[i,1]=="<INV>"){
        if(SV[i,6]=="5to5") # head to head
          # col="#FFC021"
          col = pa[3]
        if(SV[i,6]=="3to3")
          # col="#1F9C8B"
          col = pa[4]
      }
      if(SV[i,1]=="<TRA>"){
        if(SV[i,6]=="5to5")
          col="#FFC021"
        else if (SV[i,6]=="5to3")
          col="#E8501F"
        else if (SV[i,6]=="3to5")
          col="#76BF72"
        else
          col="#1F9C8B"
      }
      
      
      if(SV[i,1]=="<DUP>")
        curve(-a*(x-x1)*(x-x2)+chr1*5,x1,x2,add=TRUE, col=t_col("red3", 0.1), lwd=1)
      if(SV[i,1]=="<DEL>")
        curve(-a*(x-x1)*(x-x2)+chr1*5,x1,x2,add=TRUE, col="green4", lwd=1)
      if(SV[i,1]=="<INV>")
        if(x1==x2){
          points(x1,chr1*5+0.25, pch=21, col="red",cex=3, lwd=1)
        } else{
          if(SV[i,6]=="5to5")
            curve(-a*(x-x1)*(x-x2)+chr1*5,x1,x2,add=TRUE, col="blue3", lwd=1)
          else
            curve(-a*(x-x1)*(x-x2)+chr1*5,x1,x2,add=TRUE, col="yellow3", lwd=1)
        }
      if(SV[i,1]=="<TRA>"){
        if(SV[i,6]=="5to5")
          segments(x1,chr1*5,x2,chr2*5, lwd=1, col="blue3")
        else if (SV[i,6]=="5to3")
          segments(x1,chr1*5,x2,chr2*5, lwd=1, col="red3")
        else if (SV[i,6]=="3to5")
          segments(x1,chr1*5,x2,chr2*5, lwd=1, col="green3")
        else
          segments(x1,chr1*5,x2,chr2*5, lwd=1, col="yellow3")
        
      }
      
      if(SV[i,1]=="<SEG>"){
        curve_plot(x1,x2,chr1*5,chr2*5, "black", SV[i,7]);
        
      }           
    }else{      
      x1_y_coor=0;
      if(SV[i,7]!=0)
        #		if(!is.na(SV[i,7]))
        x1_y_coor=0.5
      #		if(!is.na(SV[i,8]))
      if(SV[i,8]!=0)
        x1_y_coor=-0.5
      
      x2_y_coor=0;
      if(SV[i,9]!=0)
        x2_y_coor=0.5
      if(SV[i,10]!=0)
        x2_y_coor=-0.5
      
      if(x1_expanded==0)
        x1_y_coor=0;
      if(x2_expanded==0)
        x2_y_coor=0;
      
      if(SV[i,1]=="<SEG>"){
        if(SV[i,1]=="<TRA>" || SV[i,2] !=SV[i,4]){
          segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=1, col=t_col("purple",SV$tag[i]))
        }else{
          curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col("purple",SV$tag[i]), SV[i,7]);
        }
      }else if(SV[i,6]=="5to5"){
        if(SV[i,2] !=SV[i,4]){
          segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=1, col=t_col("#FFC021",SV$tag[i]))
        }else{
          curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col("#FFC021",SV$tag[i]),SV[i,7]);
        }
      }
      else if (SV[i,6]=="5to3"){
        if(SV[i,2] !=SV[i,4]){
          segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=1, col=t_col("#E8501F",SV$tag[i]))
        }else{
          curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col("#E8501F",SV$tag[i]),SV[i,7]);
        }
      }
      else if (SV[i,6]=="3to5"){
        if(SV[i,2] !=SV[i,4]){
          segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=1, col=t_col("#76BF72",SV$tag[i]))
        }else{
          curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col("#76BF72",SV$tag[i]),SV[i,7]);
        }
      }
      else{
        if(SV[i,2] !=SV[i,4]){
          segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=1, col=t_col("#1F9C8B",SV$tag[i]))
        }else{
          curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col("#1F9C8B",SV$tag[i]),SV[i,7]);
        }
      }
      #	}else{
      #
      #		if(!is.na(SV[i,16])){
      #			curve_coor=0.5;
      #		}else{
      #			curve_coor=-0.5;
      #		}
      #		if(SV[i,1]=="<SEG>"){
      #                         curve(-a*(x-x1)*(x-x2)+chr1*5+curve_coor,x1,x2,add=TRUE, col="black", lwd=2)
      #		}else if(SV[i,1]=="<DUP>")
      #               		 curve(-a*(x-x1)*(x-x2)+chr1*5+curve_coor,x1,x2,add=TRUE, col="red3", lwd=2)
      #        	if(SV[i,1]=="<DEL>")
      #                	curve(-a*(x-x1)*(x-x2)+chr1*5+curve_coor,x1,x2,add=TRUE, col="green4", lwd=2)
      #        	if(SV[i,1]=="<INV>")
      #                	if(x1==x2){
      #                	points(x1,chr1*5+0.25, pch=21, col="red",cex=3, lwd=2)
      #                	} else{
      #                        	if(SV[i,6]=="5to5")
      #                                	curve(-a*(x-x1)*(x-x2)+chr1*5+curve_coor,x1,x2,add=TRUE, col="blue3", lwd=2)
      #                        	else
      #                                	curve(-a*(x-x1)*(x-x2)+chr1*5+curve_coor,x1,x2,add=TRUE, col="yellow3", lwd=2)
      #                	}
      #
      #	}
    }
  }
  
  dev.off()
}


# remove CN color, not so clear
plot_chr_noCN <- function(SV, t, fplot){
  t<-cbind("Sample.1",t)
  
  maxt=0;
  for( i in 1:22){
    if(maxt<nrow(t[t[,2]==i,])){
      maxt=nrow(t[t[,2]==i,])
    }
  }	
  maxt=maxt;	
  
  pdf(fplot,width=20, height=15)
  
  xmax=max(t[,4])
  plot(c(-10000000,xmax+(maxt+10)*10000000),c(0,120), xlab="", ylab="", col="white",cex.axis=1.5, cex.main=1.5, cex.lab=1.5, cex.sub=1.5,axes=0)
  j=0
  chr=1
  d=data.frame(chrom=0,key=0,value=0, expanded=0)
  dindex=1
  
  xlim=xmax+(maxt+10)*10000000;
  
  for(i in 1:22){
    it=i
    text(-10000000, i*5, it, cex= 1.5);
  }
  
  # for CN
  for (i in 1:nrow(t)){
    print(t[i,2])
    chrindex=as.integer(as.character(t[i,2]))*5;
    print(chrindex)
    print(chr)
    if(chrindex != chr){
      j=1;
    } else { 
      j=j+1;
    }

    segcol="gray"
    
    point1index=t[i,3]+(j-1)*5000000;
    point2index=t[i,4]+(j-1)*5000000;
    
    if(point2index-point1index<5000000){
      point1index=point1index-1500000;
      point2index=point2index+1500000;
    }
    
    # ydist = 0.5
    ydist = 1 # distance between two haplotypes
    if(i>1 && t[i-1,2] == t[i,2]){
      if(!pexpanded){
        if(t[i,"allele_1"]!=0)
          segments(px, py, point1index, chrindex+ydist, col="black", lwd=1);
        if(t[i,"allele_2"]!=0)
          segments(px, py, point1index, chrindex-ydist, col="black", lwd=1);
        
      }else{
        if(t[i-1,"allele_1"]!=0 && t[i,"allele_1"]!=0)
          segments(px1, py1, point1index, chrindex+ydist , col="black", lwd=1);
        if(t[i-1,"allele_2"]!=0 && t[i,"allele_2"]!=0)
          segments(px2, py2, point1index, chrindex-ydist , col="black", lwd=1);
      }
    }
    
      # allele A
      if(t[i,5]!=0){
        ACN_i=ydist;
        seg=t[i,5];
        segcol=colfunc1(6)[seg]
        if(i>1 && t[i-1,2]==t[i,2] && t[i-1,5] == t[i,5]){
          rect(p1_point1index, chrindex +ACN_i-0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1, border=NA);
        }else{
          rect(point1index, chrindex +ACN_i-0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1, border=NA);
        }
        p1_point1index=point1index
        px1=point2index;
        py1=chrindex+ACN_i;
      }
    
    # allele B
    if(t[i,6]!=0){
      ACN_i=-ydist;
      seg=t[i,6];
      segcol=colfunc2(6)[seg]
      if(i>1 && t[i-1,2]==t[i,2] && t[i-1,6] == t[i,6]){
        rect(p2_point1index, chrindex+ACN_i-0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1,border=NA);
      }else{
        rect(point1index, chrindex+ACN_i-0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1, border=NA);
      }
      p2_point1index=point1index
      px2=point2index;
      py2=chrindex+ACN_i;
    }
    pexpanded=1;
    expanded=1;
    
    d[dindex,1]=chrindex/5;
    d[dindex,2]=t[i,3];
    d[dindex,3]=point1index;
    d[dindex,4]=expanded;
    dindex=dindex+1;
    d[dindex,1]=chrindex/5;
    d[dindex,2]=t[i,4];
    d[dindex,3]=point2index;
    d[dindex,4]=expanded;
    dindex=dindex+1;
    
    chr=chrindex;
  }
  
  # for SV
  for(i in 1:nrow(SV)){
    SVi2=SV[i,2];
    SVi4=SV[i,4];
    
    d1=d[d$chrom==SVi2,];
    chr1=as.integer(as.character(SVi2));
    x1=d1[d1$key==SV[i,3],3];
    x1_expanded=d1[d1$key==SV[i,3],4];
    d1=d[d$chrom==SVi4,];
    chr2=as.integer(as.character(SVi4));	
    x2=d1[d1$key==SV[i,5],3];
    x2_expanded=d1[d1$key==SV[i,5],4];

    a=15/(x1-x2)^2
    
    print(x1_expanded)
    print(x2_expanded)
    
    x1_y_coor=0;
    if(SV[i,7]==0)
      x1_y_coor=ydist
    else
      x1_y_coor=-ydist
    
    x2_y_coor=0;
    if(SV[i,8]==0)
      x2_y_coor=ydist
    else
      x2_y_coor=-ydist
    
    if(x1_expanded==0)
      x1_y_coor=0;
    if(x2_expanded==0)
      x2_y_coor=0;
      
      sv_lwd = 2
      if(SV[i,1]=="<TRA>" || SV[i,2] != SV[i,4]){
        segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=sv_lwd, col=t_col(sv_color["BND"]))
        #curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col(sv_color["BND"]),SV[i,7]);
      }else if(SV[i,6]=="5to5"){
          curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col(sv_color["H2HINV"]),SV[i,7], lwd = sv_lwd);
      }
      else if (SV[i,6]=="5to3"){
          curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col(sv_color["DEL"]),SV[i,7], lwd = sv_lwd);
      }
      else if (SV[i,6]=="3to5"){
          curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col(sv_color["DUP"]),SV[i,7], lwd = sv_lwd);
      }
      else{
          curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col(sv_color["T2TINV"]),SV[i,7], lwd = sv_lwd);
      }
  }
  
  dev.off()
}



# based on code of "signals", not tested
make_copynumber_heatmap <- function(copynumber,
                                    clones,
                                    colvals = cn_colours,
                                    legendname = "Copy Number",
                                    library_mapping = NULL,
                                    clone_pal = NULL,
                                    sample_label_idx = 1,
                                    cutoff = NULL,
                                    maxf = 1.0,
                                    plotcol = "state",
                                    plotfrequency = FALSE,
                                    show_legend = TRUE,
                                    show_library_label = TRUE,
                                    show_clone_label = TRUE,
                                    show_clone_text = TRUE,
                                    chrlabels = TRUE,
                                    SV = NULL,
                                    labeladjust = -1,
                                    nticks = 4,
                                    Mb = TRUE,
                                    annotation_height = NULL, 
                                    annofontsize = 14,
                                    na_col = "white",
                                    linkheight = 5,
                                    str_to_remove = NULL,
                                    anno_width = 0.4,
                                    rasterquality = 15,
                                    ...) {
  
  if (class(colvals) == "function"){
    leg_params <- list(nrow = 3,
                       direction = "vertical",
                       labels_gp = grid::gpar(fontsize = annofontsize-1),
                       title_gp = grid::gpar(fontsize = annofontsize, fontface = "bold"),
                       legend_gp = grid::gpar(fontsize = annofontsize-1))
  } else {
    leg_params <- list(nrow = 3,
                       direction = "vertical",
                       at = names(colvals),
                       labels_gp = grid::gpar(fontsize = annofontsize-1),
                       title_gp = grid::gpar(fontsize = annofontsize, fontface = "bold"),
                       legend_gp = grid::gpar(fontsize = annofontsize-1))
  }
  
  copynumber_hm <- ComplexHeatmap::Heatmap(
    name = legendname,
    as.matrix(copynumber),
    col = colvals,
    na_col = na_col,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    bottom_annotation = make_bottom_annot(copynumber, chrlabels = chrlabels, 
                                          Mb = Mb, nticks = nticks, 
                                          annotation_height = annotation_height, 
                                          labeladjust = labeladjust,
                                          annofontsize = annofontsize, 
                                          linkheight = linkheight),
    left_annotation = make_left_annot(copynumber, clones, anno_width = anno_width,
                                      library_mapping = library_mapping, clone_pal = clone_pal, show_clone_label = show_clone_label, show_clone_text = show_clone_text,
                                      idx = sample_label_idx, show_legend = show_legend, show_library_label = show_library_label,annofontsize = annofontsize, 
                                      str_to_remove = str_to_remove
    ),
    heatmap_legend_param = leg_params,
    top_annotation = make_top_annotation_gain(copynumber,
                                              cutoff = cutoff, maxf = maxf,
                                              plotfrequency = plotfrequency, plotcol = plotcol, SV = SV
    ),
    use_raster = TRUE,
    raster_quality = rasterquality,
    ...
  )
  return(copynumber_hm)
}



