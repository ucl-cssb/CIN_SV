#!/usr/bin/env Rscript

# install 
# install.packages("devtools")
# require(devtools)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BiocGenerics")
# BiocManager::install("graph")
# BiocManager::install("S4Vectors")
# BiocManager::install("GenomicRanges")
# BiocManager::install("IRanges")
# install_github("parklab/ShatterSeek", dependencies = T)

suppressMessages(library(optparse))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(tidyverse))
suppressMessages(library(ShatterSeek))
suppressMessages(library(tools))


chr_info = readRDS("/Users/bl0033/Gdrive/git/CIN_SV/data/chr_info.rds")
chr_info_hg38 = readRDS("/Users/bl0033/Gdrive/git/CIN_SV/data/chr_info_hg38.rds")


#' Plot chromothripsis regions, based on ShatterSeek code
#' This function serves to plot chromothripsis regions
#' in order to facilitate the revision of 
#' candidate chromothripsis regions
#' 
#' @param ShatterSeek_output the output of the function shatterseek
#' @param chr chromosome for which the plot will be generated (note that only the region where there is a cluster of interleaved SVs will be shown)
#' @param BAF B allele frequencies (BAF)
#' @param sample_name name of the sample to be shown in the table
#' @param DEL_color colour to show the deletion-like SVs
#' @param DUP_color colour to show the duplication-like SVs
#' @param t2tINV_color colour to show the t2tINV SVs
#' @param h2hINV_color colour to show the h2hINV SVs
#' @param arc_size size of the arcs representing intrachromosomal SVs
#' @param genome reference genome (hg19 or hg38)
#' @return a list containing ggplot objects (chromosome ideogram, SVs, CN profile, and table with information about the region)
#' 
plot_chromothripsis_local <- function(ShatterSeek_output, chr=chr,BAF=NULL,sample_name="",
                                      DEL_color='blue1',DUP_color='red',
                                      t2tINV_color="forestgreen",h2hINV_color="black",
                                      arc_size=.2, genome = "hg19"){
  chromNames = c( paste0('chr', c(1:22, 'X')) , c(1:22, 'X'))
  if ( !(as.character(chr) %in% chromNames)){stop("Chromosome not valid")}
  
  common_ggplot2 <- theme_bw() + theme(axis.text.x=element_text(size=7,angle=0),
                                       axis.text.y=element_text(size=7),
                                       axis.title.y=element_text(size=7),
                                       axis.title.x=element_blank(),
                                       legend.position="none",
                                       legend.text = element_text(size=7),
                                       legend.key = element_blank(),
                                       plot.margin=unit(c(0.1,0.1,0,0.1),"cm"),
                                       plot.title=element_blank(),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(), 
                                       legend.title=element_blank(),
                                       plot.background = element_blank(),
                                       axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                       axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))  
  #-------------------------------------------------------------------------------------------------------
  summ <- ShatterSeek_output@chromSummary %>% filter(!is.na(start)) %>% filter(chrom == chr)
  
  cand = gsub("chr","",chr)
  chr=paste("chr",cand,sep="")
  summary = ShatterSeek_output@chromSummary
  candidate_chrs <- ShatterSeek_output@chromSummary$chrom 
  cluster_sizes <- sapply(ShatterSeek_output@detail$connComp,length)
  # get the SVs for the cluster of SVs in this chromosome
  cand_clust_size <- ShatterSeek_output@chromSummary$clusterSize[ ShatterSeek_output@chromSummary$chrom == cand]
  idx = which(cluster_sizes == cand_clust_size) 
  SVsnow <- ShatterSeek_output@detail$SV##[ as.numeric(unlist(ShatterSeek_output@detail$connComp[idx])) ,]
  CNVsnow <- ShatterSeek_output@detail$CNV  
  SVsnow <- unique(SVsnow[SVsnow$chrom1 == cand, ]) # remove if there are more
  CNVsnow <- CNVsnow[CNVsnow$chrom == cand, ] # remove if there are more
  
  
  y1=4; y2=12
  df = SVsnow
  df$y1 = rep(y1,nrow(df)); df$y2 = rep(y2,nrow(df))
  min_coord=0;max_coord=10
  if(nrow(df)!=0){
    min_coord = min(df$pos1)
    max_coord = max(df$pos2) # will show more regions than chromothripsis
    df$diff = abs(df$pos1 - df$pos2)
    #max_coord = max(summ$end)  # some SVs may be beyond this region
    #df$diff = abs(df$pos1 - summ$end)
    df$curv = 1 - (df$diff / max(df$diff))
    max_diff = max(df$diff)
    # df$curv[(df$diff / max_diff) > 0.2] <- .15
    # df$curv[(df$diff / max_diff) > 0.8] <- .08
    # df$curv[(df$diff / max_diff) < 0.2] <- 1
    df$curv[(df$diff / max_diff) > 0.2] <- .06
    df$curv[(df$diff / max_diff) > 0.8] <- .03
    df$curv[(df$diff / max_diff) < 0.2] <- 1
  }
  d = data.frame(x=c(min_coord),y=c(1),leg=c("DEL","DUP","t2tINV","h2hINV"))
  idx = c()
  
  idx1=which(CNVsnow$start <= min_coord & CNVsnow$end > min_coord) 
  idx2=which(CNVsnow$start >= min_coord & CNVsnow$end <= max_coord) 
  idx3=which(CNVsnow$end >= max_coord & CNVsnow$start < max_coord) 
  idx=c(idx1, idx2, idx3) 
  
  CNVsnow = CNVsnow[idx,]
  CNVsnow[1,]$start = min_coord
  CNVsnow[nrow(CNVsnow),]$end = max_coord
  
  SV_plot = ggplot(d,aes(x=x,y=y)) +
    geom_point(colour="white") + ylim(0,y2+5) + common_ggplot2  +
    geom_line(data=rbind(d,d),aes(x=x,y=y,colour=leg)) 
  
  #-------------------------------------------------------------------------------------------------------
  # Interchromosomal SVs
  #-------------------------------------------------------------------------------------------------------
  inter <- ShatterSeek_output@detail$SVinter   
  inter <- inter[which(inter$chrom1 == cand | inter$chrom2 == cand), ] 
  
  if (nrow(inter)>0){
    inter$SVtype = factor(inter$SVtype,levels=c("DEL","DUP","h2hINV","t2tINV"))
    inter$y = rep(0,nrow(inter))
    inter$y[which(inter$SVtype %in% c("DUP","DEL"))] = 4
    inter$y[!(inter$SVtype %in% c("DUP","DEL"))] = 12
    inter$type_SV = rep("",nrow(inter))
    inter$type_SV[which(inter$strand1 == "-" & inter$strand2 == "-")] = "t2tINV"
    inter$type_SV[which(inter$strand1 == "-" & inter$strand2 == "+")] = "DUP"
    inter$type_SV[which(inter$strand1 == "+" & inter$strand2 == "-")] = "DEL"
    inter$type_SV[which(inter$strand1 == "+" & inter$strand2 == "+")] = "h2hINV"
    inter$SVtype = inter$type_SV; inter$type_SV=NULL
    inter$colour = rep("",nrow(inter))
    inter$colour[which(inter$SVtype == "DUP")] = DUP_color
    inter$colour[which(inter$SVtype == "DEL")] = DEL_color
    inter$colour[which(inter$SVtype == "h2hINV")] = h2hINV_color
    inter$colour[which(inter$SVtype == "t2tINV")] = t2tINV_color
    
    #inter = data.frame(pos = c(inter$pos1,inter$pos2), y=c(inter$y,inter$y), SVtype=c(inter$SVtype,inter$SVtype) )
    inter = data.frame(chr = c(inter$chrom1,inter$chrom2), pos = c(inter$pos1,inter$pos2), y=c(inter$y,inter$y), SVtype=c(inter$SVtype,inter$SVtype) )
    inter = inter[inter$chr == cand, c(2:4)]
    SV_plot = SV_plot + geom_point(data=inter,size=1,alpha=1,aes(x=pos,y=as.numeric(y),colour=SVtype))
  }
  #-------------------------------------------------------------------------------------------------------
  
  SV_plot = SV_plot + geom_hline(yintercept=y1,size=0.5) + geom_hline(yintercept=y2,size=0.5) 
  
  if(nrow(df)>300){options(expressions= 100000)} 
  #to avoid "Error: evaluation nested too deeply: infinite recursion / options(expressions=)?"
  
  now = df[df$SVtype == "DUP",]
  if (nrow(now) > 0){
    for (i in 1:nrow(now)){
      SV_plot = SV_plot + geom_curve( size=arc_size,data = now[i,] , aes(x = pos1, y = y1, xend = pos2, yend = y1), curvature = now$curv[i],colour=DUP_color,ncp=8)
    }
    idx= c(idx,1)
  }
  SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y1)) + geom_point(data=now,size=.5,aes(x=pos2,y=y1))
  
  
  now = df[df$SVtype == "DEL",]
  if (nrow(now) > 0){
    for (i in 1:nrow(now)){
      SV_plot = SV_plot + geom_curve( size=arc_size,data = now[i,] , aes(x = pos1, y = y1, xend = pos2, yend = y1), curvature = -1*now$curv[i],colour=DEL_color) 
    }
    idx= c(idx,2)
  }
  SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y1)) + geom_point(data=now,size=.5,aes(x=pos2,y=y1))
  
  now = df[df$SVtype == "t2tINV",]
  if (nrow(now) > 0){
    for (i in 1:nrow(now)){
      SV_plot = SV_plot + geom_curve( size=arc_size,data = now[i,] , aes(x = pos1, y = y2, xend = pos2, yend = y2), curvature = now$curv[i],colour=t2tINV_color) 
    }
    idx= c(idx,3)
  }
  SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y2)) + 
    geom_point(data=now,size=.5,aes(x=pos2,y=y2))
  
  now = df[df$SVtype == "h2hINV",]
  if (nrow(now) > 0){
    for (i in 1:nrow(now)){
      SV_plot = SV_plot + geom_curve( size=arc_size,data = now[i,] , 
                                      aes(x = pos1, y = y2, xend = pos2, yend = y2), curvature = -1*now$curv[i],colour=h2hINV_color) 
    }
    idx= c(idx,4)
  }
  SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y2)) + geom_point(data=now,size=.5,aes(x=pos2,y=y2))
  
  
  SV_plot = SV_plot + theme(axis.ticks.x=element_blank(),panel.border = element_blank(),
                            axis.title.y=element_text(colour="white"),
                            axis.text.y=element_text(colour="white"),
                            axis.ticks.y=element_line(colour="white")) + 
    scale_x_continuous(expand = c(0.01,0.01))+ coord_cartesian(xlim=c(min_coord,max_coord))
  
  idx = c(1,2,3,4)
  vals = c(DEL_color=DEL_color,DUP_color=DUP_color,t2tINV_color=t2tINV_color,h2hINV_color=h2hINV_color)
  labs = c('DEL','DUP',"t2tINV","h2hINV")
  
  SV_plot = SV_plot +  scale_colour_manual(name = 'SV type', 
                                           values =c(DEL_color,DUP_color,t2tINV_color,h2hINV_color),
                                           labels = labs[idx]) + theme(legend.position="none")
  
  #-------------------------------------------------------------------------------------------------------
  if (max(CNVsnow$total_cn) <=10){
    CNV_plot = ggplot() +
      geom_segment(data=CNVsnow, aes(x = start, y =total_cn , xend = end, yend = total_cn),size=2) + 
      common_ggplot2 + ylab("CN")+ xlab(NULL)+
      #scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")})#,limits=c(min_coord,max_coord))
      scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")},limits=c(min_coord,max_coord))
    CNV_plot = CNV_plot + 
      scale_y_continuous(minor_breaks = NULL,breaks=c(0,sort(unique(CNVsnow$total_cn))),limits=c(min(CNVsnow$total_cn) -0.35,max(CNVsnow$total_cn)+.35)) 
  }
  
  if (max(CNVsnow$total_cn) <=100 & max(CNVsnow$total_cn) >10){
    idx = which(CNVsnow$total_cn ==0);if(length(idx)>0){CNVsnow$total_cn[which(CNVsnow$total_cn ==0)]=0.99}
    mm2 = max(CNVsnow$total_cn)
    if(mm2 <=20){mm = c(0,1,2,4,10,mm2)}
    if(mm2 <=40 & mm2 >20){mm = c(0,1,2,4,10,20,mm2)}
    if(mm2 <=100 & mm2 >40){mm = c(0,1,2,4,10,20,mm2)}
    
    CNV_plot = ggplot() + geom_segment(data=CNVsnow, aes(x = start, y =log(total_cn,base=2), xend = end, yend = log(total_cn,base=2)),size=2) + 
      common_ggplot2 + ylab("CN")+ xlab(NULL)+ 
      #scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")})#,limits=c(min_coord,max_coord))
      scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")},limits=c(min_coord,max_coord))
    
    CNV_plot = CNV_plot + scale_y_continuous(minor_breaks = NULL,breaks=log(mm,base=2), 
                                             limits=c( -.01 + log(min(CNVsnow$total_cn),base=2), 
                                                       log(max(CNVsnow$total_cn) , base=2)+.01 ), 
                                             labels=as.character(mm))
  }
  
  if (max(CNVsnow$total_cn) >100){
    idx = which(CNVsnow$total_cn ==0);if(length(idx)>0){CNVsnow$total_cn[which(CNVsnow$total_cn ==0)]=0.99}
    CNV_plot = ggplot() +
      geom_segment(data=CNVsnow, aes(x = start, y =log(total_cn,base=10) , xend = end, yend = log(total_cn,base=10)),size=2) + 
      common_ggplot2 + #scale_color_manual(values=c("forestgreen","red1","blue","brown")) +
      ylab("CN")+ xlab(NULL)+ 
      #scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")})#,limits=c(min_coord,max_coord))
      scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")},limits=c(min_coord,max_coord))
    mm2 = max(CNVsnow$total_cn)
    mm = c(0,1,2,4,20,50,100,mm2)
    CNV_plot = CNV_plot + scale_y_continuous(minor_breaks = NULL,breaks=log(mm,base=10), 
                                             limits=c( log(min(CNVsnow$total_cn),base=10), 
                                                       log(max(CNVsnow$total_cn) , base=10)), 
                                             labels=as.character(mm))
  }
  
  
  #-------------------------------------------------------------------------------------------------------
  common_ggplot2_chrom =  theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                plot.background = element_blank(),
                                axis.text=element_blank(),
                                axis.title=element_blank(),
                                plot.title=element_blank(), 
                                axis.ticks=element_blank())
  
  
  ##chr_info = readRDS("/Users/isidro/Dropbox/Park_lab/paper_chromothripsis/chr_info.rds")
  ## Added chrom info for ref genome
  if(genome == "hg19"){
    chr_info <- chr_info
  }else if(genome == "hg38"){
    chr_info <- chr_info_hg38
  }else{
    stop("Reference genome is not supported - Only use hg19 or hg38")
  }
  chr_info$color[chr_info$gieStain == "gneg"] = "white"
  chr_info$color[chr_info$gieStain == "gpos25"] = "grey75"
  chr_info$color[chr_info$gieStain == "gpos50"] = "grey50"
  chr_info$color[chr_info$gieStain == "gpos75"] = "grey25"
  chr_info$color[chr_info$gieStain == "gpos100"] = "grey0"
  chr_info$color[chr_info$gieStain == "acen"] = "red"
  chr_info = chr_info[chr_info$seqnames==chr,]
  chr_info$y = rep(1,nrow(chr_info))
  
  if(nrow(chr_info)<7){chr_info_annot=chr_info}else{
    chr_info_annot=chr_info[seq(3,(nrow(chr_info)-3),3),]}
  
  ideogram =ggplot(data=chr_info,aes(x=y,y=y)) +geom_point(colour="white")
  ideogram = ideogram +  
    geom_rect(data=chr_info,mapping = aes(xmin = start, xmax = end,   ymin = 0, ymax = 1),
              fill = chr_info$color,color="black",size=.1) 
  
  ideogram = ideogram +ylim(0,4)
  ideogram = ideogram + annotate(geom = "text",x = (chr_info_annot$start+chr_info_annot$end)/2,y=chr_info_annot$y + 1.2,
                                 label=chr_info_annot$name,vjust=.5, angle=90,size=2)
  #ideogram = ideogram +theme_bw() + common_ggplot2_chrom + #xlim(min_coord,max_coord) + 
  ideogram = ideogram +theme_bw() + common_ggplot2_chrom + #xlim(min_coord,max_coord) + 
    scale_x_continuous(expand = c(0.01,0.01)) +
    coord_cartesian(xlim=c(min_coord, max_coord))
  
  summary$pos = paste(summary$chrom,":",summary$start,"-",summary$end,sep="")
  
  summary$oscillations = paste(summary$max_number_oscillating_CN_segments_2_states,";",                               summary$max_number_oscillating_CN_segments_3_states,sep="")
  
  cols_sel = c("pos",
               #"clusterSize",
               "clusterSize_including_TRA",
               "number_SVs_sample",
               "max_number_oscillating_CN_segments_3_states",
               "number_CNV_segments",
               "chr_breakpoint_enrichment",
               "pval_exp_cluster",
               "pval_fragment_joins")
  # ,
  #              "inter_other_chroms_coords_all")
  
  table_now = summary[which(summary$chrom == cand),cols_sel]
  
  names(table_now) = c("Position",
                       #"Nb. Intrachr. SVs", 
                       "Total nb. SVs (intrachr. + transl.)",
                       "SVs in sample",
                       "Oscillating CN (2 and 3 states)","CN segments",
                       "Pval chr. breakp. enrich.",
                       "Pval exponential dist. breakpoints",
                       "Pval fragment joins")
  #,
  #                     "Links with other chrs")
  
  mytheme <- gridExtra::ttheme_minimal(padding = unit(c(1.8,1.8),"mm"),
                                       core = list(fg_params=list(cex = .6)),
                                       colhead = list(fg_params=list(cex = .6)),
                                       rowhead = list(fg_params=list(cex = .6)))
  
  for (cc in 1:ncol(table_now)){if( is.numeric(table_now[,cc])){ table_now[,cc] = round(table_now[,cc],digits=2)  }}
  
  table_now = t(table_now) 
  colnames(table_now)= sample_name
  ss <- tableGrob(table_now, theme=mytheme)
  
  
  ########################################
  gp1 <- ggplotGrob(ideogram + theme(plot.margin=unit(c(0.5,0.5,0,0),"cm")))
  gp2 <- ggplotGrob(SV_plot+ theme(legend.position = "none") +
                      theme(plot.margin=unit(c(0.5,0.5,0,0), "cm"),axis.text.x=element_blank())) 
  gp3 <- ggplotGrob(CNV_plot+theme(plot.margin=unit(c(0,0.5,0.2,0), "cm"),legend.position="none"))
  gp4 <- ss
  
  gp3$widths <- gp2$widths
  gp1$widths <- gp2$widths
  
  if (!is.null(BAF)){
    gp5 <- ggplotGrob(BAF_plot + theme(plot.margin=unit(c(0.5,0.5,0.5,0),"cm"),legend.position="none"))
    gp5$widths <- gp2$widths
    return(k=list(gp1,gp2,gp3,gp4,gp5))
  }else{
    return(k=list(gp1,gp2,gp3,gp4))
  }
  #return(k=list(ideogram,SV_plot,CN_plot,BAF_plot))
}


get_chromothripsis <- function(dir, midfix, chr = "1", to_plot = T, ref = "hg38") {
  # dir = "/Users/bl0033/data/SV/test4paper/bfb/nDSB2_Un0.5_frag50_break0_model0_ptype1_pcorr0.5/r5"
  # midfix = "c121_div8"
  # chr = "2"
  # sample_name=""
  # DEL_color='blue1'
  # DUP_color='red'
  # t2tINV_color="forestgreen"
  # h2hINV_color="black"
  # arc_size=.2
  # ref = "hg38"

  fsv <- file.path(dir, paste0("SVData_", midfix, ".tsv"))
  # cannot compute oscillating CN after merging for some cases, but reasonable values for other cases
  fcn <- file.path(dir, paste0("CNData_merged_", midfix, ".tsv"))
  #fcn <- file.path(dir, paste0("CNData_", midfix, ".tsv"))
  fplot <- file.path(dir, paste0("plot_chr", chr, "_", midfix, ".pdf"))
  print(fsv)
  print(fcn)
  SVData_svgen <- read.delim(fsv)
  SVData_svgen = SVData_svgen %>% mutate(svclass = case_when(svclass == "T2TINV" ~ "t2tINV", svclass == "H2HINV" ~ "h2hINV", .default = svclass)) %>% mutate(end2 = if_else(end2 == start1, start1 + 1, end2))
  CNData_svgen <- read.delim(fcn)
  
  # # merge regions to reduce high seg number in shatterseek plot
  # cn_ranges = GRanges(seqnames = CNData_svgen$chromosome, ranges = IRanges(start = CNData_svgen$start, end = CNData_svgen$end), total_cn = CNData_svgen$total_cn, cnA = CNData_svgen$cnA, cnB = CNData_svgen$cnB)
  # #cn_merged = reduce(split(cn_ranges, elementMetadata(cn_ranges)))
  # cn_merged = merge_consecutive_regions(CNData_svgen)
  # length(cn_ranges)
  # length(cn_merged)
  # cn_df = as.data.frame(cn_ranges)
  # cn_df = data.frame(chromosome = cn_merged@seqnames, start = cn_merged@ranges$start, total_cn = cn_merged@metadata)

  # load data
  SV_data <- SVs(
    chrom1 = as.character(SVData_svgen$chrom1),
    pos1 = as.numeric(SVData_svgen$start1),
    chrom2 = as.character(SVData_svgen$chrom2),
    pos2 = as.numeric(SVData_svgen$end2),
    SVtype = as.character(SVData_svgen$svclass),
    strand1 = as.character(SVData_svgen$strand1),
    strand2 = as.character(SVData_svgen$strand2)
  )

  CN_data <- CNVsegs(
    chrom = as.character(CNData_svgen$chromosome),
    start = CNData_svgen$start,
    end = CNData_svgen$end,
    total_cn = CNData_svgen$total_cn
  )

  # run ShatterSeek
  start_time <- Sys.time()
  chromothripsis <- shatterseek(SV.sample = SV_data, seg.sample = CN_data, genome = ref)
  end_time <- Sys.time()
  print(paste0("Running time (s): ", round(end_time - start_time, digits = 2)))

  # Plot chromothripsis per chromosome (not work on ubuntu)
  #plots_chr <- plot_chromothripsis(ShatterSeek_output = chromothripsis, chr = chr, genome = ref)
  plots_chr <- plot_chromothripsis_local(ShatterSeek_output = chromothripsis, chr = chr, genome = ref)
  p <- arrangeGrob(plots_chr[[1]],
    plots_chr[[2]],
    plots_chr[[3]],
    plots_chr[[4]],
    nrow = 4, ncol = 1, heights = c(0.2, .4, .4, .4)
  )
  
  if(to_plot){
    print("saving plot to file")
    pg = plot_grid(p)
    save_plot(fplot, pg, base_height = 6)
    # print("finish saving")
  }

  # print summary
  summ <- chromothripsis@chromSummary
  # names(chromothripsis@detail)
  # find rows with positions
  summ %>% filter(!is.na(start)) -> res
  # res %>%
  #   t() %>%
  #   View()
  fout <- file.path(dir, paste0("shatterseek_", midfix, "_chr", chr, ".tsv"))
  write_tsv(res, fout)
  
  return(chromothripsis)
}


##################### Run on simulated data ############################


option_list <- list(
  make_option(c("-d", "--dir"),
    type = "character", default = "",
    help = "input file directory [default=%default]", metavar = "character"
  ),
  make_option(c("-f", "--fsv"),
    type = "character", default = "",
    help = "The name of the file containing structural variants [default=%default]", metavar = "character"
  ),  
  make_option(c("-n", "--div_ID"),
    type = "character", default = "1",
    help = "The number of cell divisions [default=%default]", metavar = "character"
  ),
  make_option(c("-c", "--cell_ID"),
    type = "character", default = "2",
    help = "The ID of the cell [default=%default]", metavar = "character"
  ),
  make_option(c("-r", "--chr"),
              type = "integer", default = 1,
              help = "The chromosome to examine [default=%default]", metavar = "number"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# dir <- "d:/data/SV/sim/test20220425-2-checked"
# cell_ID <- 3
# div_ID <- 1
#
# dir = "d:/data/SV/real/HX13/nDSB25_nUn0.1_frag49/r2"
# cell_ID = 4
# div_ID = 2
# chr = 22

dir <- opt$dir
fsv <- opt$fsv
cell_ID <- opt$cell_ID
div_ID <- opt$div_ID
chr <- opt$chr

midfix <- paste0("c", cell_ID, "_div", div_ID)
if(fsv != ""){
  print(fsv)
  bname = tools::file_path_sans_ext(fsv)
  fields = str_split(bname, "_")[[1]]
  midfix = paste(fields[2], fields[3], sep = "_")
}
print(midfix)


########## detect chromothripsis
cm <- get_chromothripsis(dir, midfix, chr)



