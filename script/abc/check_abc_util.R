#!/usr/bin/env Rscript

# This script is used to parse the results of ABC on simulated or real data 

library(tidyverse)
library(ggsci)
library(ggpubr)
library(tools)
library(Hmisc)  # for weighted summary, correlation
library(patchwork)
library(GGally)
library(data.table)
library(doParallel)
library(scales)
library(PerformanceAnalytics) # for correlation plots
library(corrplot)
library(modi) # for weighted.var


nbin = 30
mypal = pal_npg("nrc", alpha = 0.7)(10)
model_colours = mypal[5:6] 


cmp_colors = c('#e5f5f9','#99d8c9')
data_colors = c("#5ab4ac", "#d8b365")
#unrep_colors = c('#e5f5f9','#99d8c9','#2ca25f')
unrep_colors = c('#e0f3db','#a8ddb5','#43a2ca')
wgd_colors = c('#e0ecf4','#9ebcda','#8856a7')

lbl_fusion = "mean #fusions per cycle"
lbl_ecdna = "mean #ecDNAs per cell"
lbl_dsb = "DSB rate per cycle"
lbl_unrepair = "%unrepaired DSBs per cycle"
lbl_wgd = "probability of WGD per cell"

# lbl_rate = "mean DSB rate"
# lbl_unrep = "mean fraction of unrepaired DSBs"
# lbl_wgd = "mean probability of WGD"
# lbl_fusion = "mean of average #fusions"


#################### utility functions  #################### 
#from https://stat.ethz.ch/pipermail/r-help/2008-July/168762.html
weighted.var1 <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
                                       na.rm)
}

weighted.var2 <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  (sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2))
}


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# adapted from chart.Correlation
plot_corr_chart <- function(R, histogram = TRUE, nbreak = 8, method = c("pearson", "kendall", 
                                          "spearman"), ...){
  x = checkData(R, method = "matrix")
  if (missing(method)) 
    method = method[1]
  cormeth <- method
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 
                        method = cormeth, cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 0.8/strwidth(txt)
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
                                                                              "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
    text(0.8, 0.8, Signif, cex = cex, col = 2)
  }
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  hist.panel = function(x, ... = NULL) {
    par(new = TRUE)
    # breaks = "FD", 
    hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
         main = "", breaks = nbreak)
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    rug(x)
  }
  if (histogram) 
    pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
          diag.panel = hist.panel)
  else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor)
}


#################### functions for checking ABC rejection predictions #################### 
plot_nDSB <- function(dist_all, infile, ftype = ".png", nrun = 10){
  plt_list = list()
  for(i in 1:nrun){
    # i = 1
    dist_all %>% select(nDSB, dist) %>% arrange(dist) %>% head(topn) -> p1
    
    # on simulated data
    title = paste0("number of real DSBs ", real_nDSB)
    plt = ggplot(p1, aes(x=nDSB)) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr() + geom_vline(xintercept = mean(p1$nDSB), color="blue", linetype="dashed", size=1) + ggtitle(title) + geom_vline(xintercept = real_nDSB, color="red", linetype="dashed", size=1) + xlab("number of DSBs") + scale_x_continuous(limits = c(lmax_val, rmax_val))
    
    fout = str_replace(fname, ".txt", ".jpg")
    ggsave(fout, width = 8, height = 5)
    
    plt_list[[i]] = plt
  }
  
  pnDSB = ggarrange(plotlist = plt_list, nrow = 1)
  #pnDSB = ggarrange(plotlist = plt_list, ncol = 5, nrow = 2)
  fout =  str_replace(infile, ".txt", paste0("_distr", ftype))
  #ggsave(fout, width = 15, height = 6)
  ggsave(fout, pnDSB, width = 2 * length(plt_list), height = 4)
}


plot_hist <- function(dist_all, nrun, real_rDSB, topn){
  plt_list = list()
  # select top n for each run
  for(i in 1:nrun){
    # i = 1
    dist_all %>% select(rDSB, dist=paste0("run",i)) %>% arrange(dist) %>% head(topn) -> p1
    plt = ggplot(p1, aes(x=rDSB)) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + theme_pubr() + geom_vline(xintercept = mean(p1$rDSB), color="blue", linetype="dashed", size=1) + geom_vline(xintercept = real_rDSB, color="red", linetype="dashed", size=1) + xlab("DSB rate (per division)")
    plt_list[[i]] = plt
  }
  
  ggarrange(plotlist = plt_list, ncol = 5, nrow = 2)
  fout =  str_replace(fname, ".txt", "_hist.pdf")
  ggsave(fout, width = 12, height = 6)
}



plot_box <- function(dist_all, nrun, real_rDSB, topn){
  res_sel = data.frame()
  for(i in 1:nrun){
    # i = 1
    dist_all %>% select(rDSB, dist=paste0("run",i)) %>% arrange(dist) %>% head(topn) -> p1
    p1$run = paste0("run",i)
    res_sel = rbind(res_sel, p1)
  }
  res_sel$run = factor(res_sel$run, levels = paste0("run",c(1:nrun)))
  ggplot(res_sel, aes(x=run, y = rDSB)) +
    geom_hline(yintercept = real_rDSB, color="red", linetype="dashed", size=1.5) +
    # geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + theme_pubr() + ylab("estimated DSB rate (per division)") + xlab("")
  
  #+ guides(fill = guide_legend(nrow = 1)) + scale_fill_npg(alpha = 0.7) + labs(fill = "")
  fout =  str_replace(fname, ".txt", "_box.pdf")
  ggsave(fout, width = 12, height = 6)
}



#################### functions for checking simulated data #################### 
get_sv_pos <- function(dir){
  fnames = list.files(dir, pattern = "^SVData")
  sv_all = c()
  sv_full = data.frame()
  for(fb in fnames){
    # fb = fnames[1]
    bname = file_path_sans_ext(fb)
    fields = str_split(bname, pattern = "_")[[1]]
    cell_ID = fields[[2]]

    fsv = file.path(dir, list.files(dir, pattern = paste0("SVData_", cell_ID, "_div.*.tsv")))
    svs = read_tsv(fsv, show_col_types = F)

    svs$cell = cell_ID
    sv_full = rbind(sv_full, svs)

    sv_sel1 = svs %>% unite(pos, chrom1, start1, chrom2, end2, sep="_")
    sv_sel2 = svs %>% unite(pos, chrom2, end2, chrom1, start1, sep="_")
    sv_union = union(sv_sel1$pos, sv_sel2$pos)
    sv_all = c(sv_all, sv_union)
  }
  return(list(sv_all = sv_all, sv_full = sv_full))
}


plot_cn <- function(rdir, pattern = "^CNBin", main = "CN by bin", theme = theme1){
  cn_bin = data.frame()
  fbins = list.files(rdir, pattern = pattern)
  for(fb in fbins){
    # fb = fbins[1]
    bname = file_path_sans_ext(fb)
    fields = str_split(bname, pattern = "_")[[1]]
    sid = fields[[2]]

    fbpath = file.path(rdir, fb)

    cn = read_tsv(fbpath, show_col_types = F)
    cn$sample = sid
    cn_bin = rbind(cn_bin, cn)
  }

  seg_bin = cn_bin %>% select(sample, chrom = "chromosome", start, end, cn = total_cn)
  fout = file.path(rdir, paste0("plot_", str_replace_all(pattern, "[^[:alnum:]]", ""), ".pdf"))

  pc1 = plot.cn.heatmap(seg_bin, main, theme = theme1)
  ggsave(fout, pc1, width = 12, height = 10)
}


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


plot_sv <- function(fb, fsv, fcnv){
  cn = read_tsv(fcnv, show_col_types = F)
  cn$chromosome = paste0("chr", cn$chromosome)
  cnsel = cn %>% select(chr = chromosome, start, end, value1 = cnA, value2 = cnB)
  # cnsel = cnsel %>% mutate(value1 = 1, value2 = 1)   # for normal data

  # CN colors for each haplotype
  MAXCN = 5
  # cols = c("#6283A9","#f0f0f0", "#B9574E", "#3b0107")
  cols = c("#bdd7e7","#f0f0f0", '#fdcc8a','#fc8d59','#e34a33','#b30000')
  col_fun = colorRamp2(breaks = seq(0:MAXCN)-1, colors = cols)

  # output file
  outfile <- file.path(rdir, paste0("plot_", str_replace(fb, ".tsv", ".png")))
  name <- ""
  jpeg(outfile, height=800, width=800)

  circos.initializeWithIdeogram(species="hg19", chromosome.index = paste0("chr", seq(1:22)))

  circos.genomicHeatmap(cnsel, col = col_fun, side = "inside", border = "white")

  # get SV data
  res = get_SV(fsv)
  if(!is.na(res)){
    svs = res$svs
    sv.first = res$sv.first
    sv.second = res$sv.second

    # SV colors
    pa <- c("blue1","red", "black", "forestgreen", "purple")
    svtype <- c("DEL", "DUP", "H2HINV", "T2TINV", "BND")
    names(pa) <- c("DEL", "DUP", "H2HINV", "T2TINV", "BND")
    col_fun_sv = colorRamp2(svtype, pa)
    col <- pa[svs$svclass]

    circos.genomicLink(sv.first, sv.second, col=col)
  }

  circos.clear()

  title(name)
  dev.off()
}


plot_sv_all <- function(rdir, pattern = "^SVData"){
  fnames = list.files(rdir, pattern = pattern)
  for(fb in fnames){
    # fb = fnames[1]
    bname = file_path_sans_ext(fb)
    fields = str_split(bname, pattern = "_")[[1]]
    cell_ID = fields[[2]]

    fsv = file.path(rdir, list.files(rdir, pattern = paste0("SVData_", cell_ID, "_div.*.tsv")))
    fcnv = file.path(rdir, list.files(rdir, pattern = paste0("CNData_", cell_ID, "_div.*.tsv")))
    plot_sv(fb, fsv, fcnv)
  }
}


# plot smc results on simulated data
get_plots_params_by_dsb <- function(sel_rDSB, rfracs, rwgds, smc_all, meanvals){
  plist = list()
  i = 1
  for(sel_frac in rfracs){
    for(sel_wgd in rwgds){
      # if(i == 3) break
      # sel_rDSB = 30
      # sel_frac = 0.3
      # sel_wgd = 0.3
      smc_sel = smc_all %>% filter(real_rDSB == sel_rDSB & real_frac == sel_frac & real_wgd == sel_wgd)
      meanvals_sel = meanvals %>% filter(real_rDSB == sel_rDSB & real_frac == sel_frac & real_wgd == sel_wgd)
      sel_fusion = mean(meanvals_sel$real_fusion)
      sel_ecdna = mean(meanvals_sel$real_ecdna)
      
      print(paste(sel_rDSB, sel_frac, sel_wgd, sel_fusion, sel_ecdna, sep = " "))
      
      smc_sel %>% group_by(real_rDSB, real_frac, real_wgd) %>% tally()
      p = plot_smc_all_sim(meanvals_sel, smc_sel, rsdir, odir, sel_rDSB, sel_frac, sel_wgd, sel_fusion, sel_ecdna, max_ndsb, max_unrepair, max_wgd, max_fusion, max_ecdna)  
      plist[[i]] = p
      i = i + 1
    }
  }
  return(plist)
}


#################### functions for checking ABC SMC predictions #################### 
# read summary statistics of PPC 
get_sstat_all <- function(rdir, ddir, to_plot_sstat = F, to_plot_ppc = F){
  fnames = list.files(rdir, pattern = "smc.*[0-9ab]$")
  sstat_all = data.frame()
  cstat_all = data.frame()
  
  for(i in 1:length(fnames)){
    # f = "smc_nparam3_epsilon0.5_ncell10_maxDSB50_stat24_maxWGD1.0_SA1055"
    f = fnames[i]
    fields = str_split(f, "_")[[1]]
    ncell = str_replace(fields[4], "ncell", "")
    dataset = fields[length(fields)]
    suffix = str_replace(f, "smc", "")
    rstat = check_ppc_stat(rdir, ddir, ncell, dataset, suffix, plot_stat = to_plot_sstat, plot_ppc = to_plot_ppc)
    sstat = rstat$sim_stat
    cstat = rstat$cell_stat
    sstat$dataset = dataset
    cstat$dataset = dataset
    sstat_all = rbind(sstat_all, sstat)
    cstat_all = rbind(cstat_all, cstat)
  }
  
  return(list(sstat_all = sstat_all, cstat_all = cstat_all))
}


plot_param3_sel <- function(res, fres, title, max_ndsb, max_selection, max_unrepair){
  p1 = ggplot(res, aes(x=dsb_rate, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + ggtitle(title) + xlab("rate of DSBs per cycle") + scale_x_continuous(limits = c(0, max_ndsb)) 
  
  p2 = ggplot(res, aes(x=frac_unrepaired, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + ggtitle(title) + xlab(lbl_unrepair) + scale_x_continuous(limits = c(0, max_unrepair)) 
  
  p3 = ggplot(res, aes(x=selection, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + ggtitle(title) + xlab("selection strength") + scale_x_continuous(limits = c(1, max_selection)) 
  
  ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
  fout = paste0(fres, "_hist.pdf")
  ggsave(fout, width = 13, height = 5)
}


plot_param3_frag <- function(res, fres, title, max_ndsb, max_frag, max_unrepair){
  p1 = ggplot(res, aes(x=dsb_rate, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + ggtitle(title) + xlab("rate of DSBs per cycle") + scale_x_continuous(limits = c(0, max_ndsb)) 
  
  p2 = ggplot(res, aes(x=frac_unrepaired, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + ggtitle(title) + xlab(lbl_unrepair) + scale_x_continuous(limits = c(0, max_unrepair)) 
  
  p3 = ggplot(res, aes(x=frag, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + ggtitle(title) + xlab("mean number of DSBs in local fragmentation") + scale_x_continuous(limits = c(0, max_frag)) 
  
  ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
  fout = paste0(fres, "_hist.pdf")
  ggsave(fout, width = 13, height = 5)
}


plot_param3 <- function(res, fres, title, max_ndsb, max_wgd, max_unrepair, real_rDSB = -1, real_frac = -1, real_wgd = -1){
  p1 = ggplot(res, aes(x=dsb_rate, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + ggtitle(title) + xlab("rate of DSBs per cycle") + scale_x_continuous(limits = c(0, max_ndsb)) + geom_vline(xintercept = weighted.mean(res$dsb_rate, res$weight), linetype = "dashed", color = "blue", linewidth = 1.5)
  if(real_rDSB > 0){
    p1 = p1 + geom_vline(xintercept = real_rDSB, linetype = "dashed", color = "red", linewidth = 1.5)
  }
  
  p2 = ggplot(res, aes(x=frac_unrepaired, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + xlab(lbl_unrepair) + scale_x_continuous(limits = c(0, max_unrepair)) + geom_vline(xintercept = weighted.mean(res$frac_unrepaired, res$weight), linetype = "dashed", color = "blue", linewidth = 1.5)  + ggtitle(title)
  if(real_frac > 0){
    p2 = p2 + geom_vline(xintercept = real_frac, linetype = "dashed", color = "red", linewidth = 1.5)
  }
  
  p3 = ggplot(res, aes(x=wgd, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()   + xlab(lbl_wgd) + scale_x_continuous(limits = c(0, max_wgd)) + geom_vline(xintercept = weighted.mean(res$wgd, res$weight), linetype = "dashed", color = "blue", linewidth = 1.5)  + ggtitle(title)
  if(real_wgd > 0){
    p3 = p3 + geom_vline(xintercept = real_wgd, linetype = "dashed", color = "red", linewidth = 1.5)    
  }
  
  p = ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = T)
  fout = paste0(fres, "_hist.pdf")
  ggsave(fout, p, width = 13, height = 5)
  
  return(p)
}


plot_param3_break <- function(res, fres, title, max_ndsb, max_break, max_unrepair){
  p1 = ggplot(res, aes(x=dsb_rate, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + ggtitle(title) + xlab("rate of DSBs per cycle") + scale_x_continuous(limits = c(0, max_ndsb)) + geom_vline(xintercept = weighted.mean(res$dsb_rate, res$weight), linetype = "dashed", color = "red", linewidth = 1.5)
  
  p2 = ggplot(res, aes(x=frac_unrepaired, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + xlab(lbl_unrepair) + scale_x_continuous(limits = c(0, max_unrepair)) + geom_vline(xintercept = weighted.mean(res$frac_unrepaired, res$weight), linetype = "dashed", color = "red", linewidth = 1.5)  + ggtitle(title)
  
  p3 = ggplot(res, aes(x=div_break, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()   + xlab("maximum cell cycle ID with DSBs") + scale_x_continuous(limits = c(0, max_break)) + geom_vline(xintercept = weighted.mean(res$div_break, res$weight), linetype = "dashed", color = "red", linewidth = 1.5)  + ggtitle(title)
  
  p = ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = T)
  fout = paste0(fres, "_hist.pdf")
  ggsave(fout, p, width = 13, height = 5)
  
  return(p)
}



plot_param4 <- function(res, fres, title, max_ndsb, max_wgd, max_unrepair, max_break){
  p1 = ggplot(res, aes(x=dsb_rate, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + ggtitle(title) + xlab("rate of DSBs per cycle") + scale_x_continuous(limits = c(0, max_ndsb)) + geom_vline(xintercept = weighted.mean(res$dsb_rate, res$weight), linetype = "dashed", color = "red", linewidth = 1.5)
  
  p2 = ggplot(res, aes(x=frac_unrepaired, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + xlab(lbl_unrepair) + scale_x_continuous(limits = c(0, max_unrepair)) + geom_vline(xintercept = weighted.mean(res$frac_unrepaired, res$weight), linetype = "dashed", color = "red", linewidth = 1.5)  + ggtitle(title)
  
  p3 = ggplot(res, aes(x=wgd, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()   + xlab(lbl_wgd) + scale_x_continuous(limits = c(0, max_wgd)) + geom_vline(xintercept = weighted.mean(res$wgd, res$weight), linetype = "dashed", color = "red", linewidth = 1.5)  + ggtitle(title)
 
  p4 = ggplot(res, aes(x=div_break, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()   + xlab("maximum cell cycle ID with DSBs") + scale_x_continuous(limits = c(0, max_break)) + geom_vline(xintercept = weighted.mean(res$div_break, res$weight), linetype = "dashed", color = "red", linewidth = 1.5)  + ggtitle(title)
  
  p = ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = T)
  fout = paste0(fres, "_hist.pdf")
  ggsave(fout, p, width = 15, height = 5)
  
  return(p)
}


plot_param2 <- function(res, fres, title, max_ndsb, max_unrepair){
  p1 = ggplot(res, aes(x=dsb_rate, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr() + ggtitle(title) + xlab("rate of DSBs per cycle") + scale_x_continuous(limits = c(0, max_ndsb)) + geom_vline(xintercept = weighted.mean(res$dsb_rate, res$weight), linetype = "dashed", color = "blue", linewidth = 1.5)
  
  p2 = ggplot(res, aes(x=frac_unrepaired, y=after_stat(density), weight = weight))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr() + xlab(lbl_unrepair) + scale_x_continuous(limits = c(0, max_unrepair)) + geom_vline(xintercept = weighted.mean(res$frac_unrepaired, res$weight), linetype = "dashed", color = "blue", linewidth = 1.5)  + ggtitle(title)
  
  p = ggarrange(p1, p2, ncol = 2, nrow = 1)
  fout = paste0(fres, "_hist.pdf")
  ggsave(fout, width = 9, height = 5)
  
  return(p)
}



plot_smc_hist <- function(rdir, pattern){
  fresall = list.files(rdir, pattern = pattern)
  
  # list plot not work, due to lacking of title in subplots?
  plts1 = list()
  plts2 = list()
  k = 0
  j = 0
  for(i in 1:length(fresall)){
    f = fresall[i]
    print(f)
    fields = strsplit(f, "_")[[1]]
    nf = length(fields)
    nparam = as.numeric(str_replace(fields[2], "nparam", ""))
    max_unrepair = as.numeric(str_replace(fields[5], "mp", ""))
    max_ndsb = as.numeric(str_replace(fields[6], "maxDSB", ""))
    max_wgd = as.numeric(str_replace(fields[7], "maxWGD", ""))
    fres = file.path(rdir, f)
    res <- read.table(fres, header = T)
    print(nparam)
    title = fields[nf]
    if(nparam == 3){
      names(res) <- c("dsb_rate", "frac_unrepaired", "wgd", "dist", "weight")
      plt = plot_param3(res, fres, title, max_ndsb, max_wgd, max_unrepair)
      k = k + 1
      plts1[[k]] = plt
    }else{
      names(res) <- c("dsb_rate", "frac_unrepaired", "dist", "weight")
      plt = plot_param2(res, fres, title, max_ndsb, max_unrepair)
      j = j + 1
      plts2[[j]] = plt
    }
  } 
  
  length(plts1)
  length(plts2)
  
  ncol = 3
  nrow = ceiling(length(plts1)/ncol)
  ggarrange(plotlist = plts1, ncol = ncol, nrow = nrow)
  fout = file.path(rdir, "smc_wgd_all_np3.pdf")
  ggsave(fout, width = 20, height = 3 * nrow)
  
  ncol = 3
  nrow = ceiling(j/ncol)
  ggarrange(plotlist = plts2, ncol = ncol, nrow = nrow)
  fout = file.path(rdir, "smc_wgd_all_np2.pdf")
  ggsave(fout, width = 20, height = 3 * nrow)
  
}


get_pop_epsilons <- function(fname){
  fcon = file(fname, "r")
  res = readLines(con = fcon, n = 4)
  close(fcon)
  nsim = strsplit(res[2], ":")[[1]][2]
  tol = res[4]
  fields = strsplit(tol, ":")[[1]][2]
  fields = str_replace(fields, "\\[", "")
  fields = str_replace(fields, "\\]", "")
  tvals = strsplit(fields, ",")[[1]]
  
  return(tvals)
}


plot_smc_all <- function(meanvals, rescmb, rdir, ordered = T, suffix = ""){
  meanvals %>% arrange(desc(dsb_rate)) %>% select(sample) -> ds
  # + scale_fill_manual(values = pal) + guides(fill = guide_legend(nrow = 1)) 
  pdsb = ggplot(rescmb, aes(x=sample, y=dsb_rate, fill=signature_type)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_dsb) + theme_pubr() + theme_pubr() + scale_x_discrete(limits=ds$sample) + coord_flip()  + xlab("") + scale_y_continuous(n.breaks = 10) + guides(fill = guide_legend(nrow = 1)) + scale_fill_npg()
  
  pfrac = ggplot(rescmb, aes(x=sample, y=frac_unrepaired, fill=signature_type)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_unrepair) + theme_pubr()+ theme_pubr() + coord_flip() + xlab("") + scale_y_continuous(n.breaks = 10) + guides(fill = guide_legend(nrow = 1))  + scale_fill_npg()
  if(ordered){
    meanvals %>% arrange(desc(frac_unrepaired)) %>% select(sample) -> df
    pfrac = pfrac + scale_x_discrete(limits=df$sample)
  }else{
    pfrac = pfrac + scale_x_discrete(limits=ds$sample)
  }
  
  
  # reswgd = rescmb %>% filter(wgd >= 0)
  pwgd = ggplot(rescmb, aes(x=sample, y=wgd, fill=signature_type)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_wgd) + theme_pubr() + scale_fill_npg() + theme_pubr() + coord_flip() + xlab("") + scale_y_continuous(n.breaks = 10, limits = c(0,1)) + guides(fill = guide_legend(nrow = 1)) 
  if(ordered){
    meanvals %>% arrange(desc(wgd)) %>% select(sample) -> dw
    pwgd = pwgd + scale_x_discrete(limits=dw$sample) 
  }else{
    pwgd = pwgd + scale_x_discrete(limits=ds$sample) 
  }
  
  ggarrange(pdsb, pfrac, pwgd, nrow = 1, ncol = 3, common.legend = T)
  
  fname = paste0("sum_violin_plot", suffix, ".pdf")
  if(ordered){
    fname = paste0("sum_violin_plot_ordered", suffix, ".pdf")
  }
  fout = file.path(rdir, fname)
  ggsave(fout, width = 15, height = 10)
}

# plots for figure in paper
plot_smc_all_final <- function(meanvals, rescmb, sstat1, cstat1, odir, my_colors = data_colors, nbreak = 5, fsize = 16, w = 14, h = 16){
  meanvals %>% arrange((dsb_rate)) %>% mutate(type=if_else(wgd<0,"darkgrey", "black")) %>% select(sample,type) -> ds
  # + scale_fill_manual(values = pal) + guides(fill = guide_legend(nrow = 1)) 
  pdsb = ggplot(rescmb, aes(x=sample, y=dsb_rate, fill=type)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_dsb) + theme_pubr(base_size = fsize) + scale_x_discrete(limits=ds$sample) + xlab("") + scale_y_continuous(n.breaks = nbreak) + guides(fill = guide_legend(nrow = 1)) + scale_fill_manual(values = my_colors) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + labs(fill = "") + theme(axis.text.x = element_text(colour=ds$type))  
  
  pfrac = ggplot(rescmb, aes(x=sample, y=frac_unrepaired, fill=type)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_unrepair) + theme_pubr(base_size = fsize) + xlab("") + scale_y_continuous(n.breaks = nbreak) + guides(fill = guide_legend(nrow = 1)) + scale_fill_manual(values = my_colors) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + labs(fill = "")
  pfrac = pfrac + scale_x_discrete(limits=ds$sample) + theme(axis.text.x = element_text(colour=ds$type)) 
  
  # reswgd = rescmb %>% filter(wgd >= 0)
  pwgd = ggplot(rescmb, aes(x=sample, y=wgd, fill=type)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_wgd) + theme_pubr(base_size = fsize) + scale_fill_manual(values = my_colors) + xlab("") + scale_y_continuous(n.breaks = nbreak, limits = c(0,1)) + guides(fill = guide_legend(nrow = 1)) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + labs(fill = "")
  pwgd = pwgd + scale_x_discrete(limits=ds$sample) + theme(axis.text.x = element_text(colour=ds$type))  
  
  pfusion = ggplot(sstat1, aes(x=sample, y=avg_fusion, fill=type)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_fusion) + theme_pubr(base_size = fsize) + xlab("") + scale_y_continuous(n.breaks = nbreak) + guides(fill = guide_legend(nrow = 1)) + guides(fill = guide_legend(nrow = 1)) + scale_fill_manual(values = my_colors) + scale_x_discrete(limits=ds$sample) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + labs(fill = "") + theme(axis.text.x = element_text(colour=ds$type)) 
 
   pec = ggplot(cstat1, aes(x=sample, y=nECDNA, fill=type)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_ecdna) + theme_pubr(base_size = fsize)  + xlab("") + scale_y_continuous(n.breaks = nbreak) + guides(fill = guide_legend(nrow = 1)) + guides(fill = guide_legend(nrow = 1)) + scale_fill_manual(values = my_colors) + scale_x_discrete(limits=ds$sample) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + labs(fill = "") + theme(axis.text.x = element_text(colour=ds$type)) 
  
  pres = ggarrange(pdsb, pfrac, pwgd, pfusion, pec, nrow = 5, ncol = 1, common.legend = T)
  
  fname = paste0("sum_real_plot_distr.pdf")
  fout = file.path(odir, fname)
  ggsave(fout, pres, width = w, height = h)
  
  return(pres)
}



plot_smc_all_type <- function(meanvals, rescmb, rdir, ordered = T, suffix = "", my_colors = data_colors){
  meanvals %>% arrange(desc(dsb_rate)) %>% select(sample) -> ds
  # + scale_fill_manual(values = pal) + guides(fill = guide_legend(nrow = 1)) 
  pdsb = ggplot(rescmb, aes(x=sample, y=dsb_rate, fill=type)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_dsb) + theme_pubr() + scale_x_discrete(limits=ds$sample) + coord_flip()  + xlab("") + scale_y_continuous(n.breaks = 10) + guides(fill = guide_legend(nrow = 1)) + scale_fill_manual(values = my_colors)
  
  pfrac = ggplot(rescmb, aes(x=sample, y=frac_unrepaired, fill=type)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_unrepair) + theme_pubr()+ theme_pubr() + coord_flip() + xlab("") + scale_y_continuous(n.breaks = 10) + guides(fill = guide_legend(nrow = 1)) + scale_fill_manual(values = my_colors)
  if(ordered){
    meanvals %>% arrange(desc(frac_unrepaired)) %>% select(sample) -> df
    pfrac = pfrac + scale_x_discrete(limits=df$sample)
  }else{
    pfrac = pfrac + scale_x_discrete(limits=ds$sample)
  }
  
  # reswgd = rescmb %>% filter(wgd >= 0)
  pwgd = ggplot(rescmb, aes(x=sample, y=wgd, fill=type)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_wgd) + theme_pubr() + scale_fill_manual(values = my_colors) + coord_flip() + xlab("") + scale_y_continuous(n.breaks = 10, limits = c(0,1)) + guides(fill = guide_legend(nrow = 1)) 
  if(ordered){
    meanvals %>% arrange(desc(wgd)) %>% select(sample) -> dw
    pwgd = pwgd + scale_x_discrete(limits=dw$sample) 
  }else{
    pwgd = pwgd + scale_x_discrete(limits=ds$sample) 
  }
  
  pres = ggarrange(pdsb, pfrac, pwgd, nrow = 1, ncol = 3, common.legend = T)
  
  fname = paste0("sum_violin_plot", suffix, ".pdf")
  if(ordered){
    fname = paste0("sum_violin_plot_ordered", suffix, ".pdf")
  }
  fout = file.path(rdir, fname)
  ggsave(fout, width = 15, height = 10)
  
  return(pres)
}



plot_smc_all_ctype <- function(meanvals, rescmb, rdir, ordered = T, suffix = ""){
  meanvals %>% arrange(desc(dsb_rate)) %>% select(sample) -> ds
  # + scale_fill_manual(values = pal) + guides(fill = guide_legend(nrow = 1)) 
  pdsb = ggplot(rescmb, aes(x=sample, y=dsb_rate, fill=cell_type)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_dsb) + theme_pubr()+ theme_pubr() + scale_x_discrete(limits=ds$sample) + coord_flip()  + xlab("") + scale_y_continuous(n.breaks = 10) + guides(fill = guide_legend(nrow = 1)) + scale_fill_npg()
  
  pfrac = ggplot(rescmb, aes(x=sample, y=frac_unrepaired, fill=cell_type)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_unrepair) + theme_pubr()+ theme_pubr() + coord_flip() + xlab("") + scale_y_continuous(n.breaks = 10) + guides(fill = guide_legend(nrow = 1))  + scale_fill_npg()
  if(ordered){
    meanvals %>% arrange(desc(frac_unrepaired)) %>% select(sample) -> df
    pfrac = pfrac + scale_x_discrete(limits=df$sample)
  }else{
    pfrac = pfrac + scale_x_discrete(limits=ds$sample)
  }
  
  
  # reswgd = rescmb %>% filter(wgd >= 0)
  pwgd = ggplot(rescmb, aes(x=sample, y=wgd, fill=cell_type)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_wgd) + theme_pubr() + scale_fill_npg() + theme_pubr() + coord_flip() + xlab("") + scale_y_continuous(n.breaks = 10, limits = c(0,1)) + guides(fill = guide_legend(nrow = 1)) 
  if(ordered){
    meanvals %>% arrange(desc(wgd)) %>% select(sample) -> dw
    pwgd = pwgd + scale_x_discrete(limits=dw$sample) 
  }else{
    pwgd = pwgd + scale_x_discrete(limits=ds$sample) 
  }
  
  ggarrange(pdsb, pfrac, pwgd, nrow = 1, ncol = 3, common.legend = T)
  
  fname = paste0("sum_violin_plot", suffix, ".pdf")
  if(ordered){
    fname = paste0("sum_violin_plot_ordered", suffix, ".pdf")
  }
  fout = file.path(rdir, fname)
  ggsave(fout, width = 15, height = 10)
}


plot_smc_all_fit <- function(meanvals, rescmb, rdir, ordered = T, suffix = ""){
  meanvals %>% arrange(desc(dsb_rate)) %>% select(sample) -> ds
  # + scale_fill_manual(values = pal) + guides(fill = guide_legend(nrow = 1)) 
  pdsb = ggplot(rescmb, aes(x=sample, y=dsb_rate, fill=fit)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_dsb) + theme_pubr()+ theme_pubr() + scale_x_discrete(limits=ds$sample) + coord_flip()  + xlab("") + scale_y_continuous(n.breaks = 10) + guides(fill = guide_legend(nrow = 1)) + scale_fill_npg()
  
  pfrac = ggplot(rescmb, aes(x=sample, y=frac_unrepaired, fill=fit)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_unrepair) + theme_pubr()+ theme_pubr() + coord_flip() + xlab("") + scale_y_continuous(n.breaks = 10) + guides(fill = guide_legend(nrow = 1))  + scale_fill_npg()
  if(ordered){
    meanvals %>% arrange(desc(frac_unrepaired)) %>% select(sample) -> df
    pfrac = pfrac + scale_x_discrete(limits=df$sample)
  }else{
    pfrac = pfrac + scale_x_discrete(limits=ds$sample)
  }
  
  
  # reswgd = rescmb %>% filter(wgd >= 0)
  pwgd = ggplot(rescmb, aes(x=sample, y=wgd, fill=fit)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_wgd) + theme_pubr() + scale_fill_npg() + theme_pubr() + coord_flip() + xlab("") + scale_y_continuous(n.breaks = 10, limits = c(0,1)) + guides(fill = guide_legend(nrow = 1)) 
  if(ordered){
    meanvals %>% arrange(desc(wgd)) %>% select(sample) -> dw
    pwgd = pwgd + scale_x_discrete(limits=dw$sample) 
  }else{
    pwgd = pwgd + scale_x_discrete(limits=ds$sample) 
  }
  
  ggarrange(pdsb, pfrac, pwgd, nrow = 1, ncol = 3, common.legend = T)
  
  fname = paste0("sum_violin_plot", suffix, ".pdf")
  if(ordered){
    fname = paste0("sum_violin_plot_ordered", suffix, ".pdf")
  }
  fout = file.path(rdir, fname)
  ggsave(fout, width = 15, height = 10)
}



plot_smc_all_sim <- function(meanvals, rescmb, rdir, odir, real_rDSB, real_frac, real_wgd, real_fusion, real_ecdna, max_rDSB, max_frac, max_wgd, max_fusion, max_ecdna){
  meanvals %>% arrange(desc(dsb_rate)) %>% select(sample) -> ds
  # + scale_fill_manual(values = pal) + guides(fill = guide_legend(nrow = 1)) 
  # + scale_y_continuous(n.breaks = 10)
  pdsb = ggplot(rescmb, aes(x=sample, y=dsb_rate, weight = weight)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_dsb) + theme_pubr() + theme(axis.text.y = element_blank()) + scale_x_discrete(limits=ds$sample) + scale_y_continuous(limits = c(0, max_ndsb)) + coord_flip()  + xlab("")  + guides(fill = guide_legend(nrow = 1)) + geom_hline(yintercept = real_rDSB, linetype="dashed", color = "red")
  
  # meanvals %>% arrange(desc(frac_unrepaired)) %>% select(sample) -> df
  pfrac = ggplot(rescmb, aes(x=sample, y=frac_unrepaired, weight = weight)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_unrepair) + theme_pubr() + theme(axis.text.y = element_blank()) + scale_x_discrete(limits=ds$sample) + scale_y_continuous(limits = c(0, max_unrepair)) + coord_flip() + xlab("") + guides(fill = guide_legend(nrow = 1)) + geom_hline(yintercept = real_frac, linetype="dashed", color = "red")
  
  # meanvals %>% arrange(desc(wgd)) %>% select(sample) -> dw
  reswgd = rescmb %>% filter(wgd >= 0)
  pwgd = ggplot(reswgd, aes(x=sample, y=wgd, weight = weight)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_wgd) + theme_pubr() + theme(axis.text.y = element_blank()) + scale_x_discrete(limits=ds$sample) + scale_y_continuous(limits = c(0, max_wgd)) + coord_flip() + xlab("") + guides(fill = guide_legend(nrow = 1)) + geom_hline(yintercept = real_wgd, linetype="dashed", color = "red")
 
  pfusion = ggplot(rescmb, aes(x=sample, y=avg_fusion)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_fusion) + theme_pubr() + theme(axis.text.y = element_blank()) + scale_x_discrete(limits=ds$sample) + scale_y_continuous(limits = c(0, max_fusion)) + coord_flip() + xlab("") + guides(fill = guide_legend(nrow = 1)) + geom_hline(yintercept = real_fusion, linetype="dashed", color = "red")
  
  pecdna = ggplot(rescmb, aes(x=sample, y=nECDNA)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + ylab(lbl_ecdna) + theme_pubr() + theme(axis.text.y = element_blank()) + scale_x_discrete(limits=ds$sample) + scale_y_continuous(limits = c(0, max_ecdna)) + coord_flip() + xlab("") + guides(fill = guide_legend(nrow = 1)) + geom_hline(yintercept = real_ecdna, linetype="dashed", color = "red")
  
  p = ggarrange(pdsb, pfrac, pwgd, pfusion, pecdna, nrow = 1, ncol = 5, common.legend = T)
  
  suffix = paste(real_rDSB, real_frac, real_wgd, sep="_")
  fout = file.path(odir, paste0("sum_violin_plot_", suffix, ".pdf"))
  ggsave(fout, p, width = 14, height = 3)
  
  return(p)
}


# read and plot SMC results for one set of parameters
check_one_sim_smc <- function(sdir, suffix, title, max_ndsb, max_wgd, max_unrepair, real_rDSB, real_frac, real_wgd, plot = T){
  # read abc smc results
  fres = file.path(sdir, paste0("smc", suffix))
  res <- read.table(fres, header = T)
  names(res) <- c("dsb_rate", "frac_unrepaired", "wgd", "dist", "weight")
  
  if(plot){
    plot_param3(res, fres, title, max_ndsb, max_wgd, max_unrepair, real_rDSB, real_frac, real_wgd)
  }
  
  return(res)
}


# check predictive posterior distributions of summary statistics
check_one_sim_stat <- function(sdir, ncell, suffix, midfix, header_ppc1, title, nbin = 30, nrow = 6, width = 12, height = 18){
  # read PPC summary
  fppc = file.path(sdir, paste0("resample"), paste0("check-smc", suffix))
  rpp = read_tsv(fppc, col_names = F, show_col_types = F)
  names(rpp) = header_ppc1
  # rpg = rpp %>% gather(all_of(header), key = "stat", value = "val")
  
  # read real summary stats
  freal = file.path(sdir, paste0("target_ncell", midfix))
  sreal = read.table(freal) %>% t() %>% as.data.frame()
  names(sreal) = header_ppc1
  
  check_summary_stats_boxplot(rpp, sreal, ncell, fppc, title, nrow, width, height)
  check_summary_stats_hist(rpp, sreal, ncell, fppc, title, nbin, nrow, width, height)   
}
 
 
get_one_sim_all <- function(dir, ncell, max_ndsb, max_unrepair, max_wgd, epsilon, header_ppc1, odir, check_ppc, plot, nrow = 6, width = 12, height = 18){
  bname = basename(dir)
  
  # get parameters for each run
  seed = str_replace(bname, "seed", "")
  title = seed
  freal = list.files(dir, pattern = "target_")
  fields = str_split(freal, pattern = "_")[[1]]
  real_rDSB = as.numeric(str_replace(fields[3], "dsb", "")) 
  real_frac = as.numeric(str_replace(fields[4], "frac", ""))
  real_wgd = as.numeric(str_replace(fields[5], "wgd", ""))
  
  # read real summary statistics
  tdir = file.path(dir, "target")
  tstat_sim = get_sstat_sim(tdir)
  tstat_cell = get_sstat_cell(tdir)
  real_fusion = tstat_sim$avg_fusion
  real_ecdna = tstat_cell$nECDNA
  tstat = list(ss = tstat_sim, cs = tstat_cell)
  
  midfix = paste0(ncell, "_dsb", real_rDSB, "_frac", real_frac, "_wgd", real_wgd)
  suffix = paste0("_sim_ncell", midfix, "_nparam3_epsilon", epsilon, "_maxDSB", max_ndsb, "_maxWGD1_2340225")
  
  # using seed as figure title, need res_smc to get a full result table with PPC data
  res_smc0 = check_one_sim_smc(dir, suffix, title, max_ndsb, max_wgd, max_unrepair, real_rDSB, real_frac, real_wgd, plot)
  
  if(check_ppc){
    title = paste(real_rDSB, real_frac, real_wgd, seed)
    
    # nrow = 6 
    # width = 14 
    # height = 20
    
    check_one_sim_stat(dir, ncell, suffix, midfix, header_ppc1, title, nbin, nrow, width, height)
    
    # very slow to read 10 files for 500 runs
    res_stat = get_ppc_sstat(dir)
    if(plot){
      plot_ppc_hist(res_stat, dir, suffix, title, tstat)
      plot_ppc_boxplot(res_stat, dir, suffix, title, tstat)     
    }    
  }
  
  df_real = data.table(sample = seed, real_rDSB = real_rDSB, real_frac = real_frac, real_wgd = real_wgd, real_fusion = real_fusion, real_ecdna = real_ecdna)
  
  # get merged results
  sim_stat = res_stat$sim_stat
  cell_stat = res_stat$cell_stat %>% dplyr::rename(cell_nDSB = nDSB, cell_nUnrepair = nUnrepair)
  #  not the same value: c("nDSB", "nUnrepair")
  sstat = cbind(sim_stat, cell_stat)
  res_smc1 = cbind(res_smc0, df_real)
  res_smc2 = cbind(res_smc1, sstat)
  # smc_all = rbind(smc_all, res_smc)
  
  # get mean value for each run
  # mval = data.table(dsb_rate = weighted.mean(dsb_rate, weight), frac_unrepaired = weighted.mean(res_smc$frac_unrepaired, weight), wgd = weighted.mean(res_smc$wgd, weight), fusion = mean(sim_stat$avg_fusion), ecdna = mean(cell_stat$nECDNA))
  # mval = cbind(mval, df_real)
  # meanvals = rbind(meanvals, mval)
  
  return(res_smc2)
}


get_smc_all_real <- function(rdir, max_wgd = 1, max_unrepair = 1, plot = F){
  pattern = "smc.*[a|b|0-9]$"
  fresall = list.files(rdir, pattern = pattern)
  
  meanvals = data.frame()
  resall = data.frame()
  for(i in 1:length(fresall)){
    print(i)
    f = fresall[i]
    print(f)
    fields = strsplit(f, "_")[[1]]
    nf = length(fields)
    nparam = as.numeric(str_replace(fields[2], "nparam", ""))
    # max_unrepair = as.numeric(str_replace(fields[5], "mp", ""))
    max_ndsb = as.numeric(str_replace(fields[5], "maxDSB", ""))
    # max_wgd = as.numeric(str_replace(fields[7], "maxWGD", ""))
    sample = fields[nf]
    
    fres = file.path(rdir, f)
    res <- read.table(fres, header = T)
    
    if(nparam == 3){
      names(res) <- c("dsb_rate", "frac_unrepaired", "wgd", "dist", "weight")
      if(plot) plot_param3(res, fres, sample, max_ndsb, max_wgd, max_unrepair)
    }else{
      names(res) <- c("dsb_rate", "frac_unrepaired", "dist", "weight")
      res$wgd = -1
      if(plot) plot_param2(res, fres, sample, max_ndsb, max_unrepair)
    }
    
    res$max_ndsb = max_ndsb
    res$nparam = nparam
    res$sample = sample
    ns = nrow(res)
    # print(ns)
    
    mval = data.frame(sample = sample, dsb_rate = weighted.mean(res$dsb_rate, res$weight), frac_unrepaired = weighted.mean(res$frac_unrepaired, res$weight), wgd = weighted.mean(res$wgd, res$weight), wsd_dsb = sqrt(wtd.var(res$dsb_rate, res$weight, normwt = T)), wsd_unrepair = sqrt(wtd.var(res$frac_unrepaired, res$weight, normwt = T)), wsd_wgd = sqrt(wtd.var(res$wgd, res$weight, normwt = T)))
    meanvals = rbind(meanvals, mval)
    
    resall = rbind(resall, res)
  }
  
  return(list(resall = resall, meanvals = meanvals))
}


# show points and error bar
plot_corr_smc <- function(meancmb, rdir, suffix = "_pdo"){
  p3 = ggplot(meancmb, aes(x=frac_unrepaired, y = dsb_rate, color = signature_type, label = sample)) + geom_point(size = 2) + theme_pubr()  + guides(color = guide_legend(nrow = 1)) + geom_text(hjust=0, vjust=0) + ylab("mean DSB rate") + xlab("mean fraction of unrepaired DSBs") + scale_color_npg()
  
  p1 = ggplot(meancmb, aes(x=signature_type, y = dsb_rate, color = signature_type)) + geom_point(size = 2) + theme_pubr()  + guides(color = guide_legend(nrow = 1)) + geom_errorbar(aes(ymin = dsb_rate - wse_dsb, ymax = dsb_rate + wse_dsb)) + scale_color_npg()
  
  p2 = ggplot(meancmb, aes(x=signature_type, y = frac_unrepaired, color = signature_type)) + geom_point(size = 2) + theme_pubr()  + guides(color = guide_legend(nrow = 1)) + geom_errorbar(aes(ymin = frac_unrepaired - wse_unrepair, ymax = frac_unrepaired + wse_unrepair)) + scale_color_npg()
  
  # ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = T)
  (p1 | p2) / p3
  # + plot_layout(guides = "collect") & theme(plot.tag = element_text(family = "EB Garamond", size = 8), legend.position = "top")
  fout = file.path(rdir, paste0("corr_signature", suffix, ".pdf"))
  ggsave(fout, width = 15, height = 10)
}


# show box plot
plot_corr_smc_box <- function(meancmb, rdir, suffix = "_pdo"){
  p3 = ggplot(meancmb, aes(x=frac_unrepaired, y = dsb_rate, color = signature_type, label = sample)) + geom_point(size = 2) + theme_pubr()  + guides(color = guide_legend(nrow = 1)) + geom_text(hjust=0, vjust=0) + ylab("mean DSB rate") + xlab("mean fraction of unrepaired DSBs") + scale_color_npg()
  
  p1 = ggplot(meancmb, aes(x=signature_type, y = dsb_rate, color = signature_type)) + geom_boxplot() + theme_pubr()  + guides(color = guide_legend(nrow = 1)) + scale_color_npg()
  
  p2 = ggplot(meancmb, aes(x=signature_type, y = frac_unrepaired, color = signature_type)) + geom_boxplot() + theme_pubr()  + guides(color = guide_legend(nrow = 1)) + scale_color_npg()
  
  # ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = T)
  (p1 | p2) / p3
  # + plot_layout(guides = "collect") & theme(plot.tag = element_text(family = "EB Garamond", size = 8), legend.position = "top")
  fout = file.path(rdir, paste0("corr_signature", suffix, ".pdf"))
  ggsave(fout, width = 15, height = 10)
}


# show box plot vs cell type
plot_corr_smc_ctype <- function(meancmb, rdir, suffix = "_pdo"){
  p3 = ggplot(meancmb, aes(x=frac_unrepaired, y = dsb_rate, color = cell_type, label = sample)) + geom_point(size = 2) + theme_pubr()  + guides(color = guide_legend(nrow = 1)) + geom_text(hjust=0, vjust=0) + ylab("mean DSB rate") + xlab("mean fraction of unrepaired DSBs") + scale_color_npg()
  
  p1 = ggplot(meancmb, aes(x=cell_type, y = dsb_rate, color = cell_type)) + geom_boxplot() + theme_pubr()  + guides(color = guide_legend(nrow = 1)) + scale_color_npg()
  
  p2 = ggplot(meancmb, aes(x=cell_type, y = frac_unrepaired, color = cell_type)) + geom_boxplot() + theme_pubr()  + guides(color = guide_legend(nrow = 1)) + scale_color_npg()
  
  # ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = T)
  (p1 | p2) / p3
  # + plot_layout(guides = "collect") & theme(plot.tag = element_text(family = "EB Garamond", size = 8), legend.position = "top")
  fout = file.path(rdir, paste0("corr_ctype", suffix, ".pdf"))
  ggsave(fout, width = 15, height = 10)
}


plot_diff_by_group <- function(meanvals, odir, y1, y2, y3, y4, y5, base_size = 12, ncol = 1, width = 7, height = 23){
  p1 = ggplot(meanvals, aes(x = as.factor(real_rDSB), y = diff_dsb)) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y1) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_wgd)
  p2 = ggplot(meanvals, aes(x = as.factor(real_rDSB), y = diff_frac)) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y2) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_wgd)
  p3 = ggplot(meanvals, aes(x = as.factor(real_rDSB), y = diff_wgd)) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y3) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_wgd)
  p4 = ggplot(meanvals, aes(x = as.factor(real_rDSB), y = diff_fusion)) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y4) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_wgd)
  p5 = ggplot(meanvals, aes(x = as.factor(real_rDSB), y = diff_ecdna)) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y5) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_wgd)
  ggarrange(p1, p2, p3, p4, p5, nrow = ceiling(5 / ncol), ncol = ncol)
  fout = file.path(odir, "sim_abc_diff_group.pdf")
  ggsave(fout, width = width, height = height) 
}


plot_diff_by_wgd <- function(meanvals, y1, y2, y3, y4, y5, base_size = 12, ncol = 1){
  p1 = ggplot(meanvals, aes(x = as.factor(real_frac), y = diff_dsb, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size) + xlab("") + ylab(y1) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_rDSB ~ real_wgd) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p2 = ggplot(meanvals, aes(x = as.factor(real_frac), y = diff_frac, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y2) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_rDSB ~ real_wgd) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p3 = ggplot(meanvals, aes(x = as.factor(real_frac), y = diff_wgd, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y3) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_rDSB ~ real_wgd) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p4 = ggplot(meanvals, aes(x = as.factor(real_frac), y = diff_fusion, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y4) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_rDSB ~ real_wgd) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p5 = ggplot(meanvals, aes(x = as.factor(real_frac), y = diff_ecdna, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y5) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_rDSB ~ real_wgd) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  pa = ggarrange(p1, p2, p3, p4, p5, nrow = ceiling(5 / ncol), ncol = ncol, common.legend = T)
  # fout = file.path(odir, "sim_abc_diff_wgd.pdf")
  # ggsave(fout, width = width, height = height) 
  return(pa)
}


plot_diff_by_wgd_long <- function(meanvals, y1, y2, y3, y4, y5, base_size = 12, ncol = 1){
  p1 = ggplot(meanvals, aes(x = as.factor(real_frac), y = diff_dsb, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size) + xlab("") + ylab(y1) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_wgd ~ real_rDSB) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p2 = ggplot(meanvals, aes(x = as.factor(real_frac), y = diff_frac, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y2) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_wgd ~ real_rDSB) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p3 = ggplot(meanvals, aes(x = as.factor(real_frac), y = diff_wgd, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y3) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_wgd ~ real_rDSB) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p4 = ggplot(meanvals, aes(x = as.factor(real_frac), y = diff_fusion, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y4) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_wgd ~ real_rDSB) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p5 = ggplot(meanvals, aes(x = as.factor(real_frac), y = diff_ecdna, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y5) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_wgd ~ real_rDSB) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  pa = ggarrange(p1, p2, p3, p4, p5, nrow = ceiling(5 / ncol), ncol = ncol, common.legend = T)
  # fout = file.path(odir, "sim_abc_diff_wgd.pdf")
  # ggsave(fout, width = width, height = height) 
  return(pa)
}


plot_diff_by_unrepair <- function(meanvals, y1, y2, y3, y4, y5, base_size = 12, ncol = 1){
  p1 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = diff_dsb, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y1) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_rDSB ~ real_frac) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd) + theme(legend.position = "none", axis.text.x = element_blank())
  p2 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = diff_frac, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y2) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_rDSB ~ real_frac) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd) + theme(legend.position = "none", axis.text.x = element_blank())
  p3 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = diff_wgd, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y3) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_rDSB ~ real_frac) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd) + theme(legend.position = "none", axis.text.x = element_blank())
  p4 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = diff_fusion, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y4) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_rDSB ~ real_frac) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd) + theme(legend.position = "none", axis.text.x = element_blank())
  p5 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = diff_ecdna, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y5) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_rDSB ~ real_frac) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd) + theme(legend.position = "none", axis.text.x = element_blank())
  pa = ggarrange(p1, p2, p3, p4, p5, nrow = ceiling(5 / ncol), ncol = ncol, common.legend = T)
  return(pa)
}


plot_diff_by_unrepair_long <- function(meanvals, y1, y2, y3, y4, y5, base_size = 12, ncol = 1){
  p1 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = diff_dsb, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y1) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_rDSB) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd) + theme(legend.position = "none", axis.text.x = element_blank())
  p2 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = diff_frac, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y2) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_rDSB) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd) + theme(legend.position = "none", axis.text.x = element_blank())
  p3 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = diff_wgd, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y3) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_rDSB) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd) + theme(legend.position = "none", axis.text.x = element_blank())
  p4 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = diff_fusion, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y4) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_rDSB) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd) + theme(legend.position = "none", axis.text.x = element_blank())
  p5 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = diff_ecdna, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y5) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_rDSB) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd) + theme(legend.position = "none", axis.text.x = element_blank())
  pa = ggarrange(p1, p2, p3, p4, p5, nrow = ceiling(5 / ncol), ncol = ncol, common.legend = T)
  return(pa)
}




plot_diff_prop_by_group <- function(meanvals, y1, y2, y3, y4, y5, base_size = 12){
  p1 = ggplot(meanvals, aes(x = as.factor(real_rDSB), y = prop_diff_dsb)) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y1) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_wgd)
  p2 = ggplot(meanvals, aes(x = as.factor(real_rDSB), y = prop_diff_frac)) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y2) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_wgd)
  p3 = ggplot(meanvals, aes(x = as.factor(real_rDSB), y = prop_diff_wgd)) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y3) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_wgd)
  p4 = ggplot(meanvals, aes(x = as.factor(real_rDSB), y = prop_diff_fusion)) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y4) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_wgd)
  p5 = ggplot(meanvals, aes(x = as.factor(real_rDSB), y = prop_diff_ecdna)) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y5) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_wgd)
  pa = ggarrange(p1, p2, p3, p4, p5, nrow = 1, ncol = 5, common.legend = T)
  return(pa)
}



plot_diff_prop_by_unrepair <- function(meanvals, y1, y2, y3, y4, y5, base_size = 12){
  p1 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = prop_diff_dsb, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y1) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_rDSB) + theme(legend.position = "none", axis.text.x = element_blank()) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd)
  p2 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = prop_diff_frac, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y2) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_rDSB) + theme(legend.position = "none", axis.text.x = element_blank()) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd)
  p3 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = prop_diff_wgd, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y3) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_rDSB) + theme(legend.position = "none", axis.text.x = element_blank()) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd)
  p4 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = prop_diff_fusion, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y4) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_rDSB) + theme(legend.position = "none", axis.text.x = element_blank()) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd)
  p5 = ggplot(meanvals, aes(x = as.factor(real_wgd), y = prop_diff_ecdna, fill = as.factor(real_wgd))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y5) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_frac ~ real_rDSB) + theme(legend.position = "none", axis.text.x = element_blank()) + scale_fill_manual(values = wgd_colors) + labs(fill = lbl_wgd)
  pa = ggarrange(p1, p2, p3, p4, p5, nrow = 1, ncol = 5, common.legend = T)
  return(pa)
}



plot_diff_prop_by_wgd <- function(meanvals, y1, y2, y3, y4, y5, base_size = 12, ncol = 1){
  p1 = ggplot(meanvals, aes(x = as.factor(real_frac), y = prop_diff_dsb, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size) + xlab("") + ylab(y1) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_wgd ~ real_rDSB) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p2 = ggplot(meanvals, aes(x = as.factor(real_frac), y = prop_diff_frac, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y2) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_wgd ~ real_rDSB) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p3 = ggplot(meanvals, aes(x = as.factor(real_frac), y = prop_diff_wgd, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y3) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_wgd ~ real_rDSB) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p4 = ggplot(meanvals, aes(x = as.factor(real_frac), y = prop_diff_fusion, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y4) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_wgd ~ real_rDSB) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  p5 = ggplot(meanvals, aes(x = as.factor(real_frac), y = prop_diff_ecdna, fill = as.factor(real_frac))) + geom_boxplot() + theme_pubr(base_size = base_size)  + xlab("") + ylab(y5) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey") + facet_grid(real_wgd ~ real_rDSB) + scale_fill_manual(values = unrep_colors) + labs(fill = lbl_unrepair) + theme(legend.position = "none", axis.text.x = element_blank())
  pa = ggarrange(p1, p2, p3, p4, p5, nrow = 1, ncol = 5, common.legend = T)
  # fout = file.path(odir, "sim_abc_diff_wgd.pdf")
  # ggsave(fout, width = width, height = height) 
  return(pa)
}



plot_diff <- function(meanvals, odir, y1, y2, y3){
  mean_plt = meanvals %>% unite(real, real_rDSB, real_frac, real_wgd, sep = ",")
  p1 = ggplot(mean_plt, aes(x = real, y = diff_frac)) + geom_boxplot() + theme_pubr()  + xlab("") + ylab(y1) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey")
  p2 = ggplot(mean_plt, aes(x = real, y = diff_frac)) + geom_boxplot() + theme_pubr()  + xlab("") + ylab(y2) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey")
  p3 = ggplot(mean_plt, aes(x = real, y = diff_wgd)) + geom_boxplot() + theme_pubr()  + xlab("") + ylab(y3) + geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey")
  ggarrange(p1, p2, p3, nrow = 3, ncol = 1)
  fout = file.path(odir, "sim_abc_diff.pdf")
  ggsave(fout, width = 7, height = 15)
}


#################### functions for checking summary statistics #################### 

check_summary_stats_boxplot <- function(rpp, sreal, ncell, fppc, title_all = "", nrow = 6, width = 12, height = 18){
  title = ""
  # box plot to compare with target
  p1 = ggplot(rpp, aes(y=pga)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("percentage genome altered") + geom_hline(yintercept = sreal$pga, linetype="dashed", color = "red")
  p2 = ggplot(rpp, aes(y=mean_div)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("mean divergence") + geom_hline(yintercept = sreal$mean_div, linetype="dashed", color = "red")
  p3 = ggplot(rpp, aes(y=sd_div)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("standard deviation of divergence") + geom_hline(yintercept = sreal$sd_div, linetype="dashed", color = "red")
  
  p4 = ggplot(rpp, aes(y=pgaA)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("percentage genome altered (A)") + geom_hline(yintercept = sreal$pgaA, linetype="dashed", color = "red")
  p5 = ggplot(rpp, aes(y=mean_divA)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("mean divergence (A)") + geom_hline(yintercept = sreal$mean_divA, linetype="dashed", color = "red")
  p6 = ggplot(rpp, aes(y=sd_divA)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("standard deviation of divergence (A)") + geom_hline(yintercept = sreal$sd_divA, linetype="dashed", color = "red")
  
  p7 = ggplot(rpp, aes(y=pgaB)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("percentage genome altered (B)") + geom_hline(yintercept = sreal$pgaB, linetype="dashed", color = "red")
  p8 = ggplot(rpp, aes(y=mean_divB)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("mean divergence (B)") + geom_hline(yintercept = sreal$mean_divB, linetype="dashed", color = "red")
  p9 = ggplot(rpp, aes(y=sd_divB)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("standard deviation of divergence (B)") + geom_hline(yintercept = sreal$sd_divB, linetype="dashed", color = "red")
  
  list1 = list(p1, p2, p3, p4, p5, p6, p7, p8, p9)
  
  for(i in 1:ncell){
    freq = i
    lbl_freq0 = paste0("frequency of breakpoints in ", freq)
    if(i == 1){
      lbl_freq = paste0(lbl_freq0, " cell")
    }else{
      lbl_freq = paste0(lbl_freq0, " cells")
    }
    idx = 9 + i
    hname = paste0("bp", i)
    p = ggplot(rpp, aes(y=hname)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab(lbl_freq) + geom_hline(yintercept = sreal[, idx], linetype="dashed", color = "red")
    list1[[length(list1) + 1]] = p
  }
  
  p10 = ggplot(rpp, aes(y=frac_wgd)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("fraction of cells with WGD") + geom_hline(yintercept = sreal$frac_wgd, linetype="dashed", color = "red")
  
  p11 = ggplot(rpp, aes(y=frac_del)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("fraction of deletions") + geom_hline(yintercept = sreal$frac_del, linetype="dashed", color = "red")
  
  p12 = ggplot(rpp, aes(y=frac_dup)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("fraction of duplications") + geom_hline(yintercept = sreal$frac_dup, linetype="dashed", color = "red")
  
  p13 = ggplot(rpp, aes(y=frac_inv)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("fraction of inversions") + geom_hline(yintercept = sreal$frac_inv, linetype="dashed", color = "red")
  
  p14 = ggplot(rpp, aes(y=frac_tra)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("fraction of intra-chromosomal SVs") + geom_hline(yintercept = sreal$frac_tra, linetype="dashed", color = "red")
  
  list2 = list(p10, p11, p12, p13, p14)
  plt_list = c(list1, list2)
  
  np = length(plt_list)
  ncol = ceiling(np / nrow)
  plt = ggarrange(plotlist = plt_list, nrow = nrow, ncol = ncol, common.legend = T)
  annotate_figure(plt, top = text_grob(title_all, face = "bold"))
  fout = paste0(fppc, "_boxplot.pdf")
  ggsave(fout, width = width, height = width)
}


check_summary_stats_hist <- function(rpp, sreal, ncell, fppc, title_all = "", nbin = 30, nrow = 6, width = 12, height = 18){
  title = ""
  # box plot to compare with target
  p1 = ggplot(rpp, aes(x=pga, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title) + ylab("percentage genome altered") + geom_vline(xintercept = sreal$pga, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$pga), linetype = "dashed", color = "blue", linewidth = 1.5)
  p2 = ggplot(rpp, aes(x=mean_div, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("mean divergence") + geom_vline(xintercept = sreal$mean_div, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$mean_div), linetype = "dashed", color = "blue", linewidth = 1.5)
  p3 = ggplot(rpp, aes(x=sd_div, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("standard deviation of divergence") + geom_vline(xintercept = sreal$sd_div, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$sd_div), linetype = "dashed", color = "blue", linewidth = 1.5)
  
  p4 = ggplot(rpp, aes(x=pgaA, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("percentage genome altered (A)") + geom_vline(xintercept = sreal$pgaA, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$pgaA), linetype = "dashed", color = "blue", linewidth = 1.5)
  p5 = ggplot(rpp, aes(x=mean_divA, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("mean divergence (A)") + geom_vline(xintercept = sreal$mean_divA, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$mean_divA), linetype = "dashed", color = "blue", linewidth = 1.5)
  p6 = ggplot(rpp, aes(x=sd_divA, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("standard deviation of divergence (A)") + geom_vline(xintercept = sreal$sd_divA, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$sd_divA), linetype = "dashed", color = "blue", linewidth = 1.5)
  
  p7 = ggplot(rpp, aes(x=pgaB, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("percentage genome altered (B)") + geom_vline(xintercept = sreal$pgaB, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$pgaB), linetype = "dashed", color = "blue", linewidth = 1.5)
  p8 = ggplot(rpp, aes(x=mean_divB, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("mean divergence (B)") + geom_vline(xintercept = sreal$mean_divB, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$mean_divB), linetype = "dashed", color = "blue", linewidth = 1.5)
  p9 = ggplot(rpp, aes(x=sd_divB, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("standard deviation of divergence (B)") + geom_vline(xintercept = sreal$sd_divB, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$sd_divB), linetype = "dashed", color = "blue", linewidth = 1.5)
  
  list1 = list(p1, p2, p3, p4, p5, p6, p7, p8, p9)
  
  for(i in 1:ncell){
    freq = i
    lbl_freq0 = paste0("frequency of breakpoints in ", freq)
    if(i == 1){
      lbl_freq = paste0(lbl_freq0, " cell")
    }else{
      lbl_freq = paste0(lbl_freq0, " cells")
    }
    idx = 9 + i
    hname = paste0("bp", i)
    p = ggplot(rpp, aes_string(x=hname), aes(y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab(lbl_freq) + geom_vline(xintercept = sreal[, idx], linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(unlist(rpp[, idx])), linetype = "dashed", color = "blue", linewidth = 1.5)
    # print(p)
    list1[[length(list1) + 1]] = p
  }
  
  p10 = ggplot(rpp, aes(x=frac_wgd, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("fraction of cells with WGD") + geom_vline(xintercept = sreal$frac_wgd, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$frac_wgd), linetype = "dashed", color = "blue", linewidth = 1.5)
  
  p11 = ggplot(rpp, aes(x=frac_del, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("fraction of deletions") + geom_vline(xintercept = sreal$frac_del, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$frac_del), linetype = "dashed", color = "blue", linewidth = 1.5)
  
  p12 = ggplot(rpp, aes(x=frac_dup, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("fraction of duplications") + geom_vline(xintercept = sreal$frac_dup, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$frac_dup), linetype = "dashed", color = "blue", linewidth = 1.5)
  
  p13 = ggplot(rpp, aes(x=frac_inv, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title)  + ylab("fraction of inversions") + geom_vline(xintercept = sreal$frac_inv, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$frac_inv), linetype = "dashed", color = "blue", linewidth = 1.5)
  
  p14 = ggplot(rpp, aes(x=frac_tra, y=after_stat(density))) + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue")  + xlab("") + theme_pubr() + ggtitle(title) + ylab("fraction of intra-chromosomal SVs") + geom_vline(xintercept = sreal$frac_tra, linetype="dashed", color = "red", linewidth = 1.5) + geom_vline(xintercept = mean(rpp$frac_tra), linetype = "dashed", color = "blue", linewidth = 1.5)
  
  list2 = list(p10, p11, p12, p13, p14)
  plt_list = c(list1, list2)
  
  np = length(plt_list)
  ncol = ceiling(np / nrow)
  plt = ggarrange(plotlist = plt_list, nrow = nrow, ncol = ncol, common.legend = T)
  annotate_figure(plt, top = text_grob(title_all, face = "bold"))
  fout = paste0(fppc, "_hist.pdf")
  ggsave(fout, width = width, height = height)
}


# read summary statistics from one simulation in a directory
get_sstat_cell <- function(rdir){
  fnames = list.files(rdir, full.names = T, pattern = "sumStats_total_*")
  file_content <- fnames |> map(\(x) read_tsv(x, show_col_types = F))
  # df = do.call(rbind, file_content)
  cs = bind_rows(file_content)
  cs_mean = cs %>% select(-c("cycleID", "cellID")) %>% summarise_all(mean)
  
  # return(list(ss = ss, cs = cs_mean))
  return(cs_mean)
}


get_sstat_sim <- function(rdir){
  fss = file.path(rdir, "sumStats_sim.tsv")
  ss = read_tsv(fss, show_col_types = F)
  ss$avg_fusion = ss$nTelofusion / (ss$nCell - 1)
  
  return(ss)
}


get_ppc_sstat <- function(rdir, osuffix = ""){
  dir_smc = file.path(rdir, paste0("resample", osuffix))
  dir_resample = list.dirs(dir_smc, recursive = F)
  # length(dir_resample)

  cs_all <- map(dir_resample, get_sstat_cell)
  cell_stat <- bind_rows(cs_all)
  
  ss_all <- map(dir_resample, get_sstat_sim)
  sim_stat <- bind_rows(ss_all)
  
  return(list(sim_stat = sim_stat, cell_stat = cell_stat))
}


plot_ppc_boxplot <- function(res_stat, sdir, suffix, title_all = "", tstat = list()){
  sim_stat = res_stat$sim_stat
  cell_stat = res_stat$cell_stat
  title = ""
  pt = ggplot(sim_stat, aes(y=nTelofusion / (nCell - 1))) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab(lbl_fusion) 
  pc = ggplot(cell_stat, aes(y=nECDNA)) + geom_boxplot(width = 0.05)  + xlab("") + theme_pubr() + ggtitle(title) + theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab(lbl_ecdna) 
  
  if(length(tstat) > 0){  
    ncell = unique(sim_stat$nCell)
    pt = pt + geom_hline(yintercept = tstat$ss$nTelofusion / (ncell - 1), linetype = "dashed", color = "red", linewidth = 1.5)
    pc = pc + geom_hline(yintercept = tstat$cs$nECDNA, linetype = "dashed", color = "red", linewidth = 1.5)
  }
  
  plt = ggarrange(pt, pc, nrow = 1, ncol = 2, common.legend = T)
  annotate_figure(plt, top = text_grob(title_all, face = "bold"))
  fout = file.path(sdir, paste0("ppc", suffix, ".pdf"))
  ggsave(fout, width = 9, height = 5)
}


plot_ppc_hist <- function(res_stat, sdir, suffix, title_all, tstat = list()){
  sim_stat = res_stat$sim_stat
  cell_stat = res_stat$cell_stat
  title = ""
  pt = ggplot(sim_stat, aes(x=avg_fusion, y=after_stat(density)))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr() + ggtitle(title) + xlab(lbl_fusion) + geom_vline(xintercept = mean(sim_stat$avg_fusion), linetype = "dashed", color = "blue", linewidth = 1.5)

  pc = ggplot(cell_stat, aes(x=nECDNA, y=after_stat(density)))  + geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = nbin) + geom_density(alpha=.5, fill="lightblue") + theme_pubr()  + ggtitle(title) + xlab(lbl_ecdna) + geom_vline(xintercept = mean(cell_stat$nECDNA), linetype = "dashed", color = "blue", linewidth = 1.5)  
  
  if(length(tstat) > 0){  
    ncell = unique(sim_stat$nCell)
    pt = pt + geom_vline(xintercept = tstat$ss$avg_fusion, linetype = "dashed", color = "red", linewidth = 1.5)
    pc = pc + geom_vline(xintercept = tstat$cs$nECDNA, linetype = "dashed", color = "red", linewidth = 1.5)
  }
  
  plt = ggarrange(pt, pc, nrow = 1, ncol = 2, common.legend = T)
  annotate_figure(plt, top = text_grob(title_all, face = "bold"))
  fout = file.path(sdir, paste0("ppc", suffix, "_hist.pdf"))
  ggsave(fout, width = 9, height = 5)
}


check_ppc_stat <- function(rdir, ddir, ncell, dataset, suffix, tstat = list(), plot_stat = F, plot_ppc = F){
  # check summary stats used for inference
  # header = c("pga", "mean_div", "sd_div", "pgaA", "mean_divA", "sd_divA", "pgaB", "mean_divB", "sd_divB", paste0("bp", c(1:ncell)))
  header_ppc1 = c("pga", "mean_div", "sd_div", "pgaA", "mean_divA", "sd_divA", "pgaB", "mean_divB", "sd_divB", paste0("bp", c(1:ncell)), "frac_wgd", "frac_del", "frac_dup", "frac_inv", "frac_tra")
  osuffix = paste0("_", dataset)
  
  fppc = file.path(rdir, "resample" ,paste0("resample", osuffix), paste0("check-smc", suffix))
  rpp = read_tsv(fppc, col_names = F, show_col_types = F)
  names(rpp) = header_ppc1
  # rpg = rpp %>% gather(all_of(header), key = "stat", value = "val")
  
  freal = file.path(ddir, paste0(dataset, "_sstat.tsv")) 
  sreal = read.table(freal) %>% t() %>% as.data.frame()
  names(sreal) = header_ppc1

  res_stat = get_ppc_sstat(file.path(rdir, "resample"), osuffix)

  if(plot_stat){
    check_summary_stats_boxplot(rpp, sreal, ncell, fppc)
    check_summary_stats_hist(rpp, sreal, ncell, fppc)
  }
  
  if(plot_ppc){
    title = dataset
    plot_ppc_hist(res_stat, rdir, suffix, title, tstat)
    plot_ppc_boxplot(res_stat, rdir, suffix, title, tstat)      
  }
  
  return(res_stat)
}


#################### functions for model selection #################### 

get_dic_all <- function(odir, pattern = "dic_nparam_epsilon0.2_ncell7"){
  files = list.files(odir, pattern = pattern)
  
  header = c("epsilon", "ncell", "model", "maxDSB", "dbreak" , "msample", "nsample", "DSB_target", "unrepair_target", "seed_target", "model_target", "dbreak_target", "run", "value")
  res_all = data.frame()
  for(i in 1:length(files)){
    # i = 1
    fname = file.path(odir, files[i])
    print(fname)
    res = read_tsv(fname, col_names = F, show_col_types = F)
    if(ncol(res) < length(header)) next
    # fields = strsplit(files[i], "_")[[1]]
    # type = str_replace(fields[length(fields)], "mp", "")
    # print(type)
    #
    # model = "neutral"
    # if(type == "1"){
    #   model = "selection"
    # }
    
    # print(ncol(res))
    names(res) = header
    # res$model = model
    res_all = rbind(res_all, res)
  }
  
  return(res_all)
}

get_dic_all_real <- function(odir, pattern){
  files = list.files(odir, pattern = pattern)
  
  res_all = data.frame()
  for(i in 1:length(files)){
    # i = 1
    fname = file.path(odir, files[i])
    print(fname)
    res = read_tsv(fname, col_names = F, show_col_types = F)
    
    # fields = strsplit(files[i], "_")[[1]]
    # type = str_replace(fields[length(fields)], "mp", "")
    # print(type)
    #
    # model = "neutral"
    # if(type == "1"){
    #   model = "selection"
    # }
    
    # print(ncol(res))
    if(ncol(res) == 10){
      names(res) = c("epsilon", "ncell", "model", "maxDSB", "dbreak" , "msample", "nsample", "dataset", "run", "value")     
    }else{
      names(res) = c("epsilon", "ncell", "model", "maxDSB", "dbreak" , "msample", "nsample", "run", "value")
    }
    
    # res$model = model
    res_all = rbind(res_all, res)
  }
  
  return(res_all)
}


# suffix is for name of output file
plot_msel_dic <- function(odir, files, suffix, ylimits = c(), ftype = ".png"){
  ggplot(dic_all, aes(x=as.factor(model), y=(value), fill=as.factor(model))) +
    # geom_violin(scale = "width") +
    geom_boxplot() +
    ylab("DIC")  + theme_pubr()  + xlab("") + theme(axis.text.x = element_blank()) + scale_fill_manual(breaks = c("neutral", "selection"), values = model_colours) + guides(fill = guide_legend("")) + ggtitle(title)
  
    p = ggplot(dic_all, aes(x=as.factor(model), y=(value), fill=as.factor(model))) +
      # geom_violin(scale = "width") +
      geom_boxplot() +
      ylab("DIC") +
      theme_pubr() + theme(axis.text.x = element_blank()) + guides(fill = guide_legend(nrow = 1)) + theme_pubr()  + xlab("") + theme(axis.text.x = element_blank()) + scale_fill_npg() + facet_grid(. ~ as.factor(seed_target)) + labs(fill = "")

    # if(length(ylimits) > 0){
    #   p = p + scale_y_continuous(limits = c(ylimits[1], ylimits[2]))
    # }

    ggplot(dic_all, aes(x=as.factor(dbreak), y=(value), fill=as.factor(dbreak))) +
      # geom_violin(scale = "width") +
      geom_boxplot() +
      ylab("DIC") +
      theme_pubr() + theme(axis.text.x = element_blank()) + guides(fill = guide_legend(nrow = 1)) + theme_pubr()  + xlab("") + theme(axis.text.x = element_blank()) + scale_fill_npg() + facet_grid(. ~ as.factor(seed_target)) + labs(fill = "")
    
  ftype = ".pdf"
  fout=file.path(odir, paste0("dic_boxplot", suffix, ftype))
  ggsave(fout, width = 8, height = 5)

}
