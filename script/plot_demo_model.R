#!/usr/bin/env Rscript

# This script is used to plot data simulated in Fig 1

# change path accordingly
source("parse_sim_util.R")

######################### set up #####################################
# simulate until 3 cells
dir1 = "data_svmodel/Fig1/r1"

# simulate until 2 cells
dir2 = "data_svmodel/Fig1/r2"

odir = "data_svmodel/Fig1"

######################### visualize simulated data in haplotype plot #####################################
cells = c(2)
for(c in cells){
  bdir = file.path(dir2, paste0("c", c))
  print(bdir)
  
  fsv = file.path(bdir, "SVs.CN_opt.phased")
  fcn = file.path(bdir, "copy_number.CN_opt.phased")   
  fplot = file.path(odir, paste0("ACN_c", c,".pdf"))
  
  break_thres=1000;
  
  SV=read.table(fsv,stringsAsFactors=F) %>% unique()

  t<-read.table(fcn, header=T)
  
  plot_chr_noCN(SV, t, fplot)
}


###############  check tree lineage ############### 
plot_lineage_with_cn(dir1)

# plot ancestor nodes as well, with total and A CN
ftree = file.path(dir1, "cell_lineage.tsv")
# aligning CN with tree
dtree = read_tsv(ftree, show_col_types = F)
ig = graph_from_edgelist(as.matrix(dtree[, c(1, 2)]), directed = TRUE)
mytree = as.phylo(ig)
p = ggtree(mytree) + ggtitle(title) + geom_text(aes(label=label), hjust = 1.2, vjust = -0.5) + xlim(-0.5, 2) 
layer_scales(p)$x$range$range. # help to show leftmost number, otherwise disappear
# layer_scales(p)$y$range$range # not work  

d = fortify(mytree)
slevel = d$label[order(d$y, decreasing = T)]

plot_cn_merged <- function(dir1, dir2, slevel, chr_var, ctype = "total", ncn = 2, title = "", theme1 = theme2){
  d_seg1 = get_cn_per_run(dir1, type = ctype)
  d_seg2 = get_cn_per_run(dir2, type = ctype)
  d_seg = rbind(d_seg1, d_seg2) %>% unique()  
  
  d_seg$sample = factor(d_seg$sample, levels = slevel)
  
  d_seg_sel = d_seg %>% filter(chrom %in% chr_var)
  
  d_norm = d_seg_sel %>% filter(sample == 4) %>% select(chrom, start, end) %>% mutate(sample = 1, cn = ncn)
  d_cn = rbind(d_seg_sel, d_norm)
  # plot not affected by different number of segments in different samples
  phmap = plot.cn.heatmap(d_cn, title, theme = theme1)
  
  return(phmap)
}

# different runs have different segments
d_seg1 = get_cn_per_run(dir1, type = "total")
d_seg2 = get_cn_per_run(dir2, type = "total")
d_seg = rbind(d_seg1, d_seg2) %>% unique()
chr_var = d_seg %>% filter(cn != 2) %>% select(chrom) %>% unique() %>% unlist()
chr_var
base_size = 8
phmap1 = plot_cn_merged(dir1, dir2, slevel, chr_var, "total", 2, "total CN", theme1) + theme(title = element_text(size = base_size))
phmap2 = plot_cn_merged(dir1, dir2, slevel, chr_var, "B", 1, "B CN", theme2) + theme(title = element_text(size = base_size))
phmap3 = plot_cn_merged(dir1, dir2, slevel, chr_var, "A", 1, "A CN", theme2) + theme(title = element_text(size = base_size))
#
(p | phmap1 | phmap3 | phmap2) + plot_layout(widths = c(1, 1, 1, 1))  
fout = file.path(odir, "demo_lcn.pdf")
ggsave(fout, width = 5, height = 3)


######################### visualize one run of simulated data  #####################################
plot_circos_per_run(dir1, "hg38", size = 4)

d_seg = get_cn_per_run(dir1, type = "total")
fout = file.path(dir1, "cnT_all.pdf")
main=""
pc = plot.cn.heatmap(d_seg, main)
ggsave(fout, pc, width = 12, height = 10)

