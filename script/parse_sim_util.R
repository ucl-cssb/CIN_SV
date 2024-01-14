# parse simulated summary 

#### make sure to load gUtils
library(gTrack)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(tools)
library(IRanges)
library(GenomicRanges)  
library(gGnome)   # failed to install on windows
# devtools::install_github('mskilab/fishHook') failed
# install.packages("/Users/bl0033/Downloads/mskilab-org-fishHook-c9c8587.tar", repos=NULL, type = "source")
library(patchwork) # for multiple plots with +/|
library(sqldf) # to find differences in data frame

# may need to specify the full path 
source("visualize_sv.R")
source("seismic_amplification_detection.R")

lbl_unrepair = "%unrepaired DSBs per cycle"
type_colors = c('#ece2f0','#a6bddb','#1c9099') # cell cycle
unrep_colors = c('#e0f3db','#a8ddb5','#43a2ca')
# mypal = pal_npg("nrc", alpha = 1)(10)
# model_colours = mypal[5:7]
# neutral, positive, negative
model_colours = c("#e0f3f8", "#fc8d59", "#4575b4")
model_breaks = c("neutral evolution", "positive selection", "negative selection")

model_colours2 = c("#e0f3f8", "#fee090")
model_breaks2 = c("neutral evolution", "selection")

cell_colors = c('#fee0d2','#fc9272','#de2d26')

lbl_fusion = "#chromosome fusions in a run"
lbl_ncell = "total number of cells"
lbl_scoef = "selection coefficient per cell"

############ parameters to detect chromothripsis ############ 
# High confidence: 
# at least 6 interleaved intrachromosomal SVs, 
# 7 contiguous segments oscillating between 2 CN states, 
# the fragment joins test, 
# and either the chromosomal enrichment
# or the exponential distribution of breakpoints test. 
pval = 0.05
nsv = 6
nosc_hconf = 7

# Low confidence: 
# at least 6 interleaved intrachromosomal SVs, 
# 4, 5 or 6 adjacent segments oscillating between 2 CN states, 
# the fragment joins test, 
# and either the chromosomal enrichment
# or the exponential distribution of breakpoints test.
nosc_lconf_l = 4
nosc_lconf_u = 6


############ common functions ############ 
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# Function to merge consecutive regions with the same copy number
# not work as regions with normal CN are assigned the same ID
merge_consecutive_regions <- function(df) {
  # assign regions to merge into the same group and then take the min and max distacne for each region/group
  CNData_svgen %>%
    group_by(chromosome, total_cn, cnA, cnB) %>% 
    arrange(chromosome, start) %>% mutate(dist1 = start - lag(end) - 1)  %>% mutate(dist2 = lead(start) - end - 1) %>% mutate(dist = if_else(dist1 == 0 | dist2 == 0, 0, dist1)) %>% mutate(tag = ifelse(start == 1 & is.na(dist), 0, dist)) %>%
    group_by(chromosome, total_cn, cnA, cnB, tag) %>% mutate(grp = cur_group_id()) %>%
    ungroup() %>%
    select(-diff_start, -grp)
}

############ plot functions ############ 

# plot one specific chromosome
plot_hap_cn <- function(d_seg, hap, chr, model, odir, sort_by_start = T, suffix = "", theme = theme, width = 12, height = 10){
  d_seg %>% filter(chrom == chr) -> dB7
  
  if(sort_by_start){
    dB7 %>% filter(start == 1) %>% arrange(end) -> dB_sel2
  }else{
    dB7 %>% filter(end == max(max(dB7$end))) %>% arrange(start) -> dB_sel2
  }
  
  dB7$sample = factor(dB7$sample, levels = dB_sel2$sample)
  if(hap == "T"){
    main=paste0("Total copy number")
  }else{
    main=paste0("Haplotype ", hap, " copy number")
  }
  
  p7 = plot.cn.heatmap(dB7, main, theme = theme)
  fout = file.path(odir, paste0("cn", hap, "_chr", chr, "_model", model, suffix, ".pdf"))
  ggsave(fout, p7, width = width, height = height)
}


# show node IDs to facilitate tracking
plot_lineage_with_cn_label <- function(dir, title = " ", theme = theme1){
  ftree = file.path(dir, "cell_lineage.tsv")
  
  # plot like in cnetml, aligning CN with tree
  dtree = read_tsv(ftree, show_col_types = F)
  ig = graph_from_edgelist(as.matrix(dtree[, c(1, 2)]), directed = TRUE)
  mytree = as.phylo(ig)
  p = ggtree(mytree) + ggtitle(title) + geom_text(aes(label=label), hjust=-.3)
  
  # read from txt file to keep node ID  
  
  d = fortify(mytree)
  slevel = d$label[order(d$y, decreasing = T)]
  # slevel = paste0("c", ordered_nodes)
  # length(slevel)
  
  d_seg = get_cn_per_run(dir, type = "total")
  # d_seg = d_seg %>% filter(!sample %in% c(2, 4, 5))
  d_seg$sample = factor(d_seg$sample, levels = slevel)
  
  # get the node order of the tree and reorder heatmap 
  phmap = plot.cn.heatmap(d_seg, "", theme = theme)
  
  # hard to align the height
  pc = ggarrange(p, phmap, nrow = 1, widths = c(0.5, 1))
  # pc
  
  # Error in grid.Call.graphics(C_setviewport, vp, TRUE) : variable names are limited to 10000 bytes
  # library(aplot)  
  # pc = phmap %>% insert_left(p, width = 0.3)
  # library(patchwork)
  # p + phmap + plot_layout(widths = c(1, 3))
  
  fout = file.path(dir, "cell_lineage_cn.pdf") 
  ggsave(fout, pc, width = 25, height = 16)
}


# show heatmap of selected chromosomes 
plot_lineage_with_cn_sel <- function(dir, chrs, odir, suffix, showT = T, title = " ", show_lbl = T, theme = theme1, w = 6, h = 3.5, plotB = F, tsize = 12, wsplit = c(1, 2)){
  ftree = file.path(dir, "cell_lineage.tsv")
  
  # plot like in cnetml, aligning CN with tree
  dtree = read_tsv(ftree, show_col_types = F)
  ig = graph_from_edgelist(as.matrix(dtree[, c(1, 2)]), directed = TRUE)
  mytree = as.phylo(ig)
  p = ggtree(mytree) + ggtitle(title) + theme(plot.title = element_text(size=tsize))
  if(show_lbl){
   p = p + geom_text(aes(label=label), hjust=-.3)
  }
  
  # read from txt file to keep node ID  
  
  d = fortify(mytree)
  slevel = d$label[order(d$y, decreasing = T)]
  # slevel = paste0("c", ordered_nodes)
  # length(slevel)
  
  if(showT){
    d_seg = get_cn_per_run(dir, type = "total")
    ctitle = "total CN"
  }else{
    d_seg = get_cn_per_run(dir, type = "B")
    ctitle = "B CN"
  }
  # d_seg = d_seg %>% filter(!sample %in% c(2, 4, 5))
  d_seg$sample = factor(d_seg$sample, levels = slevel)
  
  # get the node order of the tree and reorder heatmap 
  d_seg_2plot = d_seg %>% filter(chrom %in% chrs)
  phmap = plot.cn.heatmap(d_seg_2plot, ctitle, theme = theme, tsize = tsize) 
  
  if(plotB){
    # use RCK output for merged regions
    d_segB = get_cn_per_run_rck(dir, type = "B")    
    
    d_segB %>% filter(chrom %in% chrs) -> dB7
    dB7 %>% filter(start == 1) %>% arrange(end) -> dB_sel2
    dB7$sample = factor(dB7$sample, levels = dB_sel2$sample)  
    p7 = plot.cn.heatmap(dB7, "sorted B CN", theme = theme, tsize = tsize) 
    pc = p | phmap | p7
  }else{
    # hard to align the height
    # pc = ggarrange(p, phmap, nrow = 1, widths = wsplit)
    pc = p | phmap    
  }

  # Error in grid.Call.graphics(C_setviewport, vp, TRUE) : variable names are limited to 10000 bytes
  # library(aplot)  
  # pc = phmap %>% insert_left(p, width = 0.3)
  # library(patchwork)
  # p + phmap + plot_layout(widths = c(1, 3))
  
  fout = file.path(odir, paste0("cell_lineage_cn_sel", suffix,".pdf"))
  ggsave(fout, pc, width = w, height = h)
}


plot_cn_sorted <- function(odir, chr, model, dir, theme = theme, suffix = "_frag50", width = 5, height = 6){
  sort_by_start = T
  
  # the breakpoints are unique for each cell
  d_segT = get_cn_per_run_rck(dir, type = "total")
  d_segA = get_cn_per_run_rck(dir, type = "A")  
  d_segB = get_cn_per_run_rck(dir, type = "B")
  
  hap = "B"
  d_seg = d_segB
  plot_hap_cn(d_seg, hap, chr, model, odir, sort_by_start, suffix, theme, width, height)
  
  hap = "A"
  d_seg = d_segA
  plot_hap_cn(d_seg, hap, chr, model, odir, sort_by_start, suffix, theme, width, height)
  
  hap = "T"
  d_seg = d_segT
  plot_hap_cn(d_seg, hap, chr, model, odir, sort_by_start, suffix, theme, width, height)
}


plot_lineage_with_cn <- function(dir, odir, suffix = "", theme = theme1, show_tip = F, tsize = 6, w = 3.6, h = 3.5, sel_nodes = c(), wsplit = c(1, 2), htitle = ""){
  ftree = file.path(dir, "cell_lineage.tsv")
  # plot_lineage_tree(ftree)
  
  # plot like in cnetml, aligning CN with tree
  dtree = read_tsv(ftree, show_col_types = F)
  ig = graph_from_edgelist(as.matrix(dtree[, c(1, 2)]), directed = TRUE)
  mytree = as.phylo(ig)
  p = ggtree(mytree) + ggtitle(" ") 
  if(show_tip){ # for chromoplex plot in paper
    if(length(sel_nodes) > 0){
      ntip = mytree$Nnode + 1
      node_type0 = data.frame(node = 1:ntip, label = mytree$tip.label, chmplx_type = "N")
      node_type = node_type0 %>% mutate(chmplx_type = ifelse(label %in% sel_nodes, "Y", chmplx_type))
      # geom_tippoint(aes(shape = chmplx_type, size = 3))
      p = p %<+% node_type + geom_tiplab(hjust=0.1, size=tsize, aes(color = chmplx_type)) + theme(legend.position = "None") + scale_color_manual(values = c("grey", "black"))      
    }else{
      p = p + geom_tiplab(size=tsize)
    }
  }
 
  d = fortify(mytree)
  slevel = d$label[order(d$y, decreasing = T)]
  # slevel = paste0("c", ordered_nodes)
  # length(slevel)
  
  d_seg = get_cn_per_run(dir, type = "total")
  # d_seg = d_seg %>% filter(!sample %in% c(2, 4, 5))
  d_seg$sample = factor(d_seg$sample, levels = slevel)
  
  # get the node order of the tree and reorder heatmap 
  phmap = plot.cn.heatmap(d_seg, htitle, theme = theme)
  
  # hard to align the height
  #pc = ggarrange(p, phmap, nrow = 1, widths = c(0.5, 1))
  # heatmap not shown weirdly sometimes
  pc = (p | phmap) + plot_layout(widths = wsplit)
  
  # Error in grid.Call.graphics(C_setviewport, vp, TRUE) : variable names are limited to 10000 bytes
  # library(aplot)  
  # pc = phmap %>% insert_left(p, width = 0.3)
  # 
  # p + phmap + plot_layout(widths = c(1, 3))
  
  fout = file.path(odir, paste0("cell_lineage_cn", suffix, ".pdf")) 
  ggsave(fout, pc, width = w, height = h)
}


plot_circos_per_run <- function(dir, ref = "hg19", size = 5){
  fnames = list.files(dir, pattern = "SVData_c.*tsv")
  
  for(f in fnames){
    fsv = file.path(dir, f)  
    # print(fsv)
    # retrieve ID from file name
    bname = tools::file_path_sans_ext(f)
    fields = strsplit(bname, "_")[[1]]
    cell_ID = str_replace(fields[2], "c", "")
    
    fcnv = file.path(dir, list.files(dir, pattern = paste0("CNData_c", cell_ID, "_div.*.tsv")))
    print(fsv)
    print(fcnv)
    plot_CN_SV(fsv, fcnv, ref, size)
  }
}

# visualize in circos plot
plot_in_circos <- function(bdir, n_dsbs, n_unrps, n_frag, n_break, nrun, ref = "hg19", pair_type = 0, prob_corr = 0){
  for(n_dsb in n_dsbs){
    # n_dsb = 50
    print(n_dsb)
    for(n_unrp in n_unrps){
      print(n_unrp) 
      # for(i in 1:nrun){
      # i = 8 
      # n_dsb = 100
      # n_unrp = 5
      # print(i)
      suffix = paste0("/nDSB", n_dsb, "_Un", n_unrp)
      if(n_frag > -1){
        suffix = paste0(suffix, "_frag", n_frag)         
      }
      if(n_break > - 1){
        suffix = paste0(suffix, "_break", n_break)
      }
      if(pair_type == 1){
        suffix = paste0(suffix, "_ptype", pair_type, "_pcorr", prob_corr)
      }
      dir = file.path(bdir, suffix) 
      dir = file.path(dir, paste0("r", nrun))
      
      print(dir)
      # # use cell_ID to match files
      # for(cell_ID in cell_IDs){
      #   print(cell_ID)
      #   fsv = file.path(dir, list.files(dir, pattern = paste0("SVData_c", cell_ID, "_div.*.tsv")))
      #   #plot_SV(fsv)
      #   fcnv = file.path(dir, list.files(dir, pattern = paste0("CNData_c", cell_ID, "_div.*.tsv")))
      #   print(fsv)
      #   print(fcnv)
      #   plot_CN_SV(fsv, fcnv)
      # }
      plot_circos_per_run(dir, ref)
    }
    # }
  }
}


# plot CN in one chr of one cell
plot_cell_cn_chr <- function(cn_dir, fcn, odir, suffix, chrsel, xgap = 1000){
  cnsel = get_CN(file.path(cn_dir, fcn))
  seg.cn = cnsel %>% filter(chr == chrsel) %>% dplyr::select(chromosome = chr, start.pos = start, end.pos = end, A = value1, B = value2) %>% mutate(CNt = A + B) %>% as.data.frame()
  seg.cn$chromosome = str_replace(seg.cn$chromosome, "chr", "") 
  # seg.cn$chromosome = factor(seg.cn$chromosome, levels = paste0(seq(1:22)))
  fout = file.path(odir, str_replace(fcn, ".tsv", paste0("_acn", suffix, ".pdf")))
  pdf(fout, width = 8, height = 4)
  plot_genome_cn(seg.cn, info.type = "AB", xgap = xgap)
  dev.off()  
}


# plot chrsel in local regions of one cell
plot_cell_cn_local <- function(cn_dir, fcn, odir, suffix, chrsel, xgap = 100, ygap = 20){
  cnsel = get_CN(file.path(cn_dir, fcn))
  seg.cn = cnsel %>% filter(chr == chrsel) %>% dplyr::select(chromosome = chr, start.pos = start, end.pos = end, A = value1, B = value2) %>% mutate(CNt = A + B) %>% as.data.frame()
  seg.cn$chromosome = str_replace(seg.cn$chromosome, "chr", "") 
  # seg.cn$chromosome = factor(seg.cn$chromosome, levels = paste0(seq(1:22)))
  fout = file.path(odir, str_replace(fcn, ".tsv", paste0("_acn", suffix, ".pdf")))
  
  # merge regions to exclude normal ones
  seg_normal = seg.cn %>% filter(A == 1 & B == 1)
  loc_range = GRanges(seg_normal$chromosome, IRanges(seg_normal$start.pos, seg_normal$end.pos))
  loc_merged = reduce(loc_range)
  normal_ranges = ranges(loc_merged)
  print(normal_ranges)
  if(length(normal_ranges) > 0){
    nstart = normal_ranges[1]@start
  }else{
    nstart = max(seg.cn$end.pos) + 1
  }
  
  ref_start = 1
  ref_end = nstart - 1
  seg_plot = seg.cn %>% filter(start.pos < nstart)
  pdf(fout, width = 8, height = 4)
  plot_allele_cn(seg_plot, info.type = "AB", ref_start, ref_end, xgap, ygap)
  dev.off()  
}


plot_cell_with_maxcn <- function(cn_dir, odir, suffix, chrsel, top = 1, is_local = T, xgap = 100, ygap = 20){
  cn_max = read.table(file.path(cn_dir, "cn_max"))
  if(top == 1){
    mcn = max(cn_max$V2)
  }else{
    mcn = cn_max %>%  arrange(desc(V2)) %>% slice(top) %>% dplyr::select(V2) %>% unlist()
  }
  print(mcn)
  
  fcns = cn_max %>% filter(V2 == mcn) %>% dplyr::select(V1) %>% unlist()
  
  for(fcn in fcns){
    if(is_local){
      plot_cell_cn_local(cn_dir, fcn, odir, suffix, chrsel, xgap, ygap)
    }else{
      plot_cell_cn_chr(cn_dir, fcn, odir, suffix, chrsel)
    }
  }
}


################## functions related to chromothripsis ##################
get_chromothripsis <- function(dir, midfix) {
  fsv <- file.path(dir, paste0("SVData_", midfix, ".tsv"))
  fcn <- file.path(dir, paste0("CNData_", midfix, ".tsv"))
  fplot <- file.path(dir, paste0("plot_chr1_", midfix, ".pdf"))
  
  SVData_svgen <- read.delim(fsv)
  CNData_svgen <- read.delim(fcn)
  
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
  chromothripsis <- shatterseek(SV.sample = SV_data, seg.sample = CN_data, genome = "hg38")
  end_time <- Sys.time()
  print(paste0("Running time (s): ", round(end_time - start_time, digits = 2)))
  
  # plot results
  plots_chr <- plot_chromothripsis(ShatterSeek_output = chromothripsis, chr = "1", genome = "hg38")
  p <- arrangeGrob(plots_chr[[1]],
                   plots_chr[[2]],
                   plots_chr[[3]],
                   plots_chr[[4]],
                   nrow = 4, ncol = 1, heights = c(0.2, .4, .4, .4)
  )
  pdf(fplot)
  plot_grid(p)
  dev.off()
  
  return(chromothripsis)
}


# specify the directory directly 
get_chromth_count_by_dir <- function(dir, chrs_sel = c()){
  df_chromth = data.frame()
  #print(dir)
  fnames = list.files(dir, pattern = "SVData_c.*tsv")
  #print(fnames)
  
  for(f in fnames){
    fsv = file.path(dir, f)  
    #print(fsv)
    
    # retrieve ID from file name
    bname = tools::file_path_sans_ext(f)
    fields = strsplit(bname, "_")[[1]]
    cell_ID = str_replace(fields[2], "c", "")
    div_ID = str_replace(fields[3], "div", "") 
    
    svs = read_tsv(fsv, show_col_types = F) 
    
    if(length(chrs_sel) == 0){
      c1 = svs %>% group_by(chrom1) %>% tally() %>% filter(n > nsv) %>% dplyr::select(chrom1) %>% unlist()
      c2 = svs %>% group_by(chrom2) %>% tally() %>% filter(n > nsv) %>% dplyr::select(chrom2) %>% unlist()
      chrs_sel = union(c1, c2)          
    }
    
    for(chr in chrs_sel){
      #print(chr)
      midfix <- paste0("c", cell_ID, "_div", div_ID, "_chr", chr)
      fname <- file.path(dir, paste0("shatterseek_", midfix, ".tsv"))
      #print(fname)
      
      # if(!file.exists(fname)){
      #   cmd = paste0("Rscript /Users/bl0033/Gdrive/git/CIN_SV/script/detect_chromothripsis.R  -d ", dir, " -n ", div_ID, " -c ", cell_ID, " -r ", chr)
      #   stat = sys::exec_internal(cmd, error = F)
      #   if(stat$status == 1) next             
      # }
      
      if(!file.exists(fname)){
        # cat("file not exist! ", fname, "\n")
        df = data.frame(cell_ID = cell_ID, div_ID = div_ID, chrom = chr,  nhconf = 0, nlconf = 0)
        df_chromth = rbind(df_chromth, df)
        next
      }
      
      chromothripsis = read_tsv(fname, show_col_types = F)
      #print(chromothripsis)
      
      if(nrow(chromothripsis) == 0){
        # cat("empty file! ", fname)
        df = data.frame(cell_ID = cell_ID, div_ID = div_ID, chrom = chr, nhconf = 0, nlconf = 0)
        df_chromth = rbind(df_chromth, df)
        next
      } 
      
      # chromothripsis %>% t() %>% View()
      # chromothripsis %>% dplyr::select(starts_with("pval")) %>% View() 
      
      chromothripsis %>% filter(clusterSize_including_TRA - number_TRA >= nsv) %>%  filter(pval_fragment_joins < pval)  %>% filter(pval_exp_chr < pval | pval_exp_cluster < pval) -> chrom_events
      
      chrom_events %>% filter(max_number_oscillating_CN_segments_2_states >= nosc_hconf) -> hconf 
      chrom_events %>% filter(max_number_oscillating_CN_segments_2_states >= nosc_lconf_l & max_number_oscillating_CN_segments_2_states <= nosc_lconf_u) -> lconf
      
      # nDSB = n_dsb, nUnrepair = f_unrp, nfrag = frag, run = i,
      df = data.frame(cell_ID = cell_ID, div_ID = div_ID, chrom = chr, nhconf = nrow(hconf), nlconf = nrow(lconf))
      # print(df)
      df_chromth = rbind(df_chromth, df)
    }
  }

  return(df_chromth)
}


# parse shatterseek output to count chromothripsis by confidence 
# for multiple runs with fixed parameters n_dsbs, f_unrps, frags_local
get_chromth_count <- function(bdir, n_dsbs, f_unrps, frags_local, nrun, chrs_sel = c(), model = 0, dsuffix = ""){
  df_chromth = data.frame()
  
  for(n_dsb in n_dsbs){
    # n_dsb = 50
    # print(n_dsb)
    for(f_unrp in f_unrps){
      # print(f_unrp)  
      for(frag in frags_local){
        # print(frag)       
        for(i in 1:nrun){
          # print(i)
          dprefix = paste0("nDSB", n_dsb, "_Un", f_unrp, "_frag", frag)
          if(model == 1){
            dprefix = paste0(dprefix, "_model", model)
          }
          dprefix = paste0(dprefix, dsuffix)
          dir = file.path(bdir, dprefix, paste0("r", i))
          # print(dir)
          
          chromth_count = get_chromth_count_by_dir(dir, chrs_sel)
          if(nrow(chromth_count) == 0){
            stop("no output!")
          }
          df = data.frame(nDSB = n_dsb, nUnrepair = f_unrp, nfrag = frag, run = i)
          # print(chromth_count)
          # print(df)
          res_df = cbind(chromth_count, df)
          
          df_chromth = rbind(df_chromth, res_df)
        }
      }
    }
  }
  
  return(df_chromth)
}


# df_chromth = get_chromth_count(bdir, n_dsbs, f_unrps, frags_local, nsv, nrun, chrs, model, dsuffix)
write_chromth_summary <- function(bdir, df_chromth, dsuffix = ""){
  fout_count = file.path(bdir, paste0("shatterseek_count", dsuffix, ".tsv")) 
  write_tsv(df_chromth, fout_count)
 
  if("nDSB" %in% colnames(df_chromth)){
    df_chromth %>% filter(nhconf > 0 | nlconf > 0) %>% group_by(nDSB, nUnrepair, nfrag, run, cell_ID, div_ID) %>% tally() %>% mutate(nchrmth = ifelse(n >= 1, 1, 0)) -> summ  
    df_chromth %>% filter(nhconf > 0) %>%  group_by(nDSB, nUnrepair, nfrag, run, cell_ID, div_ID) %>% tally() %>% mutate(nchrmth = ifelse(n >= 1, 1, 0)) -> summh    
  }else{
    df_chromth %>% filter(nhconf > 0 | nlconf > 0) %>%  group_by(cell_ID, div_ID) %>% tally() %>% mutate(nchrmth = ifelse(n >= 1, 1, 0)) -> summ
    df_chromth %>% filter(nhconf > 0) %>%  group_by(cell_ID, div_ID) %>% tally() %>% mutate(nchrmth = ifelse(n >= 1, 1, 0)) -> summh 
  }

  fout = file.path(bdir, paste0("shatterseek_summ", dsuffix, ".tsv"))
  write_tsv(summ, fout)

  fout = file.path(bdir, paste0("shatterseek_summ_hconf", dsuffix, ".tsv"))
  write_tsv(summh, fout)
}


get_chromth_count_wrap <- function(bdir1, n_dsbs, f_unrps, frags_local, nrun, chrs_sel){
  f1 = file.path(bdir1, "df_chromth_count.tsv")
  if(!file.exists(f1)){
    df_chromth1 = get_chromth_count(bdir1, n_dsbs, f_unrps, frags_local, nrun, chrs_sel)
    write_chromth_summary(bdir1, df_chromth1)
    write_tsv(df_chromth1, f1)
  }else{
    df_chromth1 = read_tsv(f1)
  }
  return(df_chromth1)  
}
  

get_chromth_all <- function(bdir1, bdir2, bdir3, t1, t2, t3, n_dsbs, f_unrps, frags_local, nrun, chrs_sel){
  df_chromth_all = data.frame()
  
  # div_ID = 1
  # cell_IDs = c(2,3)
  df_chromth1 = get_chromth_count_wrap(bdir1, n_dsbs, f_unrps, frags_local, nrun, chrs_sel)
  df_chromth1$type = t1
  df_chromth_all = rbind(df_chromth_all, df_chromth1)
  
  # bdir = "/Users/bl0033/data/SV/test4paper/sv_type"
  # div_ID = 2
  # cell_IDs = c(4,5)
  df_chromth2 = get_chromth_count_wrap(bdir2, n_dsbs, f_unrps, frags_local, nrun, chrs_sel)
  df_chromth2$type = t2
  df_chromth_all = rbind(df_chromth_all, df_chromth2)
  
  # bdir = "/Users/bl0033/data/SV/sim/test_chromth_ncell3_break1"
  # div_ID = 2
  # cell_IDs = c(4,5)
  df_chromth3 = get_chromth_count_wrap(bdir3, n_dsbs, f_unrps, frags_local, nrun, chrs_sel)
  df_chromth3$type = t3
  df_chromth_all = rbind(df_chromth_all, df_chromth3) 
  
  return(df_chromth_all)
}


################## functions related to seismic amplifications ##################
fban = "/Users/bl0033/Gdrive/git/CIN_SV/data/cytoBand_hg38.txt"
MV_BANDS = read_tsv(fban, col_names = F, show_col_types = F) 
names(MV_BANDS) = c("chrom", "chromStart", "chromEnd", "name", "gieStain")
chrBands = makeGRangesFromDataFrame(MV_BANDS, keep.extra.columns=TRUE)


get_seismic_result <- function(odir, cfiles){
  mcns = c() # track max CN
  res_list = list()
  for(i in 1:length(cfiles)){
    f = cfiles[i]
    # print(f)
    suffix = str_replace(tools::file_path_sans_ext(f), "CNData_", "")
    # suffix = "c17_div3"
    print(suffix)
    
    fcn = file.path(odir, paste0("CNData_", suffix, ".tsv")) 
    MY_CNVs = read_tsv(fcn, show_col_types = F) 
    mcn = max(MY_CNVs$total_cn)
    mcns = c(mcns, mcn)
    if(mcn > 10){
      message("Max CN: ", mcn)
    }
    MY_CNVs$chromosome = paste0("chr", MY_CNVs$chromosome)
    MY_CNVs = MY_CNVs %>% dplyr::rename(cn = total_cn)
    cnv = makeGRangesFromDataFrame(MY_CNVs, keep.extra.columns=TRUE)
    
    fsv = file.path(odir, paste0("SVData_", suffix, ".tsv"))
    MY_SVS = read_tsv(fsv, show_col_types = F)
    sv = MY_SVS %>% dplyr::select(chr1 = chrom1, bp1 = start1, chr2 = chrom2, bp2 = end2)
    sv$chr1 = paste0("chr", sv$chr1)
    sv$chr2 = paste0("chr", sv$chr2)
    
    # minInternalSVs=14
    # ploidy=2 
    # cnvTol=5000 
    # noXY=TRUE
    sa = detect_seismic_amplification(cnv=cnv, sv=sv, chrBands=chrBands)
    # sa$amplicons[sa$amplicons$id==1,]@seqnames
    res_list[[i]] = sa
  }
  
  return(list(res_list = res_list, mcns = mcns))
}



get_seis_list <- function(fout, dir, cfiles){
  if(file.exists(fout)){
    ###### retrieve the results
    res_list = readRDS(fout)    
  }else{
    ###### compute the results
    res = get_seismic_result(dir, cfiles)
    res_list = res$res_list
    # mcns = res$mcns   
    write_rds(res_list, fout)
  }
  return(res_list)
}



get_cell_with_seismic <- function(cfiles, res_list, dir, chrsel, plot, sel_id = NA){
  nsm = 0
  cell_ids = data.frame()
  for(i in 1:length(cfiles)){
    f = cfiles[i]
    suffix = str_replace(tools::file_path_sans_ext(f), "CNData_", "")
    # if(suffix != "c169_div9") next
    sa = res_list[[i]]
    if(nrow(sa$svs) > 0){
      print(suffix)
      # print(mcns[i])
      fields = str_split(suffix, "_")[[1]]
      cid = str_replace(fields[1], "c", "")
      did = str_replace(fields[2], "div", "")
      print(cid)
      print(sa)
      
      if(!is.na(sel_id) & cid == sel_id){
        return(sa)
        break
      }

      df = data.frame(cell = cid, div = did)
      cell_ids = rbind(df, cell_ids)
      
      ####### plot of type circos and sequenza
      if(plot){
        fname_cn = paste0("CNData_", suffix, ".tsv")
        plot_cell_cn_chr(dir, fname_cn, dir, "", chrsel)
        
        fsv = file.path(dir, paste0("SVData_", suffix, ".tsv"))
        #plot_SV(fsv)
        fcnv = file.path(dir, fname_cn)
        print(fsv)
        print(fcnv)
        plot_CN_SV(fsv, fcnv)     
      }
      nsm = nsm + 1
    }
  }  
  print(nsm)
  
  return(cell_ids)
}


# get the number of cells with Seismic Amplification in each run
# allow to only consider cells at certain cycle
get_df_sm_batch <- function(bdir, chrsel, plot, div_ID = -1){
  df_sm = data.frame()
  for(ncell in 2:3){
    for(div_break in 0:1){
      dir = file.path(bdir, paste0("test_ncell", ncell, '_break', div_break))  
      print(dir)
      if(!file.exists(dir)) next
      for(ndsb in c(10, 30)){
        for(funrep in c(0, 0.1, 0.3)){
          for(local_frag in c(0, 10, 30)){
            # ncell = 2
            # div_break = 0
            # ndsb = 10
            # funrep = 0.1
            # local_frag = 10
            ddir = file.path(dir, paste0("nDSB", ndsb, "_Un", funrep, "_frag", local_frag))
            for(i in 1:50){
              #i = 1
              subdir = file.path(ddir, paste0("r", i))
              
              if(div_ID >= 0){
                if(ncell == 2){
                  cfiles = list.files(subdir, pattern = "^CNData_c.*_div1.tsv")
                }else{
                  if(ncell != 3){
                    stop("Wrong number of cells!")
                  }
                  cfiles = list.files(subdir, pattern = "^CNData_c.*_div2.tsv")
                }
              }else{
                cfiles = list.files(subdir, pattern = "^CNData_c.*tsv")
              }
              #length(cfiles)
              
              fields = strsplit(subdir, "/")[[1]]
              nf = length(fields)
              s = nf - 2
              suffix = paste(fields[s:nf], collapse = "_")
              fout = file.path(subdir, paste0("seismic_", suffix, ".rds"))
              res_list = get_seis_list(fout, subdir, cfiles)
              cell_seis = get_cell_with_seismic(cfiles, res_list, subdir, chrsel, plot)
              df = data.frame(ncell = ncell, div_break = div_break, ndsb = ndsb, funrep = funrep, local_frag = local_frag, run = i, nsm = nrow(cell_seis), nfile = length(cfiles))
              df_sm = rbind(df_sm, df)
            }
          }
        }
      }
    }
  }
  return(df_sm)
}

################## functions related to ecDNA ##################
plot_ecdna_distr <- function(df_ecnda_all, nrun, odir, suffix, nrow = 2, ncol = 5){
  plts = list()
  # 5 / 10 has selection shifting toward higher values, 3 neutral, 2 mixed
  for(i in 1:nrun){
    df_ecnda_sel = df_ecnda_all %>% filter(run == i)
    # p = ggplot(df_ecnda_sel, aes(x=(nECDNA), fill = as.factor(model)))  + geom_histogram(aes(y=..density..), alpha=0.5, position = 'identity') + geom_density(alpha=.5) + theme_pubr() + scale_fill_manual(breaks = c("neutral", "selection"), values = model_colours) + guides(fill = guide_legend(""))  + xlab("ecDNA copy number") + ylab("Frequency")
    p = ggplot(df_ecnda_sel, aes(x=(nECDNA), fill = as.factor(model)))  + geom_histogram(aes(y=..density..), alpha=0.5, position = 'identity', bins = 30) + geom_density(alpha=.5) + theme_pubr(base_size = 7) + scale_fill_manual(breaks = c("neutral", "positive selection", "negative selection"), values = model_colours) + guides(fill = guide_legend(""))  + xlab("#ecDNAs") + ylab("Frequency")  + theme(axis.text.x = element_text(angle = 45, hjust=1))
    # + scale_y_continuous(limits = c(0, 0.055))
    plts[[i]] = p
  }
  ggarrange(plotlist = plts, nrow = nrow, ncol = ncol, common.legend = T)
  
  fout = file.path(odir, paste0("ecDNA_cn_density_", suffix, ".pdf"))
  ggsave(fout, width = 7.08, height = 4)
}


get_ec_count <- function(dir){
  files = list.files(dir, pattern = "genome_.*tsv")
  nf = length(files)
  df_eccount = data.frame()
  for(i in 1:nf){
    fname = files[i]
    fields = strsplit(fname, "_")[[1]]
    cid = str_replace(fields[2], "c", "")
    fpath = file.path(dir, fname)
    gnome = read_tsv(fpath, show_col_types = F)
    unique(gnome$shape)
    ecdnas = gnome %>% filter(shape == "circular")
    ecdna_count = ecdnas %>% group_by(nodes, size) %>% tally() %>% mutate(cell = cid)
    df_eccount = rbind(df_eccount, ecdna_count)
    # hist(ecdna_count$n)
  }
  return(df_eccount)
}

plot_selcoef_distr <- function(df_sstat_selection, nrun, odir, suffix){
  plts = list()
  for(i in 1:nrun){
    df_sel = df_sstat_selection %>% as.data.frame() %>% filter(run == i)
    p = ggplot(df_sel, aes(x=(selection_coef), fill = as.factor(model)))  + geom_histogram(aes(y=..density..), alpha=0.5, position = 'identity', bins = 30) + geom_density(alpha=.5) + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend(""))  + xlab("selection coefficient") + ylab("Frequency") 
    # + scale_y_continuous(limits = c(0, 8.5)) + scale_x_continuous(limits = c(-2, 5))
    plts[[i]] = p
  }
  ggarrange(plotlist = plts, nrow = 2, ncol = 5, common.legend = T)
  
  fout = file.path(odir, paste0("ecDNA_selcoef_", suffix, ".pdf"))
  ggsave(fout, width = 7.08, height = 4)
}


# dir: selection results, dir2: neutral results 
get_ecdna_both <- function(dir, dir2, sel){
  df_eccount = get_ec_count(dir)
  df_eccount_neu = get_ec_count(dir2)
  
  # df_eccount %>% filter(n==0) %>% select(cell) %>% unique()
  # summary(df_eccount$n)
  
  cell_pos = sel %>% filter(selection_coef > 0)
  cell_neg = sel %>% filter(selection_coef < 0)
  # nrow(cell_pos)
  # nrow(cell_neg)
  # summary(sel)
  if(nrow(cell_pos) + nrow(cell_neg) != nrow(sel)){
    stop("Wrong number of cells!")
  }
  ec_sel = merge(df_eccount, sel, by = c("cell")) 
  ec_pos = ec_sel %>% filter(selection_coef > 0)
  ec_neg = ec_sel %>% filter(selection_coef < 0)
  
  ec_neu = df_eccount_neu %>% mutate(type = "neutral evolution", prob_survival = 1, selection_coef = 0) %>% select(colnames(ec_sel))
  
  ec_all = rbind(ec_sel, ec_neu)
  
  return(ec_all)
}


plot_ecdna_ncopy <- function(ec_all, necdna_all, ec_sel, odir, suffix, mcn, maxn, mins, maxs, w = 7.5, h = 3, title = ""){
  # 3 groups, to be distinguished by color 
  # pn = ggplot(ec_all, aes(x = n, group = cell, color = type, alpha = 0.2)) + geom_density() + theme(legend.position = "none") + theme_pubr() + xlab("copy number of different ecDNAs per cell") + facet_grid(. ~ type) + scale_color_manual(breaks = model_breaks, values = model_colours) + theme_pubr(base_size = 7)  + theme(legend.position = "none")  + scale_x_continuous(limits = c(0, mcn)) 
  # + scale_y_continuous(expand = c(0, 0))
  # library(bayesplot)
  # ec_cn = ec_all %>% select(cell, n, type)
  # # ec_cn %>% group_by(cell, type) %>% tally()
  # # find cell with max CN for each type
  # topn = ec_cn %>% group_by(type) %>% top_n(1, n)
  # top_cell = topn %>% select(cell, type)
  # top_cn = merge(top_cell, ec_cn)
  # y = top_cn$n
  # group = top_cn$type
  # yrep = y
  #ppc_dens_overlay_grouped(y, yrep[1, ], group = group)
  pn = ggplot(ec_all, aes(x = log10(n), group = cell, fill = type, alpha = 0.2)) + geom_histogram(bins = 30, position = 'identity')  + theme(legend.position = "none") + theme_pubr() + xlab("copy number of different ecDNAs per cell") + facet_grid(. ~ type) + scale_fill_manual(breaks = model_breaks, values = model_colours) + theme_pubr(base_size = 7)  + theme(legend.position = "none")  + scale_x_continuous(limits = c(0, mcn)) 
 #  ggplot(ec_grp, aes(x = (n), y = nn, group = cell, fill = type, alpha = 0.2)) + geom_point()  + theme(legend.position = "none") + theme_pubr() + xlab("copy number of different ecDNAs per cell") + facet_grid(. ~ type) + scale_fill_manual(breaks = model_breaks, values = model_colours) + theme_pubr(base_size = 7)  + theme(legend.position = "none") + scale_x_continuous(trans = 'log10') + scale_y_continuous(trans = 'log10')
 # # , limits = c(0, log10(mcn))
 #  mcn = max(ec1$nn)
 #  ec1 = ec_all %>% filter(cell == 9999) %>% group_by(type, cell, n) %>% tally()
 #  ec_grp = ec_all %>% group_by(type, cell, n) %>% tally()
 #  ec1 %>% group_by(type, cell) %>% tally() %>% View()
  # for bfb local frag results
  # + geom_hline(yintercept = 0.2, linetype = "dashed", color = "darkgrey") + geom_hline(yintercept = 0.6, linetype = "dashed", color = "darkgrey") + geom_vline(xintercept = 15, linetype = "dashed", color = "darkgrey")  + geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkgrey") 
  
 # ps = ggplot(ec_sel, aes(x=(selection_coef), fill = as.factor(type)))  + geom_histogram(aes(y=..density..), alpha=0.5, position = 'identity', bins = 30) + geom_density(alpha=.5) + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend(""))  + xlab(lbl_scoef) + xlab("frequency") + theme(legend.position = "none") + scale_x_continuous(limits = c(mins, maxs))  + scale_y_continuous(expand = c(0, 0))
  ps = ggplot(ec_sel, aes(x = type, y = selection_coef, fill = type)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend("")) + ylab(lbl_scoef) + theme_pubr(base_size = 7) + theme(legend.position = "none") + xlab("") + theme(axis.text.x = element_blank()) + scale_y_continuous(limits = c(mins, maxs)) 
  
  # check different number of ecDNAs per cell, colored by type 
  # ec_all %>% select(cell, nodes, type) %>% unique() %>% dim()
  # overlap, not informative
  # pd = ggplot(necdna_all, aes(x = n, fill = type))  + geom_histogram(aes(y=..density..), alpha=0.5, position = 'identity', bins = 30) + geom_density(alpha=.5) + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend("")) + xlab("#different ecDNAs per cell") + theme(legend.position = "none") + ggtitle(title) + scale_x_continuous(limits = c(0, maxn)) + scale_y_continuous(expand = c(0, 0))
  
  pd = ggplot(necdna_all, aes(x = type, y = n, fill = type)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend("")) + ylab("#different ecDNAs per cell") + theme_pubr(base_size = 7) + theme(legend.position = "none") + xlab("") + theme(axis.text.x = element_blank()) + scale_y_continuous(limits = c(0, maxn)) + ggtitle(title)
  
  # ggarrange(pn, ps, nrow = 1) 
  # , guides = "collect"
  pall = pd + pn + ps  + plot_layout(widths = c(1, 3, 1))
  
  fout = file.path(odir, paste0("ncopy_", suffix, ".pdf"))
  ggsave(fout, pall, width = w, height = h)
  
  return(pall)
}


parse_ecdna_batch <- function(bdir, ncell, mcn, maxn, mins, maxs, ptitle = ""){
  dir = file.path(bdir, paste0("selection"))  # selection
  dir2 = file.path(bdir, paste0("neutral"))  # neutral evolution
  fsel = file.path(dir, "selection.tsv")
  sel = read.table(fsel, header = T) %>% mutate(type = if_else(selection_coef > 0, "positive selection", "negative selection"))
  # take a while for 1000 cells
  feall = file.path(bdir, paste0("ecdna_all_ncell", ncell, ".tsv"))
  if(file.exists(feall)){
    ec_all = read_tsv(feall)
  }else{
    ec_all = get_ecdna_both(dir, dir2, sel)  
    write_tsv(ec_all, feall)
  }
  ec_all = ec_all %>% mutate(tcell = ncell)
  ec_sel = ec_all %>% filter(selection_coef != 0) %>% mutate(tcell = ncell)
  # necdna_all = ec_all %>% group_by(cell, type) %>% tally() %>% mutate(tcell = ncell)
  
  print(max(ec_all$n))
  # print(max(necdna_all$n))
  print(min(sel$selection_coef))
  print(max(sel$selection_coef))
  
  # count #cells for info
  nec = ec_all %>% select(cell, type) %>% unique() %>% group_by(type) %>% tally()
  ntype = sel %>% group_by(type) %>% tally()
  ncount = merge(nec, ntype, by = c("type"), all = T) %>% mutate(n.y = if_else(is.na(n.y), ncell, n.y)) %>% mutate(prop = n.x / n.y, ncell = ncell)
  
  return(list(ncount = ncount, ec_sel = ec_sel, ec_all = ec_all))
}


# all ecdnas are on chr2
find_ecdna_gene_overlap <- function(ec_all, target.gr){
  ec_nodes = ec_all %>% select(cell, nodes, type, n) %>% separate_rows(nodes, sep = ",") %>% filter(nodes != "")
  ec_pos = ec_nodes %>% mutate(nodes = str_replace_all(nodes, "2:B:", "")) %>% separate(nodes, sep = "-", into = c("start", "end")) %>% mutate(start = as.numeric(start), end = as.numeric(end))  %>% filter(abs(start - end) > 1) %>% mutate(new_start = if_else(start > end, end, start), new_end = if_else(start > end, start, end)) %>% select(cell, s = (new_start), e = (new_end), type, n) %>% unique()
  query.gr = GRanges(seqnames = "2", ranges = IRanges(ec_pos$s, ec_pos$e))
  nolp = sum(countOverlaps(query.gr, target.gr))
  
  gene_overlap = data.frame()
  if(nolp > 0){
    overlaps = findOverlaps(query.gr, target.gr)
    # print(overlaps)
    olp_query = ec_pos[overlaps@from,]
    olp_target = csep_sel[overlaps@to,]
    df = bind_cols(olp_query, olp_target)
    gene_overlap = rbind(gene_overlap, df)
  }
  return(gene_overlap)
}



############ functions to read input ############ 
# read CNs of all cells in a run
get_cn_per_run <- function(dir, type = "total"){
  fcns = list.files(dir, pattern = "CNData_c.*_div.*.tsv")
  
  d_seg = data.frame()
  for(f in fcns){
    # f = fcns[1]
    fname = file.path(dir, f)
    cns = read_tsv(fname, show_col_types = F)  
    
    bname = file_path_sans_ext(f)
    fields = str_split(bname, pattern = "_")[[1]]
    # sid = str_replace(bname, "CNData_", "")
    # sid = paste0(fields[[3]], fields[[2]], sep = "_")
    sid = str_replace(fields[[2]], "c", "")
    cns$sample = sid
    
    
    if(type == "total"){
      seg = cns %>% dplyr::select(sample, chrom = "chromosome", start, end, cn = total_cn)
    }else if(type == "A"){
      seg = cns %>% dplyr::select(sample, chrom = "chromosome", start, end, cn = cnA)
    }else{
      # type == "B"
      seg = cns %>% dplyr::select(sample, chrom = "chromosome", start, end, cn = cnB)
    }
    
    # NORM_PLOIDY = 2
    # mode_cn = Mode(seg$cn)
    # ploidy = ifelse(mode_cn < NORM_PLOIDY, NORM_PLOIDY, mode_cn)
    # seg$ploidy = ploidy
    
    d_seg = rbind(d_seg, seg)
  }
  
  return(d_seg)
}


# read CNs of all cells in a run in RCK format with merged segments
get_cn_per_run_rck <- function(dir, type = "total"){
  subdirs = list.files(dir, pattern = "^c[0-9]")
  d_seg = data.frame()
  for(d in subdirs){
    # d = subdirs[1]
    # print(d)
    dir_rck = file.path(dir, d)
    fname = file.path(dir_rck, "rck.scnt.tsv")
    cns_orig = read_tsv(fname, show_col_types = F)
   
    cns_orig %>% rowwise() %>% mutate(cnA = as.numeric(strsplit(strsplit(strsplit(extra, "=")[[1]][2],":")[[1]][3],",")[[1]][1]), cnB = as.numeric(strsplit(strsplit(strsplit(extra, "=")[[1]][2],":")[[1]][4],"}")[[1]][1])) %>% mutate(total_cn = cnA + cnB) %>% as.data.frame() -> cns
    sid = str_replace(d, "c", "")
    # get division ID to show time progression
    fsv = file.path(dir, list.files(dir, pattern = paste0("SVData_c", sid, "_div.*.tsv")))
    bname = tools::file_path_sans_ext(basename(fsv))
    fields = str_split(bname, pattern = "_")[[1]]
    # sid = str_replace(bname, "CNData_", "")
    # sid = paste0(fields[[3]], fields[[2]], sep = "_")
    did = str_replace(fields[[3]], "div", "")    
    cns$sample = paste(sid, did, sep = "_")
    
    if(type == "total"){
      seg = cns %>% dplyr::select(sample, chrom = "chr", start, end, cn = total_cn)
    }else if(type == "A"){
      seg = cns %>% dplyr::select(sample, chrom = "chr", start, end, cn = cnA)
    }else{
      # type == "B"
      seg = cns %>% dplyr::select(sample, chrom = "chr", start, end, cn = cnB)
    }
    
    d_seg = rbind(d_seg, seg)
  }
  
  return(d_seg)
}


# get SV size (DUP/DEL)
get_sv <- function(bdir, n_dsbs, n_unrps, n_frag, cell_IDs, nrun){
  df_sv = data.frame()
  for(n_dsb in n_dsbs){
    # n_dsb = 50
    print(n_dsb)
    for(n_unrp in n_unrps){
      print(n_unrp)  
      for(i in 1:nrun){
        # print(i)
        if(n_frag > -1){
          dir = file.path(bdir, paste0("/nDSB", n_dsb, "_Un", n_unrp, "_frag", n_frag, "/r", i))         
        }else{
          dir = file.path(bdir, paste0("/nDSB", n_dsb, "_Un", n_unrp, "/r", i)) 
        }

        for(cell_ID in cell_IDs){
          midfix <- paste0("c", cell_ID, "_div", div_ID)
          
          # fname = "/Users/bl0033/data/SV/sim/test20220503/nDSB50_nUn0/r1/SVData_c2_div1.tsv"
          fname <- file.path(dir, paste0("SVData_", midfix, ".tsv"))
          if(!file.exists(fname)){
            stop("file not exist!")
          }
          df = read_tsv(fname, show_col_types = F)
          df$run = i
          df$cell = cell_ID
          df$nDSB = n_dsb
          df$nUnrep = n_unrp
          df_sv = rbind(df_sv, df)
        }
      }
    }
  }
  return(df_sv)
}

# read all CNs (total, A, B) of all cells in a run in matrix format: row -- cell, column -- segment
get_cn_mat <- function(dir){
  fcns = list.files(dir, pattern = "CNData*")
  fcns
  nf = length(fcns)
  # get CN matrix for all cells, each row represents a CN vector for a cell
  cn_mat = data.frame()
  cnA_mat = data.frame()
  cnB_mat = data.frame()
  for(i in 1:nf){
    # i = 1
    # print(i)
    fname = fcns[i]
    # extract cell ID 
    fields = str_split(tools::file_path_sans_ext(fname), "_")[[1]]
    # cell = str_replace(fields[2], "c", "")  
    cell = fields[2]
    fpath = file.path(dir, fname)
    cns_orig = read_tsv(fpath, show_col_types = F)
    cns_orig$cell = cell
    # print(dim(cns_orig))
    
    cns_orig %>% dplyr::select(cell, chromosome, start, end, total_cn) %>% spread(cell, total_cn) %>% dplyr::select(-c("chromosome", "start", "end")) %>% t() -> dcn
    cn_mat = rbind(cn_mat, dcn)
    
    cns_orig %>% dplyr::select(cell, chromosome, start, end, cnA) %>% spread(cell, cnA) %>% dplyr::select(-c("chromosome", "start", "end")) %>% t() -> dcnA
    cnA_mat = rbind(cnA_mat, dcnA)
    
    cns_orig %>% dplyr::select(cell, chromosome, start, end, cnB) %>% spread(cell, cnB) %>% dplyr::select(-c("chromosome", "start", "end")) %>% t() -> dcnB
    cnB_mat = rbind(cnB_mat, dcnB)  
  }
  
  return(list(cn_mat = cn_mat, cnA_mat = cnA_mat, cnB_mat = cnB_mat, cns_orig = cns_orig))
}


# find regions with the same CN != 2 across the genome
# ncutoff: maximum types of CNs in a column, fewer types, more likely to be clonal
get_var_seg <- function(cn_mat, cnA_mat, cnB_mat, ncutoff){
  nuniq = cn_mat %>% summarise(across(everything(), n_distinct)) 
  cuniq = names(nuniq[which(nuniq > 1 & nuniq <  ncutoff)])
  # haplotype A
  nuniqA = cnA_mat %>% summarise(across(everything(), n_distinct)) 
  cuniqA = names(nuniqA[which(nuniqA > 1 & nuniq < ncutoff)])
  # haplotype B
  nuniqB = cnB_mat %>% summarise(across(everything(), n_distinct)) 
  cuniqB = names(nuniqA[which(nuniqB > 1 & nuniq <  ncutoff) ])
  # find intersection
  luniq0 = intersect(cuniq, cuniqA)
  luniq = intersect(luniq0, cuniqB)
  
  return(luniq)
}


# check what these regions are on each haplotype
# diffAB = cnA_sel - cnB_sel
# for each loci, find parallel events, use count to sort
# find rows with the same TCN (index of the same value in a vector) but cnA != cnB
get_npar <- function(cn_sel, cnA_sel, cnB_sel){
  nloc = ncol(cn_sel)
  npars = c()
  for(i in 1:nloc){
    # i = 2
    np = 0
    cn = cn_sel[, i]
    cnA = cnA_sel[, i]
    cnB = cnB_sel[, i]
    ucn = setdiff(unique(cn), c(2))
    for(c in ucn){
      # c = 1
      # print(c)
      idx = which(cn == c)
      ucnA = cnA[idx]
      ucnB = cnB[idx]
      # check that cnA + cnB = cn
      eqls = (ucnA + ucnB == cn[idx])
      if(!all(eqls)){
        stop("cnA + cnB != cn")
      } 
      # count unique combinations of cnA and cnB
      ucnAB = data.frame(cnA = ucnA, cnB = ucnB)
      # print(ucnAB)
      ncmb = unique(ucnAB)
      # print(ncmb)
      if(nrow(ncmb) > 1) np = np + length(idx)
    }
    npars = c(npars, np)
  }
  names(npars) = colnames(cn_sel)
  sorted_npars = sort(npars)
  return(sorted_npars)
}


get_edge_type <- function(d_parA2t, stree){
  # one cell can have multiple types due to positional variety 
  stype = d_parA2t %>% group_by(sample) %>% summarise(type = Mode(type))
  # need to map node labels in the tree by renaming sampe in stype
  # edge a two-column matrix where each row represents a branch (or edge) of the tree; the nodes and the tips are symbolized with integers; the n tips are numbered from 1 to n, and the m (internal) nodes from n+1 to n+m (the root being n + 1). For each row, the first column gives the ancestor.
  # edge.length (optional) a numeric vector giving the lengths of the branches given by edge.
  # tip.label a vector of mode character giving the labels of the tips; the order of these labels corresponds to the integers 1 to n in edge.
  # Nnode an integer value giving the number of nodes in the tree (m). node.label (optional) a vector of mode character giving the labels of the
  # nodes (ordered in the same way than tip.label).
  tip_name = data.frame(tip = 1:length(stree$tip.label), sample = stree$tip.label)
  stype_cmb = merge(stype, tip_name, by = c("sample")) %>% arrange(tip) %>% select(X2 = tip, type)
  edges = data.frame(stree$edge, eid = 1:nrow(stree$edge))
  # %>% mutate(type = if_else(is.na(type), "None", type))
  edge_type = merge(edges, stype_cmb, by = c("X2"), all.x = T) %>% select(node = X2, type, eid) %>% arrange(eid)
  
  return(edge_type)
}


get_subtree <- function(dir, cell_par){
  # get whole tree
  ftree = file.path(dir, "cell_lineage.tsv")
  
  dtree = read_tsv(ftree, show_col_types = F)
  ig = graph_from_edgelist(as.matrix(dtree[, c(1, 2)]), directed = TRUE)
  mytree = as.phylo(ig)
  
  # title = " "
  # p = ggtree(mytree) + ggtitle(title) + geom_text(aes(label=label), hjust=-.3)
  # extract selected cells
  
  tips = as.character(cell_par)
  stree = keep.tip(mytree, tips)
  
  return(stree)
}

############ other functions ############ 
# include events in one run of the simulation
get_sumstats_sim <- function(bdir, n_dsbs, f_unrps, frags_local, nrun, model = 0, dsuffix = ""){
  df_sstat = data.frame()
  
  for(n_dsb in n_dsbs){
    # n_dsb = 50
    # print(n_dsb)
    for(f_unrp in f_unrps){
      # print(f_unrp)  
      for(frag in frags_local){
        # print(f_unrp)       
        for(i in 1:nrun){
          # print(i)
          dprefix = paste0("nDSB", n_dsb, "_Un", f_unrp, "_frag", frag)
          if(model == 1){
            dprefix = paste0(dprefix, "_model", model)
          }
          dprefix = paste0(dprefix, dsuffix)
          dir = file.path(bdir, dprefix, paste0("r", i))
          
          fname <- file.path(dir, paste0("sumStats_sim.tsv"))
          if(!file.exists(fname)) next
          
          df = read_tsv(fname, show_col_types = F)
          df$f_unrp = f_unrp
          df$nfrag = frag
          df$run = i
          
          df_sstat = rbind(df_sstat, df)          
        }
      }
    }
  }
  
  return(df_sstat)
}


get_sumstats_total_by_dir <- function(dir){
  df_sstat = data.frame()
  fnames = list.files(dir, pattern = "sumStats_total_*")
  # for(cell_ID in cell_IDs){
  # print(cell_ID)
  for(f in fnames){
    # midfix <- paste0("c", cell_ID, "_div", div_ID)
    # fname <- file.path(dir, paste0("sumStats_total_", midfix, ".tsv"))
    fname <- file.path(dir, f)
    # print(fname)
    if(!file.exists(fname)) next
    
    df = read_tsv(fname, show_col_types = F)

    # compute ratio of DUP, DEL, INV,
    df$ratio1 = (df$nDup) / (df$nDel)
    df$ratio2 = (df$nH2HInv + df$nT2TInv) / (df$nDel)
    # # add 1 to avoid division by 0?
    # df$ratio1 = (df$nDup + 1) / (df$nDel + 1)
    # df$ratio2 = (df$nH2HInv + df$nT2TInv + 2) / (df$nDel + 1)
    
    df_sstat = rbind(df_sstat, df)          
  }  
  return(df_sstat)
}


# include events in each cell
# hard to specify cell_IDs and div_ID when there are many cells, so extracting from file names
get_sumstats_total <- function(bdir, n_dsbs, f_unrps, frags_local, nrun, model, model_suffix = F, dsuffix = ""){
  df_sstat = data.frame()
  
  for(n_dsb in n_dsbs){
    # n_dsb = 50
    # print(n_dsb)
    for(f_unrp in f_unrps){
      # print(f_unrp)  
      for(frag in frags_local){
        # print(f_unrp)       
        for(i in 1:nrun){
          # print(i)
          dprefix = paste0("nDSB", n_dsb, "_Un", f_unrp, "_frag", frag)
          if(model_suffix){
            dprefix = paste0(dprefix, "_model", model)
          }
          dprefix = paste0(dprefix, dsuffix)
          dir = file.path(bdir, dprefix, paste0("r", i))
          print(dir)
          df = get_sumstats_total_by_dir(dir)
          df$f_unrp = f_unrp
          df$nfrag = frag
          df$run = i
          df_sstat = rbind(df_sstat, df)     
        }
      }
    }
  }

  return(df_sstat)
}



# summary for each cell
get_sim_stat <- function(bdir, n_dsbs, n_unrps, cell_IDs, nrun){
  df_sstat = data.frame()
  for(n_dsb in n_dsbs){
    # n_dsb = 50
    print(n_dsb)
    for(n_unrp in n_unrps){
      print(n_unrp)  
      for(i in 1:nrun){
        # print(i)
        dir = file.path(bdir, paste0("/nDSB", n_dsb, "_Un", n_unrp, "/r", i)) 
        for(cell_ID in cell_IDs){
          midfix <- paste0("c", cell_ID, "_div", div_ID)
          
          fname <- file.path(dir, paste0("sumStats_total_", midfix, ".tsv"))
          if(!file.exists(fname)) next
          df = read_tsv(fname, show_col_types = F)
          df$run = i
          # compute ratio of DUP, DEL, INV,
          df$ratio1 = (df$nDup) / (df$nDel)
          df$ratio2 = (df$nH2HInv + df$nT2TInv) / (df$nDel)
          # # add 1 to avoid division by 0?
          # df$ratio1 = (df$nDup + 1) / (df$nDel + 1)
          # df$ratio2 = (df$nH2HInv + df$nT2TInv + 2) / (df$nDel + 1)
          
          df_sstat = rbind(df_sstat, df)
        }
      }
    }
  }
  return(df_sstat)
}


# read previous result to save time when reading many files
get_sumstats_total_wrap <- function(bdir1, n_dsbs, f_unrps, frags_local, nrun){
  f1 = file.path(bdir1, "df_sstat_total.tsv")
  if(!file.exists(f1)){
    df_sstat1 = get_sumstats_total(bdir1, n_dsbs, f_unrps, frags_local, nrun)
    write_tsv(df_sstat1, f1)
  }else{
    df_sstat1 = read_tsv(f1)
  }
  return(df_sstat1)
}


get_ecdna_all <- function(bdir1, bdir2, bdir3, t1, t2, t3, n_dsbs, f_unrps, frags_local, nrun){
  df_ecnda_all = data.frame()
  
  # 3600 records for each case
  # div_ID = 1
  # cell_IDs = c(2,3)
  df_sstat1 = get_sumstats_total_wrap(bdir1, n_dsbs, f_unrps, frags_local, nrun)
  df_sstat1$type = t1
  df_ecnda_all = rbind(df_ecnda_all, df_sstat1)
  
  # div_ID = 2
  # cell_IDs = c(4,5)
  df_sstat2 = get_sumstats_total_wrap(bdir2, n_dsbs, f_unrps, frags_local, nrun)
  df_sstat2$type = t2
  df_ecnda_all = rbind(df_ecnda_all, df_sstat2)
  
  # div_ID = 2
  # cell_IDs = c(4,5)
  df_sstat3 = get_sumstats_total_wrap(bdir3, n_dsbs, f_unrps, frags_local, nrun)
  df_sstat3$type = t3
  df_ecnda_all = rbind(df_ecnda_all, df_sstat3)
  
  return(df_ecnda_all)
}


get_sumstats_sim_wrap <- function(bdir1, n_dsbs, f_unrps, frags_local, nrun){
  f1 = file.path(bdir1, "df_sstat_sim.tsv")
  if(!file.exists(f1)){
    df_bfb1 = get_sumstats_sim(bdir1, n_dsbs, f_unrps, frags_local, nrun)
    write_tsv(df_bfb1, f1)
  }else{
    df_bfb1 = read_tsv(f1)
  }
  return(df_bfb1)  
}


get_bfb_all <- function(bdir1, bdir2, bdir3, t1, t2, t3, n_dsbs, f_unrps, frags_local, nrun){
  df_bfb_all = data.frame()
  
  df_bfb1 = get_sumstats_sim_wrap(bdir1, n_dsbs, f_unrps, frags_local, nrun)
  df_bfb1$type = t1
  df_bfb_all = rbind(df_bfb_all, df_bfb1)
  
  df_bfb2 = get_sumstats_sim_wrap(bdir2, n_dsbs, f_unrps, frags_local, nrun)
  df_bfb2$type = t2
  df_bfb_all = rbind(df_bfb_all, df_bfb2)
  
  df_bfb3 = get_sumstats_sim_wrap(bdir3, n_dsbs, f_unrps, frags_local, nrun)
  df_bfb3$type = t3
  df_bfb_all = rbind(df_bfb_all, df_bfb3)
  
  return(df_bfb_all)
}


# summary for each cell division
# nCell	nComplex	nMbreak	nTelofusion
get_stat_cell <- function(bdir, n_dsbs, n_unrps, cell_IDs, nrun){
  df_sstat_sim = data.frame()
  for(n_dsb in n_dsbs){
    # n_dsb = 50
    print(n_dsb)
    for(n_unrp in n_unrps){
      print(n_unrp)  
      for(i in 1:nrun){
        # n_dsb = 50
        # n_unrp = 5
        # i = 1
        print(i)
        dir = file.path(bdir, paste0("/nDSB", n_dsb, "_Un", n_unrp, "/r", i)) 
        # fname = "/Users/bl0033/data/SV/sim/test20220503/nDSB50_nUn0/r1/shatterseek_c2_div1.tsv"
        fname <- file.path(dir, paste0("sumStats_sim.tsv"))
        df = read_tsv(fname, show_col_types = F)
        df$run = i
        df_sstat_sim = rbind(df_sstat_sim, df)
      }
    }
  } 
  return(df_sstat_sim)
}


# read summary statistics for all cells in a run
get_sumstats_total_by_dir <- function(dir){
  df_sstat = data.frame()
  sfiles = list.files(dir, pattern = "sumStats_total_.*.tsv")
  for(i in 1:length(sfiles)){
    fname = file.path(dir, sfiles[i])
    # print(fname)
    if(!file.exists(fname)) next
    
    # extract division ID
    fields = str_split(tools::file_path_sans_ext(fname), "_")[[1]]
    sdiv = fields[length(fields)]
    divID = str_replace(sdiv, "div", "")
    
    df = read_tsv(fname, show_col_types = F)
    df$division = as.numeric(divID)
    # compute ratio of DUP, DEL, INV,
    df$ratio1 = (df$nDup) / (df$nDel)
    df$ratio2 = (df$nH2HInv + df$nT2TInv) / (df$nDel)
    # # add 1 to avoid division by 0?
    # df$ratio1 = (df$nDup + 1) / (df$nDel + 1)
    # df$ratio2 = (df$nH2HInv + df$nT2TInv + 2) / (df$nDel + 1)
    
    df_sstat = rbind(df_sstat, df)          
  }
  
  # check selection
  fsel = file.path(dir, "selection.tsv")
  # print(fsel)
  if(file.exists(fsel)){
    sel = read.table(fsel, header = T) %>% dplyr::rename(cellID = cell)
    cell_pos = sel %>% filter(selection_coef > 0)
    cell_neg = sel %>% filter(selection_coef < 0)
    if(nrow(cell_pos) + nrow(cell_neg) != nrow(sel)){
      stop("no selection under selection model!")
    }
    df_all = merge(df_sstat, sel, by = "cellID") 
    df_all = df_all %>% rowwise() %>% mutate(model = ifelse(cellID %in% cell_pos$cell, "positive selection", "negative selection"))    
  }else{
    df = data.frame(prob_survival = 1, selection_coef = 0, model = "neutral")
    df_all = cbind(df_sstat, df)
  }

  return(df_all)
}


# get overlap with cosmic
get_gene_overlap <- function(csep, d_segT, chrsel = "7"){
  csep %>% filter(chr == chrsel) %>% filter(!is.na(disease)) -> csep_sel
  d_segT  %>% filter(chrom == chrsel) %>% filter(cn != 2) -> dseg_alt
  # max(dseg_alt$end)
  # max(csep_sel$end)
  # find overlap
  # csep_sel %>% filter(end < dend) %>% filter(grepl("TSG", onc) | grepl("oncogene", onc)) -> ccg_olp
  target.gr = GRanges(seqnames = csep_sel$chr, ranges = IRanges(start = csep_sel$start, end = csep_sel$end), names = csep_sel$gname, desc = csep_sel$onc, disease = csep_sel$disease)
  query.gr = GRanges(seqnames = dseg_alt$chrom, ranges = IRanges(start = dseg_alt$start, end = dseg_alt$end), scores = dseg_alt$cn)
  nolp = sum(countOverlaps(query.gr, target.gr))
  
  gene_overlap = data.frame()
  if(nolp > 0){
    overlaps = findOverlaps(query.gr, target.gr)
    # print(overlaps)
    olp_query = dseg_alt[overlaps@from,]
    olp_target = csep_sel[overlaps@to,]
    df = bind_cols(olp_query, olp_target)
    gene_overlap = rbind(gene_overlap, df)
  }
  
  return(gene_overlap)
}




# count LOH/MASI events
count_loh_masi <- function(bdir, mu, dt, nseg){
  idir = file.path(bdir, paste0("run-5-", mu, "-", mu, "-", dt, "-mode1-size1-", nseg, "-1-6"))
  df_count = data.frame()
  
  for(i in 1:nsim){
    # i = 1
    # print(i)
    
    fsim = file.path(idir, paste0("sim-data-", i, "-allele-cn.txt.gz"))
    fsim_total = file.path(idir, paste0("sim-data-", i, "-cn.txt.gz"))
    
    cns_wide = get_cns_wide(fsim, T)
   
    nloh = cns_wide %>% filter_all(any_vars(grepl("2_0|0_2" , .))) %>% nrow()
    
    cns_wide$rname = rownames(cns_wide)
    
    cns_wide_total = get_cns_wide(fsim_total, F)
    max_cn = max(cns_wide_total)
    
    cns_wide_total$rname = rownames(cns_wide_total)
    
    # find locations with variable CNs on each row 
    # cns_wide_total %>% rowwise() %>% filter(sum(c_across(where(is.numeric))) != 2 * ncol(cns_wide_total))
    cns_var = cns_wide_total %>% rowwise() %>% filter(n_distinct(c_across(where(is.numeric))) > 1) 
    # check allele-specific CNs at those sites 
    cns_var_allele = cns_wide %>% filter(rname %in% cns_var$rname) 
    
    cns_var_mat = cns_var  %>% dplyr::select(-rname) %>% as.matrix()
    cns_var_allele_mat  = cns_var_allele  %>% dplyr::select(-rname) %>% as.matrix()
    mat_all <- matrix(paste(cns_var_mat,cns_var_allele_mat,sep="-"), dim(cns_var_mat))
    
    var_all = mat_all %>% as.data.frame() %>% rowwise() %>% filter(n_distinct(c_across()) > 2) 
    masi1 = var_all %>% filter_all(any_vars(grepl("1-1_0" , .))) %>% filter_all(any_vars(grepl("1-0_1" , .))) 
    masi2 = var_all %>% filter_all(any_vars(grepl("2-2_0" , .))) %>% filter_all(any_vars(grepl("2-0_2" , .))) 
    masi3 = var_all %>% filter_all(any_vars(grepl("3-2_1" , .))) %>% filter_all(any_vars(grepl("3-1_2" , .))) 
    masi41 = var_all %>% filter_all(any_vars(grepl("4-3_1" , .))) %>% filter_all(any_vars(grepl("4-1_3" , .))) 
    masi42 = var_all %>% filter_all(any_vars(grepl("4-4_0" , .))) %>% filter_all(any_vars(grepl("4-0_4" , .))) 
    masi51 = var_all %>% filter_all(any_vars(grepl("5-4_1" , .))) %>% filter_all(any_vars(grepl("5-1_4" , .))) 
    masi52 = var_all %>% filter_all(any_vars(grepl("5-3_2" , .))) %>% filter_all(any_vars(grepl("5-2_3" , .))) 
    masi53 = var_all %>% filter_all(any_vars(grepl("5-5_0" , .))) %>% filter_all(any_vars(grepl("5-0_5" , .))) 
    masi61 = var_all %>% filter_all(any_vars(grepl("6-5_1" , .))) %>% filter_all(any_vars(grepl("6-1_5" , .))) 
    masi62 = var_all %>% filter_all(any_vars(grepl("6-4_2" , .))) %>% filter_all(any_vars(grepl("6-2_4" , .))) 
    masi63 = var_all %>% filter_all(any_vars(grepl("6-6_0" , .))) %>% filter_all(any_vars(grepl("6-0_6" , .))) 

    masi = bind_rows(masi1, masi2, masi3, masi41, masi42, masi51, masi52, masi53, masi61, masi62, masi63)
    # for each row, check the values to count MASI events
    nmasi = nrow(masi)
    
    df = data.frame(run = i, nloh = nloh, nmasi = nmasi, mu = as.factor(mu), nseg = as.factor(nseg), max_cn = max_cn)
    df_count = rbind(df_count, df)
  }
  
  return(df_count)
}
