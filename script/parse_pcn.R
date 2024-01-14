
#!/usr/bin/env Rscript

# This script is used to parse the simulated data with parallel evolution of CNAs (Fig 2)

# change path accordingly
source("parse_sim_util.R")


################## basic settings to read input ################## 
# use this for paper, same params as chromoplexy simulations
dir = "data_svmodel/Fig2/chromplexy/r3"

# read all CNs, should have the same dimension
# parse shatterseek data, as the breakpoints are aligned 
# taking lots of memory for 1000 cells
d_seg = get_cn_per_run(dir, type = "total")
d_segA = get_cn_per_run(dir, type = "A")
d_segB = get_cn_per_run(dir, type = "B")

res_mat = get_cn_mat(dir)  # take some time
cn_mat = res_mat$cn_mat
cnA_mat = res_mat$cnA_mat
cnB_mat = res_mat$cnB_mat
cns_orig = res_mat$cns_orig   # original value for one cell, used to get exact location
dim(cn_mat)
ncell = nrow(cn_mat)
ncell

# To ensure there are just a few values for clarity
# ncutoff = ncell / 30
ncutoff = 6 
luniq = get_var_seg(cn_mat, cnA_mat, cnB_mat, ncutoff)


cn_sel = cn_mat[,luniq]
cnA_sel = cnA_mat[,luniq]
cnB_sel = cnB_mat[,luniq]
sorted_npars = get_npar(cn_sel, cnA_sel, cnB_sel)

# plot regions with more parallel events
idx_par = names(sorted_npars[which(sorted_npars > 0)])
cn_par = cn_sel[idx_par]
cnA_par = cnA_sel[idx_par]
cnB_par = cnB_sel[idx_par]

# find the locations to plot (from one cell, same across all cells)
cns_orig = cns_orig %>% as.data.frame()
rownames(cns_orig) = paste0("V", seq(1:nrow(cns_orig)))
loc_sel = cns_orig[idx_par, ] %>% select(chrom = chromosome, start, end)

# get local regions as the pattern is not clear on a whole chromosome
d_par = merge(loc_sel, d_seg, by = c("chrom", "start", "end"), all.x = T)
d_par %>% group_by(sample) %>% tally()
d_parA = merge(loc_sel, d_segA, by = c("chrom", "start", "end"), all.x = T)
d_parB = merge(loc_sel, d_segB, by = c("chrom", "start", "end"), all.x = T)
d_parAB = merge(d_parA, d_parB, by = c("sample", "chrom", "start", "end")) 


chr_sel = 22
cell_sel = d_parAB %>% filter(chrom == chr_sel) %>% filter(!(cn.x == 1 & cn.y == 1)) %>% select(sample) %>% unique() %>% unlist()
length(cell_sel)

d_par1 = d_par  %>% filter(chrom == chr_sel) %>% filter(sample %in% cell_sel)
d_parA1 = d_parA  %>% filter(chrom == chr_sel) %>% filter(sample %in% cell_sel)
d_parB1 = d_parB  %>% filter(chrom == chr_sel) %>% filter(sample %in% cell_sel)


################## plot as in Nature paper for manuscript  ################## 
# get lineage tree of selected cells  
cell_par = unique(d_par1$sample)
stree = get_subtree(dir, cell_par)

# get levels to maintain consistency with CN heatmap
d = fortify(stree)
slevel = d$label[order(d$y, decreasing = T)]
d_par1$sample = factor(d_par1$sample, levels = slevel)
d_parA1$sample = factor(d_parA1$sample, levels = slevel)

# check positions on focused chr
loc_sel_chr = loc_sel %>% filter(chrom == chr_sel)
# merge regions to confirm 
loc_range = GRanges(loc_sel_chr$chrom, IRanges(loc_sel_chr$start, loc_sel_chr$end))
loc_merged = reduce(loc_range)
ranges(loc_merged)

# plot whole regions, will be gap in between
p1 = plot.cn.heatmap(d_par1, "", theme = theme0)
p2 = plot.cn.heatmap(d_parA1, "", theme = theme0)
ggarrange(p1, p2, nrow = 1, ncol = 2)

# select largest regions on the chr when there are >1 regions
# for chr22
loc_cutoff = 27298459
loc_cutoff2 = 19444803
d_par2 = d_par1 %>% filter(end < loc_cutoff & start > loc_cutoff2)
d_parA2 = d_parA1 %>% filter(end < loc_cutoff & start > loc_cutoff2)

d_par2 = d_par1 
d_parA2 = d_parA1

# use a different color code for HCN for clarity
# For A CN, 0: B Hom if T > 0, 1: A Hom if T = 1, balanced if T = 2, 
# need all CNs to decide type
# d_ctype = d_parAB %>% mutate(type = case_when(cn.x == cn.y ~ "Balanced", cn.x > cn.y & cn.y == 0 ~ "A (hom)", cn.x > cn.y & cn.y > 0 ~ "A (gained)", cn.y > cn.x & cn.x == 0 ~ "B (hom)", cn.y > cn.x & cn.x > 0 ~ "B (gained)")) %>% select(-c("cn.x", "cn.y"))
d_ctype = d_parAB %>% mutate(type = case_when(cn.x == cn.y ~ hcn_types[5], cn.x > cn.y & cn.y == 0 ~ hcn_types[2], cn.x > cn.y & cn.y > 0 ~ hcn_types[1], cn.y > cn.x & cn.x == 0 ~ hcn_types[4], cn.y > cn.x & cn.x > 0 ~ hcn_types[3])) %>% select(-c("cn.x", "cn.y"))
d_parA2t = merge(d_parA2, d_ctype, by = c("sample", "chrom", "start", "end"))

edge_type = get_edge_type(d_parA2t, stree)
edge_colors = c(hcn_colors)
edge_types = c(hcn_types)
# node refers to all nodes, including root without incoming edges
stree_full = full_join(stree, edge_type, by = "node")
# when forgeting to set title, Error in as.character(x$label) : cannot coerce type 'closure' to vector of type 'character'
title = ""
# weird colors when using more colors and values than those in data
pt = ggtree(stree_full, aes(color = type)) + scale_color_manual(values = edge_colors, labels = edge_types) + ggtitle(title) + theme(legend.position = "none")
# pt = ggtree(stree) %<+% edge_type + aes(color = type) + scale_color_manual(values = edge_colors, labels = edge_types) + ggtitle(title) + theme(legend.position = "none")
# + geom_tiplab()

# for cases when types are not complete (selection model for r1)
edge_type %>% group_by(type) %>% tally()
my_types = sort(unique(edge_type$type))
my_colors = c("#1b7837", "#c2a5cf", "#d5d5d4")
pt = ggtree(stree_full, aes(color = type)) + scale_color_manual(values = my_colors, labels = my_types) + ggtitle(title) + theme(legend.position = "none")

p1 = plot.cn.heatmap(d_par2, "", theme = theme_notxt)
p2 = plot.hcn.heatmap(d_parA2t, "", hcn_colors, hcn_types, theme = theme_notxt)

# pdemo = ggarrange(pt, p1, p2, nrow = 1, ncol = 3).  # heatmap bottom slightly lower
pdemo = pt | p1 | p2
pdemo
suffix = ""
fout = file.path(dir, paste0("cn_par_demo_chr", chr_sel, suffix, ".pdf"))
ggsave(fout, pdemo, width = 3.7, height = 5)



