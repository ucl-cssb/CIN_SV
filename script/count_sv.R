#!/usr/bin/env Rscript

# This script is used to count different types of SVs generated in simulated data (Fig 3)

# change path accordingly
source("parse_sim_util.R")

########## basic setting ########## 
# n_dsbs = c(10, 50)
# f_unrps = c(0, 0.1, 0.5)
# frags_local = c(0, 10, 50)
n_dsbs = c(10, 30)
f_unrps = c(0, 0.1, 0.3)
frags_local = c(0, 10, 30)

nrun = 50
chrs_sel = c("1")

 # change path accordingly
odir = "data_svmodel/Fig3"
rdir = "data_svmodel/Fig3"
bdir = "data_svmodel/Fig3" 

t1 = "1 cycle"
t2 = "2 cycles (interphase DSB at 1st cycle)"
t3 = "2 cycles"
tlvl = c(t1, t2, t3)


bdir1 = file.path(bdir, "test_ncell2_break0")
bdir2 = file.path(bdir, "test_ncell3_break0")
bdir3 = file.path(bdir, "test_ncell3_break1")

# plot size, consistent among panels
w = 8
h = 4

xll = "mean number of DSBs in local fragmentation"
ylb_ecdna = "#ecDNAs in a cell"
# ylb_ecdna = "number of ecDNAs in a cell"
ylb_sm = "#cells with seismic amplifications"

########## count chromothripsis (only consider chr1) under neutral model ########## 
# 18 (params) * 50 (runs) * 2 (cells) * 3 (sets) = 5400
df_chromth_all = get_chromth_all(bdir1, bdir2, bdir3, t1, t2, t3, n_dsbs, f_unrps, frags_local, nrun, chrs_sel)
unique(df_chromth_all$type)
df_chromth_all$type = factor(df_chromth_all$type, levels = tlvl)

# take only two cells for each run
df_chromth_sel = df_chromth_all %>% filter(!(type != t1 & div_ID == 1))
dim(df_chromth_sel)
df_chromth_sel %>% group_by(cell_ID) %>% tally()
# 18 * 2 * 3
df_chromth_sel %>% select(nDSB, nUnrepair, nfrag, type, cell_ID) %>% unique() -> grps
dim(grps)

# rename cell IDs to facilitate counting
df_chromth_rname = df_chromth_sel %>% mutate(cell_ID = case_when(cell_ID == 4 ~ 2, cell_ID == 5 ~ 3, .default = cell_ID))
df_chromth_rname %>% group_by(cell_ID) %>% tally()
df_chromth_rname %>% group_by(type) %>% tally()

# %>% filter(nhconf > 0 | nlconf > 0), only consider high confidence calls to reduce FPs
# ignore those with n = 0 to avoid summarizing issues
df_chromth_count = df_chromth_rname %>% filter(nhconf > 0) %>% group_by(nDSB, nUnrepair, nfrag, run, cell_ID, type, .drop = T) %>% tally() 
# only count presence of chromothripsis, 684 / 5400 with detected events
df_chromth_count %>% mutate(nchrmth = ifelse(n >= 1, 1, 0)) %>% select(-n) -> sum_per_cell
dim(sum_per_cell)
unique(sum_per_cell$nchrmth)

# not use percentage, a bit confusing
sum_per_cell %>% group_by(nDSB, nUnrepair, nfrag, type, .drop = T) %>% count(nchrmth) %>% mutate(perc = n / 100) -> summ
summ_all = merge(grps, summ, by = c("nDSB", "nUnrepair", "nfrag", "type"), all.x = T)

pct = ggplot(summ_all, aes(y = n, x=as.factor(nfrag), fill = as.factor(nUnrepair))) + geom_bar(stat = "identity", position = position_dodge()) + geom_hline(aes(yintercept = 10), linetype = "dashed", color = "darkgrey") + theme_pubr() + ylab("#cells with chromothripsis on chr1") + facet_grid(as.factor(nDSB) ~ as.factor(type)) + xlab(xll) + guides(fill = guide_legend(lbl_unrepair)) + scale_y_continuous(n.breaks = 5) + scale_fill_manual(values = unrep_colors) 

pct2 = ggplot(summ_all, aes(y = n, x=as.factor(nfrag), fill = as.factor(type))) + geom_bar(stat = "identity", position = position_dodge()) + geom_hline(aes(yintercept = 10), linetype = "dashed", color = "darkgrey") + theme_pubr() + ylab("#cells with chromothripsis on chr1") + facet_grid(as.factor(nDSB) ~ as.factor(nUnrepair)) + xlab(xll) + guides(fill = guide_legend("")) + scale_y_continuous(n.breaks = 5) + scale_fill_manual(values = type_colors) 



########## check circles ########## 
df_ecnda_all = get_ecdna_all(bdir1, bdir2, bdir3, t1, t2, t3, n_dsbs, f_unrps, frags_local, nrun)
df_ecnda_all$type = factor(df_ecnda_all$type, levels = tlvl)

df_ecnda_sel = df_ecnda_all %>% filter(!(type != t1 & cycleID == 1))

# 2 cells for run, 100 cells for each case
df_ecnda_sel %>% group_by(nDSB, f_unrp, nfrag, type, cycleID) %>% tally() %>% View()
df_ecnda_sel %>% group_by(type) %>% tally()

pec = ggplot(df_ecnda_sel, aes(y = (nECDNA), x=as.factor(nfrag), fill = as.factor(f_unrp))) + geom_boxplot() + theme_pubr() + ylab(ylb_ecdna) + facet_grid(as.factor(nDSB) ~ as.factor(type)) + xlab(xll) + guides(fill = guide_legend(lbl_unrepair)) + scale_fill_manual(values = unrep_colors) + geom_hline(aes(yintercept = 20), linetype = "dashed", color = "darkgrey") 

pec2 = ggplot(df_ecnda_sel, aes(y = (nECDNA), x=as.factor(nfrag), fill = as.factor(type))) + geom_boxplot() + theme_pubr() + ylab(ylb_ecdna) + facet_grid(as.factor(nDSB) ~ as.factor(f_unrp)) + xlab(xll) + guides(fill = guide_legend("")) + scale_fill_manual(values = type_colors) + geom_hline(aes(yintercept = 20), linetype = "dashed", color = "darkgrey") 



########## check BFBs ########## 
df_bfb_all = get_bfb_all(bdir1, bdir2, bdir3, t1, t2, t3, n_dsbs, f_unrps, frags_local, nrun)
# by run, no need to exclude cells
df_bfb_all$type = factor(df_bfb_all$type, levels = tlvl)

df_bfb_all %>% group_by(nDSB, nUnrepair, f_unrp, nfrag, type) %>% tally() %>% select(n) %>% unique()

pbfb = ggplot(df_bfb_all, aes(y = nTelofusion, x=as.factor(nfrag), fill = as.factor(f_unrp))) + geom_boxplot() + theme_pubr() + ylab(lbl_fusion) + facet_grid(as.factor(nDSB) ~ as.factor(type)) + xlab(xll) + guides(fill = guide_legend(lbl_unrepair)) + scale_fill_manual(values = unrep_colors) + geom_hline(aes(yintercept = 25), linetype = "dashed", color = "darkgrey") 

pbfb2 = ggplot(df_bfb_all, aes(y = nTelofusion, x=as.factor(nfrag), fill = as.factor(type))) + geom_boxplot() + theme_pubr() + ylab(lbl_fusion) + facet_grid(as.factor(nDSB) ~ as.factor(f_unrp)) + xlab(xll) + guides(fill = guide_legend("")) + scale_fill_manual(values = type_colors) + geom_hline(aes(yintercept = 25), linetype = "dashed", color = "darkgrey") 


########## count Seismic Amplification ############ 
chrsel = "chr1"
plot = T

fout = file.path(odir, "nsm_count.tsv")

if(!file.exists(fout)){
  df_sm = get_df_sm_batch(bdir, chrsel, plot, 0)
  write_tsv(df_sm, fout)
}else{
  df_sm = read_tsv(fout)
}

df_sm %>% group_by(nfile) %>% tally()
df_sm %>% group_by(ncell, div_break) %>% tally()

df_smn = df_sm %>% group_by(ncell, div_break, ndsb, funrep, local_frag) %>% summarize(n = sum(nsm > 0))

# find cases with Seismic Amplification
df_with_sm = df_sm %>% filter(nsm > 0) 

df_with_smn = df_with_sm %>% group_by(ncell, div_break, ndsb, funrep, local_frag) %>% tally() 

df_smn_annot = df_smn %>% mutate(type = t1) %>% mutate(type = case_when(ncell == 2 ~ t1, ncell == 3 & div_break == 0 ~ t2, .default = t3))
df_smn_annot$type = factor(df_smn_annot$type, levels = tlvl)
df_smn_annot %>% group_by(type) %>% tally()


psm = ggplot(df_smn_annot, aes(y = n, x=as.factor(local_frag), fill = as.factor(funrep))) + geom_bar(stat = "identity", position = position_dodge()) + theme_pubr() + ylab(ylb_sm) + facet_grid(as.factor(ndsb) ~ as.factor(type)) + xlab(xll) + guides(fill = guide_legend(lbl_unrepair)) + scale_fill_manual(values = unrep_colors) + geom_hline(aes(yintercept = 2), linetype = "dashed", color = "darkgrey") 

psm2 = ggplot(df_smn_annot, aes(y = n, x=as.factor(local_frag), fill = as.factor(type))) + geom_bar(stat = "identity", position = position_dodge()) + theme_pubr() + ylab(ylb_sm) + facet_grid(as.factor(ndsb) ~ as.factor(funrep)) + xlab(xll) + guides(fill = guide_legend("")) + scale_fill_manual(values = type_colors) + geom_hline(aes(yintercept = 2), linetype = "dashed", color = "darkgrey") 


############ combine the plots #############

p4 = ggarrange(pbfb, pct, pec, psm, nrow = 2, ncol = 2, common.legend = T)
fplot = file.path(rdir, "sim_summary.pdf")
ggsave(fplot, p4, width = w * 2.1, height = h * 2.5)


p4 = ggarrange(pbfb2, pct2, pec2, psm2, nrow = 2, ncol = 2, common.legend = T)
fplot = file.path(rdir, "sim_summary_type.pdf")
ggsave(fplot, p4, width = w * 2.1, height = h * 2.5)






