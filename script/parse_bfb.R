#!/usr/bin/env Rscript

# This script is used to parse data with one simulated DSB and simple/complex break (Fig 2)

# change path accordingly
source("parse_sim_util.R")

################# read data and plot all cells for simulations with one unrepaired DSB ################# 
odir = "data_svmodel/Fig2"
odir_break = file.path(odir, "simple_break")
odir_frag = file.path(odir, "local_fragmentation")

######## simple break
chrs = c("7")
chr = 7
dir1 = file.path(odir_break, "neutral")
plot_lineage_with_cn_sel(dir1, chrs, odir_break, "_model0", F,  "neutral model", F, theme2, 5, 4, T)
plot_lineage_with_cn_sel(dir1, chrs, odir_break, "_model0_lbl",F,  "neutral model", F, theme1, 8, 7, T)
plot_cn_sorted(odir_break, chr, 0, dir1, theme0, "", 2.55, 6) # plot one chromosome

dir2 = file.path(odir_break, "selection")
plot_lineage_with_cn_sel(dir2, chrs, odir_break, "_model1", F, "selection model", F, theme2, 5, 4, T)
plot_lineage_with_cn_sel(dir2, chrs, odir_break, "_model1_lbl", F, "selection model", F, theme1, 8, 7, T)
plot_cn_sorted(odir_break, chr, 1, dir2, theme0, "", 2.55, 6)


######## local fragmentation
chrs = c("2")
chr = 2
dir1_frag = file.path(odir_frag, "neutral")
plot_lineage_with_cn_sel(dir1_frag, chrs, odir_frag, "_model0_frag", T, "neutral model", F, theme2, 5, 4, F)
plot_lineage_with_cn_sel(dir1_frag, chrs, odir_frag, "_model0_frag_lbl", T, "neutral model", F, theme1, 9, 8, F)
plot_cn_sorted(odir_frag, chr, 0, dir1_frag, theme0, "_frag50")

dir2_frag = file.path(odir_frag, "selection")
plot_lineage_with_cn_sel(dir2_frag, chrs, odir_frag, "_model1_frag", T, "selection model", F, theme2, 5, 4, F)
plot_lineage_with_cn_sel(dir2_frag, chrs, odir_frag, "_model1_frag_lbl", T, "selection model", F, theme1, 9, 8, F)
plot_cn_sorted(odir_frag, chr, 1, dir2_frag, theme0, "_frag50")


############ check CN plot for specific cells as in sequenza ############ 
subdir = "simple_break"
chrsel = "chr7"

# neutral
model = 0
cn_dir = dir1
suffix = paste0("_model", model)
codir = file.path(odir, subdir)
# select cells with highest CN under each model
plot_cell_with_maxcn(cn_dir, codir, suffix, chrsel, 1, T, 100, 20)

# selection
model = 1
cn_dir = dir2
suffix = paste0("_model", model)
codir = file.path(odir, subdir)
# select cells with highest CN under each model
plot_cell_with_maxcn(cn_dir, codir, suffix, chrsel, 1, T, 100, 5)


subdir = "local_frag"
chrsel = "chr2"

# neutral
model = 0
cn_dir = dir1_frag
suffix = paste0("_model", model)
plot_cell_with_maxcn(cn_dir, odir, subdir, suffix, chrsel, 3)

# for cell with seismic events
fcn = "CNData_c85_div5.tsv"
chrsel = "chr2"
plot_cell_cn_chr(cn_dir, fcn, odir_frag, suffix, chrsel, xgap = 1000)

# selection
model = 1
cn_dir = dir2_frag
suffix = paste0("_model", model)
plot_cell_with_maxcn(cn_dir, file.path(odir, subdir), suffix, chrsel, 3)


############ find Seismic Amplification ############ 
fban = "/Users/bl0033/Gdrive/git/CIN_SV/data/cytoBand_hg38.txt"
MV_BANDS = read_tsv(fban, col_names = F, show_col_types = F) 
names(MV_BANDS) = c("chrom", "chromStart", "chromEnd", "name", "gieStain")
chrBands = makeGRangesFromDataFrame(MV_BANDS, keep.extra.columns=TRUE)

#odir = "/Users/bl0033/Gdrive/git/cnv_analysis/cin_sv/bfb/result/local_frag/pcorr0/seismic"
odir = "/Users/bl0033/Gdrive/git/cnv_analysis/cin_sv/bfb/result/local_frag/pcorr0.5/seismic"
bdir = "/Users/bl0033/data/SV/test4paper/bfb"
chrsel = "chr2"

plot = F
ssdir = "r5"
#ssdir = "r2"

subdir1 = paste0("nDSB2_Un0.5_frag50_break0_model", 0, "_ptype1_pcorr0.5")
dir1 = file.path(bdir, subdir1, ssdir)
# subdir1 = paste0("nDSB2_Un0.5_frag50_break0_model", 0, "_ptype1_pcorr0")
# dir1 = file.path(bdir, subdir1, "r1")
cfiles1 = list.files(dir1, pattern = "^CNData_c.*tsv")
length(cfiles1)
# must ensure the out dir exists, or waste a long run
subodir = file.path(odir, "neutral")
dir.create(subodir, showWarnings = F)
fout1 = file.path(subodir, paste0("seismic_", ssdir, "_", subdir1, ".rds"))
res_list1 = get_seis_list(fout1, dir1, cfiles1)
cell_seis_neu = get_cell_with_seismic(cfiles1, res_list1, dir1, chrsel, plot)
#sa = get_cell_with_seismic(cfiles1, res_list1, dir1, chrsel, plot, "183")


subdir2 = paste0("nDSB2_Un0.5_frag50_break0_model", 1, "_ptype1_pcorr0.5")
dir2 = file.path(bdir, subdir2, ssdir)
# subdir2 = paste0("nDSB2_Un0.5_frag50_break0_model", 1, "_ptype1_pcorr0")
# dir2 = file.path(bdir, subdir2, "r1")
cfiles2 = list.files(dir2, pattern = "^CNData_c.*tsv")
length(cfiles2)
subodir = file.path(odir, "selection")
dir.create(subodir, showWarnings = F)
fout2 = file.path(subodir, paste0("seismic_", ssdir, "_", subdir2, ".rds"))
res_list2 = get_seis_list(fout2, dir2, cfiles2)
cell_seis_sel = get_cell_with_seismic(cfiles2, res_list2, dir2, chrsel, plot)


############ check cancer genes overlapping the altered region  ############ 
fcmc_sel = "/Users/bl0033/data/cmc_export_sel.tsv" 
if(!file.exists(fcmc_sel)){
  fname = "/Users/bl0033/data/cmc_export.tsv"  
  ccg = read_tsv(fname)
  csel = ccg[c(1:4, 21, 24)]
  names(ccg[c(1:4, 21, 24)])
  # 21 Mutation genome position GRCh38
  names(csel) = c("gname", "accession", "onc", "cgc", "pos", "disease")
  unique(csel$cgc)
  csel %>% separate(pos, c("chr", "coord"), sep=":") %>% separate(coord, c("start", "end"), sep="-") -> csep
  csep$start = as.numeric(csep$start)
  csep$end = as.numeric(csep$end)
  write_tsv(csep, fcmc_sel)
}else{
  csep = read_tsv(fcmc_sel)
}


########## for simple break 
# the breakpoints are unique for each cell
d_segT1 = get_cn_per_run_rck(dir1, type = "total")  # neutral
d_segT2 = get_cn_per_run_rck(dir2, type = "total")

chrsel = "7" # simple break
gene_overlap1 = get_gene_overlap(csep, d_segT1, chrsel)
gene_overlap2 = get_gene_overlap(csep, d_segT2, chrsel)


olp_gene1 = gene_overlap1 %>% select(sample, gname, cn, onc, disease) %>% group_by(sample) %>% distinct()
# write genes for reference 
gcount1 = olp_gene1 %>% filter(!is.na(onc)) %>% group_by(sample) %>% count(gname, onc, cn) 
gcount1 %>% group_by(sample) %>% tally()
onc1 = sort(unique(gcount1$gname))
gn1 = gcount1 %>% group_by(gname, onc, cn) %>% tally() 
fout = file.path(odir, "oncgene_count_neutral.tsv")
write_tsv(gn1, fout)

olp_gene2 = gene_overlap2 %>% select(sample, gname, cn, onc, disease) %>% group_by(sample) %>% distinct()
# write genes for reference 
gcount2 = olp_gene2 %>% filter(!is.na(onc)) %>% group_by(sample) %>% count(gname, onc, cn) 
gcount2 %>% group_by(sample) %>% tally()
onc2 = sort(unique(gcount2$gname))
gn2 = gcount2 %>% group_by(gname, onc, cn) %>% tally()
fout = file.path(odir, "oncgene_count_selection.tsv")
write_tsv(gn2, fout)

gn1_sel = gn1 %>% filter(n > 1) %>% filter((onc == "oncogene" & cn > 2) | (onc == "TSG" & cn < 2))
fout = file.path(odir, "oncgene_count_neutral_cna.tsv")
write_tsv(gn1_sel, fout)

gn2_sel = gn2 %>% filter(n > 1) %>% filter((onc == "oncogene" & cn > 2) | (onc == "TSG" & cn < 2))
fout = file.path(odir, "oncgene_count_selection_cna.tsv")
write_tsv(gn2_sel, fout)


# Plot onco gene copy numbers, not useful
ggplot(gn1_sel, aes(x = gname, y = n, fill = as.factor(cn))) + geom_bar(position = position_dodge(), stat = "identity") + theme_pubr() + scale_fill_npg()
ggplot(gn2_sel, aes(x = gname, y = n, fill = as.factor(cn))) + geom_bar(position = position_dodge(), stat = "identity") + theme_pubr() + scale_fill_npg()


########## for local fragmentation
# the breakpoints are unique for each cell
d_segT1 = get_cn_per_run_rck(dir1_frag, type = "total")  # neutral
d_segT2 = get_cn_per_run_rck(dir2_frag, type = "total")

chrsel = "2" # simple break
gene_overlap1 = get_gene_overlap(csep, d_segT1, chrsel)
gene_overlap2 = get_gene_overlap(csep, d_segT2, chrsel)

olp_gene1 = gene_overlap1 %>% select(sample, gname, cn, onc, disease) %>% group_by(sample) %>% distinct()
# write genes for reference 
gcount1 = olp_gene1 %>% filter(!is.na(onc)) %>% group_by(sample) %>% count(gname, onc, cn) 
gcount1 %>% group_by(sample) %>% tally()
onc1 = sort(unique(gcount1$gname))
gn1 = gcount1 %>% group_by(gname, onc, cn) %>% tally() 
fout = file.path(odir, "oncgene_count_neutral_frag.tsv")
write_tsv(gn1, fout)

olp_gene2 = gene_overlap2 %>% select(sample, gname, cn, onc, disease) %>% group_by(sample) %>% distinct()
# write genes for reference 
gcount2 = olp_gene2 %>% filter(!is.na(onc)) %>% group_by(sample) %>% count(gname, onc, cn) 
gcount2 %>% group_by(sample) %>% tally()
onc2 = sort(unique(gcount2$gname))
gn2 = gcount2 %>% group_by(gname, onc, cn) %>% tally()
fout = file.path(odir, "oncgene_count_selection_frag.tsv")
write_tsv(gn2, fout)

gn1_sel = gn1 %>% filter(n > 1) %>% filter((onc == "oncogene" & cn > 2) | (onc == "TSG" & cn < 2)) %>% arrange(desc(onc))
fout = file.path(odir, "oncgene_count_neutral_cna_frag.tsv")
write_tsv(gn1_sel, fout)

gn2_sel = gn2 %>% filter(n > 1) %>% filter((onc == "oncogene" & cn > 2) | (onc == "TSG" & cn < 2)) %>% arrange(desc(onc))
fout = file.path(odir, "oncgene_count_selection_cna_frag.tsv")
write_tsv(gn2_sel, fout)

# merge two tables
gn_sel = merge(gn1_sel, gn2_sel, by = c("onc", "gname", "cn"), all = T) %>% dplyr::rename(`#cells (neutral model)` = n.x, `#cells (selection model)` = n.y, `total copy number` = cn, `name` = gname, `type` = onc) %>% replace(is.na(.), 0)
fout = file.path(odir, "oncgene_count_frag.tsv")
write_tsv(gn_sel, fout)


############ check chromothripsis ############ 
chrs = c(2)

df_chromth1 = get_chromth_count_by_dir(dir1_frag, chrs)
write_chromth_summary(odir, df_chromth1, dsuffix = dsuffix1)

df_chromth2 = get_chromth_count_by_dir(dir2_frag, chrs)
write_chromth_summary(odir, df_chromth2, dsuffix = dsuffix2)

# select cells with chromothripsis
cell_chrmoth_neu = df_chromth1 %>% filter(nhconf > 0) %>% select(cell_ID) %>% unlist()
cell_chrmoth_sel = df_chromth2 %>% filter(nhconf > 0) %>% select(cell_ID) %>% unlist()


############ check ecDNA ############ 
bdir = odir
ncell = 100
dir = file.path(bdir, paste0("ncell", ncell), paste0("selection"))  # selection
dir2 = file.path(bdir, paste0("ncell", ncell), paste0("neutral"))  # neutral evolution

res0 = parse_ecdna_batch(bdir, ncell, mcn, maxn, mins, maxs, paste0(ncell, " cells"))
ec_sel0 = res0$ec_sel
ncount0 = res0$ncount
ec_all0 = res0$ec_all

necdna_all0 = ec_all0 %>% group_by(cell, type) %>% tally() %>% mutate(tcell = ncell)
mcn_all0 = ec_all0 %>% group_by(cell, type) %>% summarise(mcn = mean(n), maxcn = max(n)) %>% mutate(tcell = ncell)

fsstat0 = file.path(dir2_frag, "sumStats_sim.tsv")
ss_sel0 = read_tsv(fsstat0) %>% mutate(ncell = ncell) %>% mutate(type = "selection")
fsstat0 = file.path(dir1_frag, "sumStats_sim.tsv")
ss_neu0 = read_tsv(fsstat0) %>% mutate(ncell = ncell) %>% mutate(type = "neutral evolution")

# size, #ecdna, selection
psz = ggplot(ec_all0, aes(x = type, y = log10(size), fill = type)) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.1, fill="white") + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend(""))  + theme(legend.position = "top") + ylab("log10(ecDNA size)") + ggtitle("")  + xlab("") + theme(axis.text.x = element_blank()) + geom_hline(yintercept = 5, linetype = "dashed", color = "darkgrey")

pd = ggplot(necdna_all0, aes(x = type, y = n, fill = type)) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.1, fill="white") + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend("")) + ylab("#different ecDNAs per cell") + theme_pubr(base_size = 7) + theme(legend.position = "none") + xlab("") + theme(axis.text.x = element_blank()) + geom_hline(yintercept = 15, linetype = "dashed", color = "darkgrey")

pmc = ggplot(mcn_all0, aes(x = type, y = mcn, fill = type)) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.1, fill="white") + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend(""))  + theme(legend.position = "top") + ylab("mean ecDNA copy number per cell") + xlab("") + theme(axis.text.x = element_blank()) + geom_hline(yintercept = 1.8, linetype = "dashed", color = "darkgrey") 

pmmc = ggplot(mcn_all0, aes(x = type, y = maxcn, fill = type)) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.1, fill="white") + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend(""))  + theme(legend.position = "top") + ylab("maximum ecDNA copy number per cell") + xlab("") + theme(axis.text.x = element_blank()) + geom_hline(yintercept = 4, linetype = "dashed", color = "darkgrey")
  
ps = ggplot(ec_sel0, aes(x = type, y = selection_coef, fill = type)) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.1, fill="white") + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend("")) + ylab(lbl_scoef) + theme_pubr(base_size = 7) + theme(legend.position = "none") + xlab("") + theme(axis.text.x = element_blank()) + geom_hline(yintercept = 0.003, linetype = "dashed", color = "darkgrey")


ggarrange(pd, psz, pmmc, pmc, ps, nrow =1, common.legend = T)
fout = file.path(odir, "ecdna_n100.pdf")
ggsave(fout, width = 7, height = 3)



