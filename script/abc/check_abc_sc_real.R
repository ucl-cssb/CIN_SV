#!/usr/bin/env Rscript

# This script is used to parse the results of ABC on real data (Fig 6)


# change path accordingly
source("parse_sim_util.R")
source("abc/check_abc_util.R")


################### read meta information ###################
# save Table S6 to TSV file 
fmeta = "data_svmodel/Fig6/patient_nclone_nbp.tsv"
dmeta = read_tsv(fmeta, show_col_types = F)
# nWGD_clone: #clones with WGD; n: #clones
dmeta$fracWGD = dmeta$nWGD_clone/dmeta$n


################# basic setting for real data ################# 
max_ndsb = 50
max_unrepair = 1
max_wgd = 1
base_size = 7

rdir = "data_svmodel/Fig6/smcres"  # need to run ABC SMC to get the results for each dataset
odir = "data_svmodel/Fig6"
ddir = "data_svmodel/Fig6/stat"

pattern = "smc_nparam.*_epsilon0.15.*_maxWGD0.5_.*[a|b|0-9]$"
pdx_types = c("FBI")

################ read results in a batch #################
to_plot_stat = T 
to_plot_ppc = F
res = get_smc_all_real(rdir, max_wgd, max_unrepair, to_plot)

resall = res$resall
meanvals = res$meanvals

ressel = merge(resall, dmeta, by = c("sample"))
meansel0 = merge(meanvals, dmeta, by = c("sample"))


####### get  results #######
# take long time to process, easy to stuck, save for future processing
fpps = file.path(odir, "pps_real.rds")
if(file.exists(fpps)){
  print("reading stored results")
  res_pps = readRDS(fpps)
}else{
  res_pps = get_sstat_all(rdir, ddir, to_plot_stat, to_plot_ppc)
  saveRDS(res_pps, fpps)
  res_pps1 = readRDS(fpps)
  setequal(res_pps, res_pps1)
}

sstat_all = res_pps$sstat_all %>% rename(sample = dataset)
cstat_all = res_pps$cstat_all %>% rename(sample = dataset)

sstat = merge(sstat_all, dmeta, by = c("sample"))
cstat = merge(cstat_all, dmeta, by = c("sample"))
 


####### compute mean ####### 
sstat_avg = sstat %>% group_by(sample) %>% mutate(mean_avg_fusion = mean(avg_fusion), sd_fusion = sd(avg_fusion)) %>% select(sample, mean_avg_fusion, sd_fusion) %>% unique()
cstat_avg = cstat %>% group_by(sample) %>% mutate(mean_nECDNA = mean(nECDNA), sd_nECDNA = sd(nECDNA)) %>% select(sample, mean_nECDNA, sd_nECDNA) %>% unique()

mean_all0 = merge(meansel0, sstat_avg, by = c("sample"))
mean_all = merge(mean_all0, cstat_avg, by = c("sample"))
mean_sel = mean_all %>% select(sample, n, nSV_subclonal, ntc, rbp, dsb_rate, frac_unrepaired, wgd, wsd_dsb, wsd_unrepair, wsd_wgd, mean_avg_fusion, mean_nECDNA, sd_fusion, sd_nECDNA, type, signature_type) %>% mutate(wgd = ifelse(wgd < 0, 0, wgd))

# used for correlation plots
mean4corr = mean_sel %>% select(`DSB rate` = dsb_rate, `%unrepaired DSBs` = frac_unrepaired, `probability of WGD` = wgd, `#fusions` = mean_avg_fusion, `#ecDNAs` = mean_nECDNA)


################### plot results and check correlation by parameters to group datasets ###################
ordered = F
nbreak = 5
fsize = 6
w = 7
h = 6
# for fig 6 in paper
pres = plot_smc_all_final(meansel0, ressel, sstat, cstat, odir, data_colors, nbreak, fsize, w, h)

############# check noval correlations
corr_type = "spearman" # used in main text

pf = ggplot(mean_sel, aes(x = wgd, y = mean_avg_fusion)) + geom_point() + theme_pubr(base_size = fsize) + geom_smooth(method = 'lm') + xlab(lbl_wgd) + ylab(lbl_fusion) + stat_cor(method = corr_type, size = fsize/2)
# + geom_errorbar(aes(ymin=mean_avg_fusion-sd_fusion, ymax=mean_avg_fusion+sd_fusion), width=.2, position=position_dodge(0.05))
pe = ggplot(mean_sel, aes(x = wgd, y = mean_nECDNA)) + geom_point() + theme_pubr(base_size = fsize) + geom_smooth(method = 'lm') + xlab(lbl_wgd) + ylab(lbl_ecdna) + stat_cor(method = corr_type, size = fsize/2)
pfe = ggplot(mean_sel, aes(x = mean_avg_fusion, y = mean_nECDNA)) + geom_point() + theme_pubr(base_size = fsize) + geom_smooth(method = 'lm') + xlab(lbl_fusion) + ylab(lbl_ecdna) + stat_cor(method = corr_type, size = fsize/2)

pd = ggplot(mean_sel, aes(x = wgd, y = dsb_rate)) + geom_point() + theme_pubr(base_size = fsize) + geom_smooth(method = 'lm') + xlab(lbl_wgd) + ylab(lbl_dsb) + stat_cor(method = corr_type, size = fsize/2)

pcorr = ggarrange(pd, pf, pe, pfe, nrow = 1)

ggarrange(pres, pcorr, nrow = 2, heights = c(4.5, 1))
fname = paste0("sum_real_plot_final.pdf")
fout = file.path(odir, fname)
w = 7
h = 7.8
ggsave(fout, width = w, height = h)


chart.Correlation(mean4corr, histogram=TRUE, pch=15, method = corr_type)
fout = file.path(odir, "corr_param_chart.pdf")
pdf(fout, width = 7, height = 7)
plot_corr_chart(mean4corr, histogram=TRUE, nbreak = 8, pch=15, method = corr_type)
dev.off()



############  compare htert and fbi (for extended figure) ###########
mean4cmp = mean_sel %>% select(type, `DSB rate` = dsb_rate, `%unrepaired DSBs` = frac_unrepaired, `probability of WGD` = wgd, `#fusions` = mean_avg_fusion, `#ecDNAs` = mean_nECDNA)
p1 = ggplot(mean4cmp, aes(x = type, y = `DSB rate`)) + geom_boxplot() + theme_pubr(base_size = base_size) + xlab("") + ylab(lbl_dsb)
p2 = ggplot(mean4cmp, aes(x = type, y = `%unrepaired DSBs`)) + geom_boxplot() + theme_pubr(base_size = base_size) + xlab("") + ylab(lbl_unrepair)
p3 = ggplot(mean4cmp, aes(x = type, y = `probability of WGD`)) + geom_boxplot() + theme_pubr(base_size = base_size) + xlab("") + ylab(lbl_wgd)
p4 = ggplot(mean4cmp, aes(x = type, y = `#fusions`)) + geom_boxplot() + theme_pubr(base_size = base_size) + xlab("") + ylab(lbl_fusion)
p5 = ggplot(mean4cmp, aes(x = type, y = `#ecDNAs`)) + geom_boxplot() + theme_pubr(base_size = base_size) + xlab("") + ylab(lbl_ecdna)

ggarrange(p1, p2, p3, p4, p5, nrow = 1, ncol = 5)

fout = file.path(odir, "box_htert_pdx.pdf")
ggsave(fout, width = 7, height = 3)



################### check inferred DSB rate and WGD probability ###################
# compare the difference of weighted mean and fracWGD
rdsb = ressel %>% group_by(sample) %>% summarise(mdsb = weighted.mean(dsb_rate,  weight), wsd_dsb = sqrt(wtd.var(dsb_rate, weight, normwt = T))) 
rdsb_all = merge(rdsb, dmeta, by = c("sample")) %>% mutate(diff = abs(mdsb - rbp)) %>% mutate(diff_frac = diff / mdsb) %>% mutate(xcolor = ifelse(fracWGD<=0, "darkgrey", "black"))


rdsb_sel = rdsb_all %>% gather(mdsb, rbp, key = "type", value = "val") %>% mutate(wsd_dsb = if_else(type == "rbp", 0, wsd_dsb))
rdsb_sel$type = factor(rdsb_sel$type, levels = c("rbp", "mdsb"))
length(unique(rdsb_sel$sample))
# arrange by mwgd 
sorted_spl = rdsb_all %>% arrange(rbp) %>% select(sample) %>% unlist()
rdsb_sel$sample = factor(rdsb_sel$sample, levels = sorted_spl)

lbl_r1 = "empirical DSB rate per cycle"
lbl_r2 = "estimated DSB rate per cycle"
pcd = ggplot(rdsb_sel, aes(x = sample, y = val, fill = type)) + geom_bar(stat="identity", position=position_dodge()) + theme_pubr(base_size = base_size) + theme(axis.text.x = element_text(angle = 90), legend.title = element_blank()) + xlab("") + ylab("value") + scale_fill_manual(breaks =  c("rbp", "mdsb"), values = cmp_colors, labels = c(lbl_r1, lbl_r2)) + theme(axis.text.x = element_text(angle = 45, hjust=1))  + theme(axis.text.x = element_text(colour=rdsb_sel$xcolor))  + geom_errorbar(aes(ymin=val-wsd_dsb, ymax=val+wsd_dsb), width=.2, position=position_dodge(.9))


ggarrange(pcd, pcw, nrow = 2, ncol = 1)

fout = file.path(odir, "cmp_dsb_wgd.pdf")
ggsave(fout, width = 7, height = 5)



