#!/usr/bin/env Rscript

# This script is used to parse the results of ABC on simulated data (Fig 5)

# change path accordingly
source("parse_sim_util.R")
source("abc/check_abc_util.R")

# take a long time to parse all data, need parallel computing

###################  check SMC results on simulated data ###################  
ncell = 10
max_ndsb = 50
max_unrepair=1
max_wgd = 1
epsilon = 0.2
header_ppc1 = c("pga", "mean_div", "sd_div", "pgaA", "mean_divA", "sd_divA", "pgaB", "mean_divB", "sd_divB", paste0("bp", c(1:ncell)), "frac_wgd", "frac_del", "frac_dup", "frac_inv", "frac_tra")

odir = "data_svmodel/Fig5"
rdir = "data_svmodel/Fig5/smcres"  # need to run ABC SMC to get the results for each dataset

dirs = list.dirs(rsdir, recursive = F)
length(dirs)

###################  get inferred parameters ###################  
# save to reduce load time 
fsmc = file.path(odir, "smc_sim.rds")

if(file.exists(fmc)){
  smc_all = readRDS(fsmc)
}else{
  c1<-makeForkCluster(6, outfile = "debug.txt")
  CE<-clusterEvalQ(c1, .libPaths())
  # clusterEvalQ(c1, library(doParallel)) # there is no package called ‘doParallel’
  registerDoParallel(c1)
  print(paste0("Cores = ",detectCores()))
  
  # for each combination of real_rDSB, real_frac, real_wgd
  check_ppc = T
  plot = T
  smc_all = data.table()

  # take ~1h
  t1 = Sys.time()
  smc_all = foreach(i=1:length(dirs), .combine='rbind') %dopar% {
    # print(i)
    dir = dirs[i]
    # print(dir)
    res_smc = get_one_sim_all(dir, ncell, max_ndsb, max_unrepair, max_wgd, epsilon, header_ppc1, odir, check_ppc, plot, width = 14, height = 20)
    res_smc
  }
  t2 = Sys.time()
  print(t2 - t1)
  stopCluster(c1)
  saveRDS(smc_all, fsmc)
}

# 10 runs for each parameter setting 
meanvals = smc_all %>% group_by(real_rDSB, real_frac, real_wgd, real_fusion, real_ecdna, sample) %>% summarise(dsb_rate = weighted.mean(dsb_rate, weight), frac_unrepaired = weighted.mean(frac_unrepaired, weight), wgd = weighted.mean(wgd, weight), fusion = mean(avg_fusion), ecdna = mean(nECDNA)) %>% mutate(diff_dsb = dsb_rate - real_rDSB, diff_frac = frac_unrepaired - real_frac, diff_wgd = wgd - real_wgd, diff_fusion = fusion - real_fusion, diff_ecdna = ecdna - real_ecdna) 


############# check parameter correlation #############
corr_type = "spearman" # used in main text
fsize = 6

p1 = ggplot(meanvals, aes(x = wgd, y = dsb_rate)) + geom_point() + theme_pubr(base_size = fsize) + geom_smooth(method = 'lm') + xlab(lbl_wgd) + ylab(lbl_dsb) + stat_cor(method = corr_type, size = fsize/2)  + scale_color_manual(values = unrep_colors) + labs(color = lbl_unrepair)

p2 = ggplot(meanvals, aes(x = wgd, y = fusion)) + geom_point() + theme_pubr(base_size = fsize) + geom_smooth(method = 'lm') + xlab(lbl_wgd) + ylab(lbl_fusion) + stat_cor(method = corr_type, size = fsize/2)  + scale_color_manual(values = unrep_colors) + labs(color = lbl_unrepair)

p3 = ggplot(meanvals, aes(x = wgd, y = ecdna)) + geom_point() + theme_pubr(base_size = fsize) + geom_smooth(method = 'lm') + xlab(lbl_wgd) + ylab(lbl_ecdna) + stat_cor(method = corr_type, size = fsize/2)  + scale_color_manual(values = unrep_colors) + labs(color = lbl_unrepair)

p4 = ggplot(meanvals, aes(x = fusion, y = ecdna)) + geom_point() + theme_pubr(base_size = fsize) + geom_smooth(method = 'lm') + xlab(lbl_fusion) + ylab(lbl_ecdna) + stat_cor(method = corr_type, size = fsize/2)  + scale_color_manual(values = unrep_colors) + labs(color = lbl_unrepair)

ggarrange(p1, p2, p3, p4, nrow = 1)

fout = file.path(odir, "sim_abc_corr.pdf")
ggsave(fout, width = 7, height = 3) 


############## summary results  ##############
y1 = "((estimated - real) / real) (DSB rate per cycle)"
y2 = "((estimated - real) / real) (%unrepaired DSBs per cycle)"
y3 = "((estimated - real) / real) (probability of WGD per cell)"
y4 = "((estimated - real) / real) (mean #fusions per cycle)"
y5 = "((estimated - real) / real) (mean #ecDNAs per cell)"
# there may be no ecDNA in the data, / 0 -> INF
meanvals2 = meanvals %>% mutate(prop_diff_dsb = diff_dsb / real_rDSB, prop_diff_frac = diff_frac / real_frac, prop_diff_wgd = diff_wgd / real_wgd, prop_diff_fusion = diff_fusion / real_fusion, prop_diff_ecdna = ifelse(real_ecdna==0, diff_ecdna, diff_ecdna / real_ecdna))


pp = plot_diff_prop_by_wgd(meanvals2, y1, y2, y3, y4, y5, base_size)

fout = file.path(odir, "sim_abc_prop_diff_wgd.pdf")
ggsave(fout, pp, width = width, height = 6) 


######### plot results of individual parameter settings (Suppl Figures) ###########
ngrp = meanvals %>% group_by(real_rDSB, real_frac, real_wgd) %>% tally() %>% arrange(real_rDSB, real_frac, real_wgd, n) %>% nrow()
rDSBs = sort(unique(smc_all$real_rDSB))
rfracs = sort(unique(smc_all$real_frac))
rwgds = sort(unique(smc_all$real_wgd))

max_fusion = ceiling(max(smc_all$avg_fusion))
max_ecdna = ceiling(max(smc_all$nECDNA))


for(sel_rDSB in rDSBs){
  #sel_rDSB = rDSBs[1]
  plist = get_plots_params_by_dsb(sel_rDSB, rfracs, rwgds, smc_all, meanvals)

  # batch violin plot to see the trend
  ggarrange(plotlist = plist, ncol = 1, nrow = ngrp / 2)
  
  fout = file.path(odir, paste0("sim_abc_all_rdsb", sel_rDSB, ".pdf"))
  ggsave(fout, width = 12, height = 14) 
}


