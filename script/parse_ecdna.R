#!/usr/bin/env Rscript

# This script is used to parse the simulated data with ecDNAs (Fig 4)

# change path accordingly
source("parse_sim_util.R")


############## check long term ecDNA dynamics ############## 
############# 5000 cells
# run first to get axis limits
ncell = 5000
dir = file.path(bdir, paste0("ncell", ncell), paste0("selection"))  # selection
dir2 = file.path(bdir, paste0("ncell", ncell), paste0("neutral"))  # neutral evolution

res5 = parse_ecdna_batch(bdir, ncell, mcn, maxn, mins, maxs, paste0(ncell, " cells"))
ec_all5 = res5$ec_all
ncount5 = res5$ncount
ec_sel5 = res5$ec_sel

necdna_all5 = ec_all5 %>% group_by(cell, type) %>% tally() %>% mutate(tcell = ncell)
mcn_all5 = ec_all5 %>% group_by(cell, type) %>% summarise(mcn = mean(n), maxcn = max(n)) %>% mutate(tcell = ncell)

mcn = max(ec_all5$n)
maxn = max(necdna_all5$n)
mins = min(ec_sel5$selection_coef)
maxs = max(ec_sel5$selection_coef)

fsstat1 = file.path(dir, "sumStats_sim.tsv")
ss_sel5 = read_tsv(fsstat1) %>% mutate(ncell = ncell) %>% mutate(type = "selection")
fsstat1 = file.path(dir2, "sumStats_sim.tsv")
ss_neu5 = read_tsv(fsstat1) %>% mutate(ncell = ncell) %>% mutate(type = "neutral evolution")


############## 3000 cells
ncell = 3000
dir = file.path(bdir, paste0("ncell", ncell), paste0("selection"))  # selection
dir2 = file.path(bdir, paste0("ncell", ncell), paste0("neutral"))  # neutral evolution

res3 = parse_ecdna_batch(bdir, ncell, mcn, maxn, mins, maxs, paste0(ncell, " cells"))
ec_sel3 = res3$ec_sel
ncount3 = res3$ncount
ec_all3 = res3$ec_all

necdna_all3 = ec_all3 %>% group_by(cell, type) %>% tally() %>% mutate(tcell = ncell)
mcn_all3 = ec_all3 %>% group_by(cell, type) %>% summarise(mcn = mean(n), maxcn = max(n)) %>% mutate(tcell = ncell)

fsstat1 = file.path(dir, "sumStats_sim.tsv")
ss_sel3 = read_tsv(fsstat1) %>% mutate(ncell = ncell) %>% mutate(type = "selection")
fsstat1 = file.path(dir2, "sumStats_sim.tsv")
ss_neu3 = read_tsv(fsstat1) %>% mutate(ncell = ncell) %>% mutate(type = "neutral evolution")


######## 1000 cells 
ncell = 1000
dir = file.path(bdir, paste0("ncell", ncell), paste0("selection"))  # selection
dir2 = file.path(bdir, paste0("ncell", ncell), paste0("neutral"))  # neutral evolution

res1 = parse_ecdna_batch(bdir, ncell, mcn, maxn, mins, maxs, paste0(ncell, " cells"))
ec_sel1 = res1$ec_sel
ncount1 = res1$ncount
ec_all1 = res1$ec_all

necdna_all1 = ec_all1 %>% group_by(cell, type) %>% tally() %>% mutate(tcell = ncell)
mcn_all1 = ec_all1 %>% group_by(cell, type) %>% summarise(mcn = mean(n), maxcn = max(n)) %>% mutate(tcell = ncell)

fsstat1 = file.path(dir, "sumStats_sim.tsv")
ss_sel1 = read_tsv(fsstat1) %>% mutate(ncell = ncell) %>% mutate(type = "selection")
fsstat1 = file.path(dir2, "sumStats_sim.tsv")
ss_neu1 = read_tsv(fsstat1) %>% mutate(ncell = ncell) %>% mutate(type = "neutral evolution")



########### check ecDNA types  ###########
necdna_cmb = rbind(necdna_all1, necdna_all3, necdna_all5)

pe1 = ggplot(necdna_cmb, aes(x = as.factor(tcell), y = log10(n), fill = type)) +
  geom_violin(scale = "width", position=position_dodge(width = 0.9)) +
  geom_boxplot(aes(group=interaction(type, as.factor(tcell))), width=0.1, fill="white", position=position_dodge(width = 0.9)) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend("")) + ylab("log10(#different ecDNAs per cell)") + theme_pubr(base_size = 7) + ggtitle("") + xlab(lbl_ncell) + geom_hline(yintercept = 1, linetype = "dashed", color = "darkgrey") + geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgrey") + theme(legend.position = "none")

necdna_mean = necdna_cmb %>% group_by(tcell, type) %>% summarise(men = mean(n))
# necdna_mean = necdna_cmb %>% group_by(tcell, type) %>% summarise(men = mean(log10(n)))
pe2 = ggplot(necdna_mean, aes(y = men, x = as.factor(tcell), color = type, group = type)) + geom_point() + geom_line() + theme_pubr(base_size = 7) + scale_color_manual(breaks = model_breaks, values = model_colours)  + guides(color = guide_legend("")) + xlab(lbl_ncell) + ylab("mean #different ecDNAs per cell") + theme(legend.position = "none")

pe = (pe1 | pe2) + plot_layout(widths = c(3,1))


################## check ecDNA size  ###########
ec_cmb = rbind(ec_all1, ec_all3, ec_all5)

ec_cmb %>% select(tcell, cell, type) %>% unique() %>% group_by(tcell, type) %>% tally()

pz1 = ggplot(ec_cmb, aes(x = as.factor(tcell), y = log10(size), fill = type))  +
  geom_violin(scale = "width", position=position_dodge(width = 0.9)) +
  geom_boxplot(aes(group=interaction(type, as.factor(tcell))), width=0.1, fill="white", position=position_dodge(width = 0.9))  + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend(""))  + theme(legend.position = "none") + ylab("log10(ecDNA size)") + ggtitle("")  + xlab(lbl_ncell) + geom_hline(yintercept = 4.5, linetype = "dashed", color = "darkgrey")

size_mean = ec_cmb %>% group_by(tcell, type) %>% summarise(msize = mean(size))
# lbl_esize = "mean log10(ecDNA size)"
lbl_esize = "mean ecDNA size"
pz2 = ggplot(size_mean, aes(y = (msize), x = as.factor(tcell), color = type, group = type)) + geom_point() + geom_line() + theme_pubr(base_size = 7) + scale_color_manual(breaks = model_breaks, values = model_colours)  + guides(color = guide_legend("")) + xlab(lbl_ncell) + ylab(lbl_esize) + theme(legend.position = "none")

psize = (pz1 | pz2) + plot_layout(widths = c(3,1))


########## check max CN ##########
mcn_cmb = rbind(mcn_all1, mcn_all3, mcn_all5)

pm1 = ggplot(mcn_cmb, aes(x = as.factor(tcell), y = log10(maxcn), fill = type)) +
  geom_violin(scale = "width", position=position_dodge(width = 0.9)) +
  geom_boxplot(aes(group=interaction(type, as.factor(tcell))), width=0.1, fill="white", position=position_dodge(width = 0.9)) + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend(""))  + theme(legend.position = "none") + ylab("log10(maximum ecDNA copy number per cell)")  + xlab(lbl_ncell) + geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkgrey")

mcn_mean = mcn_cmb %>% group_by(tcell, type) %>% summarise(mmcn = mean(maxcn))
pm2 = ggplot(mcn_mean, aes(y = mmcn, x = as.factor(tcell), color = type, group = type)) + geom_point() + geom_line() + theme_pubr(base_size = 7) + scale_color_manual(breaks = model_breaks, values = model_colours)  + guides(color = guide_legend("")) + xlab(lbl_ncell) + ylab("mean of maximum ecDNA copy number per cell") + theme(legend.position = "none")

pmax = (pm1 | pm2) + plot_layout(widths = c(3,1))

################## check mean ecDNA copy number per cell ###########
pm1 = ggplot(mcn_cmb, aes(x = as.factor(tcell), y = log10(mcn), fill = type)) +
  geom_violin(scale = "width", position=position_dodge(width = 0.9)) +
  geom_boxplot(aes(group=interaction(type, as.factor(tcell))), width=0.1, fill="white", position=position_dodge(width = 0.9)) + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend(""))  + theme(legend.position = "none") + ylab("log10(mean ecDNA copy number per cell)")  + xlab(lbl_ncell) + geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgrey")

mcn_mean = mcn_cmb %>% group_by(tcell, type) %>% summarise(mmcn = mean(mcn))
pm2 = ggplot(mcn_mean, aes(y = mmcn, x = as.factor(tcell), color = type, group = type)) + geom_point() + geom_line() + theme_pubr(base_size = 7) + scale_color_manual(breaks = model_breaks, values = model_colours)  + guides(color = guide_legend("")) + xlab(lbl_ncell) + ylab("mean of mean ecDNA copy number per cell") + theme(legend.position = "none")

pmean = (pm1 | pm2) + plot_layout(widths = c(3,1))

############# check simulation summary stats ############
#ss_all = rbind(ss_sel0, ss_sel1, ss_sel3, ss_sel5, ss_neu0, ss_neu1, ss_neu3, ss_neu5)
ss_all = rbind(ss_sel1, ss_sel3, ss_sel5, ss_neu1, ss_neu3, ss_neu5)
ss_all = ss_all %>% mutate(type = if_else(type == "neutral", "neutral evolution", type))
pf2 = ggplot(ss_all, aes(y = nTelofusion, x = as.factor(ncell), fill = type, group = type)) + geom_bar(stat = "identity", position = position_dodge()) + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks2, values = model_colours2)  + guides(fill = guide_legend("")) + xlab(lbl_ncell) + ylab(lbl_fusion) + theme(legend.position = "top")

#ncount_all = rbind(ncount0, ncount1, ncount3, ncount5) %>% mutate(prop2 = n.y / ncell)
ncount_all = rbind(ncount1, ncount3, ncount5) %>% mutate(prop2 = n.y / ncell)
ncell_sel = ncount_all %>% filter(type != "neutral evolution")

pf0 = ggplot(ncell_sel, aes(y = n.y, x = as.factor(ncell), fill = type)) + geom_bar(stat = "identity", position = position_dodge()) + theme_pubr(base_size = 7) + scale_fill_manual(breaks = model_breaks[2:3], values = model_colours[2:3])  + guides(fill = guide_legend("")) + xlab(lbl_ncell) + ylab("#cells under selection") + theme(legend.position = "top")

pf1 = ggplot(ncount_all, aes(y = prop, x = as.factor(ncell), color = type, group = type)) + geom_point() + geom_line()  + theme_pubr(base_size = 7) + scale_color_manual(breaks = model_breaks, values = model_colours)  + guides(color = guide_legend("")) + xlab(lbl_ncell) + ylab("fraction of cells with ecDNAs") + theme(legend.position = "top")

pstat = (pf2 | pf0 | pf1) + plot_layout(widths = c(1, 1, 1)) 


############# check CNAs of a specific ecDNA ############# 
top_node = ec_cmb %>% select(nodes, type, cell) %>% unique() %>% group_by(nodes) %>% tally() %>% arrange(desc(n)) 

# select the most common one
top_n1 = top_node %>% head(1) %>% select(nodes) %>% unlist()
ec_sel = ec_cmb %>% filter(nodes == top_n1)
# ec_sel %>% group_by(nodes) %>% tally()
ec_seln = ec_sel %>% group_by(type, cell, n, tcell) %>% tally()
ec_selnn = ec_seln %>% group_by(type, n, tcell) %>% summarise(ncopy = n()) %>% mutate(freq = ncopy / tcell)

pcndc = ggplot(ec_sel, aes(x = (n), color = as.factor(tcell))) + geom_density() + theme(legend.position = "none") + xlab("copy number of most common ecDNA per cell") + facet_grid(. ~ type) + theme_pubr(base_size = 7) + scale_color_manual(values = cell_colors) + labs(color = lbl_ncell) 


############  combine plots for main figure ############

ggarrange(pstat, pcndc, pe, psize, pmax, pmean, nrow = 3, ncol = 2) 
fout = file.path(odir, "ecdna_all.pdf")
ggsave(fout, width = 8, height = 7)


########### check selection coef ############
ec_sel_cmb = rbind(ec_sel1, ec_sel3, ec_sel5)
ps1 = ggplot(ec_sel_cmb, aes(x = as.factor(tcell), y = selection_coef, fill = type)) +
  geom_violin(scale = "width", position=position_dodge(width = 0.9)) +
  geom_boxplot(aes(group=interaction(type, as.factor(tcell))), width=0.1, fill="white", position=position_dodge(width = 0.9)) + scale_fill_manual(breaks = model_breaks, values = model_colours) + guides(fill = guide_legend("")) + ylab(lbl_scoef) + theme_pubr(base_size = 7) + ggtitle("") + xlab(lbl_ncell)  + geom_hline(yintercept = 0.005, linetype = "dashed", color = "darkgrey") + geom_hline(yintercept = -0.0035, linetype = "dashed", color = "darkgrey")

sel_mean = ec_sel_cmb %>% group_by(tcell, type) %>% summarise(msel = mean(selection_coef))
ps2 = ggplot(sel_mean, aes(y = msel, x = as.factor(tcell), color = type, group = type)) + geom_point() + geom_line() + theme_pubr(base_size = 7) + scale_color_manual(breaks = model_breaks, values = model_colours)  + guides(color = guide_legend("")) + xlab(lbl_ncell) + ylab("mean selection coefficient per cell") + theme(legend.position = "none")

(ps1 | ps2) + plot_layout(widths = c(3,1))
fout = file.path(odir, "sel_all.pdf")
ggsave(fout, width = 7, height = 4)





