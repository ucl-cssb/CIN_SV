#!/usr/bin/env Rscript

# This script is used to detect and plot chromplexy events simulated in Fig 2

####################  check Jabba output on simulated data ##################### 

bdir = "data_svmodel/Fig2/chromplexy/r1/"
# bdir = "data_svmodel/Fig2/chromplexy/r2/"

# subdirs = list.files(bdir, pattern = "^c[0-9]$")  # for c0 to c9
subdirs = list.files(bdir, pattern = "^c[0-9]*[0-9]$")
ndir = length(subdirs)
ndir 

for(i in 1:ndir){
  # i = 5
  d = subdirs[i]
  print(d)
  dir_rck = file.path(bdir, d)
  print(dir_rck)
  
  res_rck = try(gGnome::gG(rck = dir_rck))
  if(inherits(res_rck, "try-error")){
    #error handling code, maybe just skip this iteration using
    next
  }

  saveRDS(res_rck, file.path(dir_rck, 'rck.rds'))
  
  # may fail with error: Error in aggregate.formula(formula = subject.id ~ grl.id, data = m, FUN = function(x) numwin -  : 
#  argument 'x' is  missing -- it has been renamed from 'formula'
  evs = try(events(res_rck, verbose = TRUE))
  if(inherits(evs, "try-error")){
    #error handling code, maybe just skip this iteration using
    next
  }

  if(nrow(evs$meta$event) == 0){
    next    
  }
  
  ntype = (evs$meta$event[, table(type)])

  if(!"chromoplexy" %in% names(ntype)) next
  ncpx = ntype["chromoplexy"]
  ncpx
  
  ####### check chromoplexy
  if(ncpx > 0){
    chmplx = chromoplexy(res_rck)
    print(chmplx)
    cms = unique(chmplx$edges$dt$chromoplexy)
    cms_val = cms[!is.na(cms)]
    print(cms_val)
    
    for(i in 1:length(cms_val)){
      cmi = cms_val[i]
      fout = file.path(dir_rck, paste0("chromoplexy_e", cmi, ".pdf"))
      pdf(fout)
      ed = chmplx$edges[which(chromoplexy == cmi)]$shadow
      rl = 5e7  # plot error when too much compared with the range
      plot(chmplx$gt,  ed %>% streduce(rl))
      # chmplx$edges$dt$chromoplexy
      title(paste0("Chromoplexy ", cmi))
      dev.off()      
    }
  }
  
  res_rck$json(filename = file.path(dir_rck, 'test.json'))
}


########## plot lineage tree #############
tsize = 2.5
w = 7.5
h = 4
wsplit = c(1, 1.5)
odir = "/Users/bl0033/Gdrive/git/cnv_analysis/cin_sv/chromoplexy/result"
htitle = "total CN"

# highlight node with chromoplexy
bdir = "data_svmodel/Fig2/chromplexy/r1"
flin = file.path(bdir, "cell_lineage.tsv")
sel_nodes = c(4, 6, 7)
# plot_lineage_tree(flin, sel_nodes)
# plot lineage with CN heatmap
plot_lineage_with_cn(bdir, odir, "_n5", theme0, T, tsize, w, h, sel_nodes, wsplit, htitle)


bdir = "data_svmodel/Fig2/chromplexy/r2"
flin = file.path(bdir, "cell_lineage.tsv")
sel_nodes = c(7, 10, 12, 13, 16, 17)
# plot_lineage_tree(flin, sel_nodes)
plot_lineage_with_cn(bdir, odir, "_n10", theme0, T, tsize, w, h, sel_nodes, wsplit, htitle)





