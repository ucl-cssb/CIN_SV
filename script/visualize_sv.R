# Visualize the results of SVGen with circlize
library(circlize)
library(tidyverse)


get_SV <- function(fname){
  svs = read_tsv(fname)
  svs %>% rowwise() %>%  mutate(SVTYPE= strsplit(extra, "=")[[1]][2]) %>% as.data.frame() -> svs

  # sv.first and sv.second are dataframes with columns: "chr","start","end"
  sv.first <- svs %>% select(chr=chr1, start=coord1) %>% as.data.frame()
  sv.first$chr = paste0("chr", sv.first$chr)
  sv.first$end = sv.first$start

  sv.second <- svs %>% select(chr=chr2, start=coord2) %>% as.data.frame()
  sv.second$chr = paste0("chr", sv.second$chr)
  sv.second$end = sv.second$start

  return(list(svs = svs, sv.first=sv.first, sv.second=sv.second))
}


get_CN <- function(fname, cutoff = 5*10e4){
  cns = read_tsv(fname)

  cns %>% rowwise() %>% mutate(cnA = as.numeric(strsplit(strsplit(strsplit(extra, "=")[[1]][2],":")[[1]][3],",")[[1]][1]), cnB = as.numeric(strsplit(strsplit(strsplit(extra, "=")[[1]][2],":")[[1]][4],"}")[[1]][1])) -> cns

  cns$chr = paste0("chr", cns$chr)

  cns %>% select(chr, start, end, value1 = cnA, value2 = cnB) -> cnbed
  # unique(cns$cnA)
  # unique(cns$cnB)
  # summary(cnbed$end-cnbed$start)

  cnbed -> cnsel # %>% filter(value1 != 1 | value2 != 1)  %>% filter(end-start > cutoff)
  # str(cnsel)
  #
  # unique(cnsel$value1)
  # unique(cnsel$value2)

  return(cnsel)
}

fname = fsv
plot_SV <- function(fname){
  res = get_SV(fname)
  svs = res$svs
  sv.first = res$sv.first
  sv.second = res$sv.second

  # pa <- c("#377eb8","#d95f02","#33a02c")
  # names(pa) <- c("DEL","INS","INV")
  pa <- c("#1F78B4","#FF7F00","#33A02C","#E31A1C")
  names(pa) <- c("DEL","DUP","H2HINV","T2TINV")
  col <- pa[svs$SVTYPE]


    outfile <- str_replace(fname, ".tsv", ".png")
    name <- ""
    png(outfile, height=600, width=600)

    print(sv.first)
    circos.initializeWithIdeogram(species="hg38")
    circos.genomicLink(sv.first, sv.second, col=col)
    circos.clear()

    title(name)
    dev.off()
}


plot_CN_SV <- function(fsv, fcnv){

  # get CN data
  cnsel = get_CN(fcnv)

  # CN colors
  cols = c("#6283A9","#f0f0f0", "#B9574E", "#3b0107")
  col_fun = colorRamp2(breaks = seq(0:3)-1, colors = cols)

  # get SV data
  res = get_SV(fsv)
  svs = res$svs
  sv.first = res$sv.first
  sv.second = res$sv.second

  # SV colors
  # pa <- c("#377eb8","#d95f02","#33a02c")
  # names(pa) <- c("DEL","INS","INV")
  pa <- c("#1F78B4","#FF7F00","#33A02C","#E31A1C")
  names(pa) <- c("DEL","DUP","H2HINV","T2TINV")
  col_fun_sv = colorRamp2(c("DEL","DUP","H2HINV","T2TINV"), pa)
  col <- pa[svs$SVTYPE]


  # lgd_heatmap = Legend(at = seq(0:3)-1, col_fun = col_fun, legend_gp = gpar(col = 0:3), title_position = "topleft", title = "Copy number")
  # lgd_links = Legend(at = c("DEL","DUP","H2HINV","T2TINV"), type = "lines", itle_position = "topleft", title = "Links")
  # lgd_list_vertical = packLegend(lgd_heatmap, lgd_links)
  # lgd_list_vertical

  fname = fsv
  outfile <- str_replace(fname, ".tsv", ".png")
  name <- ""
  png(outfile, height=600, width=600)

  circos.initializeWithIdeogram(species="hg38", chromosome.index = paste0("chr", seq(1:22)))

  cnsel_cycle = cnsel

  circos.genomicHeatmap(cnsel_cycle, col = col_fun, side = "inside", border = "white")
  circos.genomicLink(sv.first, sv.second, col=col)
  circos.clear()

  #draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

  title(name)
  dev.off()
}

dir = "/Users/ucbtlux/Gdrive/git/sv-model/test/"

for(cellID in c(2:3)){
  # cellID = 2
  fsv = file.path(dir, paste0("sv_data_c", cellID, ".tsv"))
  #plot_SV(fsv)
  fcnv = file.path(dir, paste0("cn_data_c", cellID, ".tsv"))
  plot_CN_SV(fsv, fcnv)
}
