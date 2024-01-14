#!/usr/bin/env Rscript

# Get reference centromere and telomere positions

hgv = "hg19"
# hgv = "hg38"

fname = paste0("../data/cytoBand_", hgv, ".txt")
fout = paste0("../data/", hgv, "_size.tsv")

# convert into "chr length boundary (centromere)"
# ref = read.table(fname, header = F)
# names(ref) <- c("chrom","start","end","cytoband", "stain")
# unique(ref$stain)


get_cyto_band <- function(fcyto){
  cyto=read.table(fcyto, header=FALSE)
  names(cyto) <- c("chrom","start","end","cytoband", "stain")
  cyto$chrom = sub("chr", "", cyto$chrom)
  chroms2incl <- c(seq(1:22), "X", "Y")
  cyto %>% filter(chrom %in% chroms2incl) -> cyto

  # Find the boundary between two arms
  arm_boudaries=data.frame()
  chr_ends=data.frame()
  for (i in 1:(nrow(cyto)-1)) {
    chr1 = cyto[i,]$chrom
    arm1 = substr(cyto[i,]$cytoband, 1, 1)
    chr2 = cyto[i+1,]$chrom
    arm2 = substr(cyto[i+1,]$cytoband, 1, 1)

    if(chr1 == chr2 && arm1 != arm2){
      item=data.frame(chrom=chr1, boundary = cyto[i+1,]$start)
      arm_boudaries = rbind(arm_boudaries, item)
      print(paste("chrom", chr1, "arm boundary", cyto[i+1,]$start))
    }
    if(chr1 != chr2){
      item=data.frame(chrom=chr1, cend = cyto[i,]$end)
      chr_ends=rbind(chr_ends, item)
    }
    if(i==nrow(cyto)-1){
      item=data.frame(chrom=chr2, cend = cyto[i+1,]$end)
      chr_ends=rbind(chr_ends, item)
    }
  }

  # Find centromere positions
  centro = data.frame()
  cyto %>% filter(stain == "acen") -> cyto_cen
  rows = seq(1, (nrow(cyto_cen)-1), by = 2)
  for (i in rows) {
    chr1 = cyto_cen[i,]$chrom
    start = cyto_cen[i,]$start
    chr2 = cyto_cen[i+1,]$chrom
    end = cyto_cen[i+1,]$end

    if(chr1 != chr2){
      stop("wrong row!")
    }else{
      item=data.frame(chrom=chr1, centro_start = start, centro_end = end)
      centro=rbind(centro, item)
    }
  }

  # Find telomere positions (take first and last band for convenience)
  telo = data.frame()
  cyto %>% filter(start == 0 | lead(start) == 0) -> cyto_telo
  rows = seq(1, (nrow(cyto_telo)-1), by = 2)
  for (i in rows) {
    chr1 = cyto_telo[i,]$chrom
    end1 = cyto_telo[i,]$end
    chr2 = cyto_telo[i+1,]$chrom
    end2 = cyto_telo[i+1,]$start

    if(chr1 != chr2){
      stop("wrong row!")
    }else{
      item=data.frame(chrom=chr1, telo_end1 = end1, telo_end2 = end2)
      telo=rbind(telo, item)
    }
  }

  merge(chr_ends, arm_boudaries, by="chrom") -> chr_end_arm
  merge(chr_end_arm, centro, by="chrom") -> chr_cent
  merge(chr_cent, telo, by="chrom") -> chr_cent_telo

  return(chr_cent_telo)
}


get_fragile_sites <- function(dir){
  files = list.files(dir, "gff3")
  pos_all = data.frame()
  for(i in 1:length(files)){
    # i = 23
    # print(i)
    fname = files[i]
    if(grepl("chrx", fname)) next
    # print(fname)
    fpath = file.path(dir, fname)
    pos = read.gff(fpath)
    pos_all = rbind(pos_all, pos)
  }
  pos_all %>% select(chrom = seqid, start, end, attributes) %>% mutate(chrom = str_replace(chrom, "chr", "")) %>% arrange(as.integer(chrom)) -> pos_sel
  return(pos_sel)
}

fcyto = fname
ref_size = get_cyto_band(fname)

ref_size %>% filter(chrom != 'X') %>% filter(chrom != 'Y') %>% arrange(as.numeric(chrom)) -> ref_sel
write_tsv(ref_sel, fout)


# dir = "data/humcfs_fragile.gff3"
# fs_autosome = get_fragile_sites(dir) 
# fout = "data/fragile_sites.tsv"
# write_tsv(fs_autosome, fout)