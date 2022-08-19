# install
# require(devtools)
# install_github("parklab/ShatterSeek")


suppressMessages(library(optparse))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(tidyverse))
suppressMessages(library(ShatterSeek))


get_chromothripsis <- function(dir, midfix, to_plot = T) {
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

  # Plot chromothripsis per chromosome (not work on ubuntu)
  plots_chr <- plot_chromothripsis(ShatterSeek_output = chromothripsis, chr = "1", genome = "hg38")
  p <- arrangeGrob(plots_chr[[1]],
    plots_chr[[2]],
    plots_chr[[3]],
    plots_chr[[4]],
    nrow = 4, ncol = 1, heights = c(0.2, .4, .4, .4)
  )
  
  if(to_plot){
    print("saving plot to file")
    pg = plot_grid(p)
    save_plot(fplot, pg, base_height = 6)
    # print("finish saving")
  }

  # print summary
  summ <- chromothripsis@chromSummary
  # print(head(summ))
  # names(chromothripsis@detail)
  # find rows with positions
  summ %>% filter(!is.na(start)) -> res
  # res %>%
  #   t() %>%
  #   View()
  fout <- file.path(dir, paste0("shatterseek_", midfix, ".tsv"))
  write_tsv(res, fout)
  
  return(chromothripsis)
}


##################### Run on simulated data ############################


option_list <- list(
  make_option(c("-d", "--dir"),
    type = "character", default = "",
    help = "input file directory [default=%default]", metavar = "character"
  ),
  make_option(c("-n", "--div_ID"),
    type = "integer", default = 1,
    help = "The number of cell divisions [default=%default]", metavar = "number"
  ),
  make_option(c("-c", "--cell_ID"),
    type = "integer", default = 2,
    help = "The ID of the cell [default=%default]", metavar = "number"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# dir <- "d:/data/SV/sim/test20220425-2-checked"
# cell_ID <- 3
# div_ID <- 1
#
# dir = "d:/data/SV/sim/test20220503/nDSB50_nUn0/r1"
# dir = "/mnt/d/data/SV/sim/test20220503/nDSB50_nUn0/r1"
# cell_ID = 2
# div_ID = 1

dir <- opt$dir
cell_ID <- opt$cell_ID
div_ID <- opt$div_ID

midfix <- paste0("c", cell_ID, "_div", div_ID)


########## detect chromothripsis
cm <- get_chromothripsis(dir, midfix)



