# Introduction
This program simulates cell divisions with a stochastic birth-death branching process.
It generates double strand breaks and repair them in each cell cycle.
It can output copy number alterations (CNAs) and structural variants (SVs) for all cells in the final population.

![The stochastic cell-cycle model](model.jpeg "The stochastic cell-cycle model of SV generation from DNA repair and replication.")


The full description of the method and its applications are described in: \
Bingxin Lu, Samuel Winnall, William Cross, Chris P. Barnes (2023). Cell-cycle dependent DNA repair and replication unifies patterns of chromosome instability. [bioRxiv](https://doi.org/10.1101/2024.01.03.574048).

The simulated and real data used for the analysis in the above manuscript are available at
[Zenodo](https://doi.org/10.5281/zenodo.10114638).


# Installation
The program has been tested on *nix systems.

## Required library
* [GSL](https://www.gnu.org/software/gsl/)
* [BOOST](https://www.boost.org)

## Compilation
Download the source code and then run `make` in folder "src" to generate the simulation program 'simsv' in a few seconds.


# Usage
Please run `./simsv -h` for all options.

Please run `bash script/run_sv.sh` for an example.



## Input
The program takes an input file containing the size of the reference genome.
Please see hg38_size.tsv in folder "data" for example.


### Important input parameters
* `-n `: number of cells in the simulated population
* `--div_break `: maximum ID of cell division when DSBs occurs
* `--dsb_rate`: mean number or rate of double strand breaks per cell division
* `--n_dsb`: number of double strand breaks introduced in G1, which will be overridden by dsb_rate > 0
* `--frac_unrepaired`: fraction of unrepaired double strand breaks in G1
* `--n_local_frag`: mean number of double strand breaks introduced by local fragmentation during mitosis
* `--prob_wgd`: probability of whole genome doubling
* `--chr_prob`: the types of assigning probability of double strand breaks
* `--fchr_prob`: the file containing the probability of double strand breaks on each chromosome


## Output 
The program will generate a list of summary statistics to the standard output for inference with Approximate Bayesian Computation (ABC), including
* the percentage of genome altered (PGA) for total and haplotype A/B copy number alterations 
* the mean and standard deviation of pairwise divergence for total and haplotype A/B copy number alterations 
* the frequency distribution of breakpoints present in different numbers of cells
* the fraction of cells with WGD
* the fraction of different types of SVs (deletion, duplication, inversions, and intra-chromosomal SVs)


In addition, several files may be generated depending on the output options specified by the user.
* sumStats_sim* : files including summary statistics for the whole simulation 
* sumStats_total_c* : files including summary statistics for each cell 
* sumStats_chrom_c* : files including summary statistics for each chromosome in each cell 
* CNData_*: files including copy numbers for each cell, format compatible with [ShatterSeek](https://github.com/parklab/ShatterSeek)
* SVData_*: files including SVs for each cell, format compatible with [ShatterSeek](https://github.com/parklab/ShatterSeek)
* c*/rck.acnt.tsv: files including copy numbers for each cell, format compatible with [RCK](https://github.com/aganezov/rck)
* c*/rck.scnt.tsv: files including SVs for each cell, format compatible with [RCK](https://github.com/aganezov/rck)

### Data fitting with ABC

Please see scripts in folder "script/abc" for how to run ABC on simulated and real data.


### Visualize output
The script `visualize_sv.R` provides several visualization functions of the data.

Please see the following scripts in folder "script" for how to generated the plots in the manuscript.
* Fig 1: plot_demo_model.R
* Fig 2: parse_bfb.R, check_cmplxy.R, parse_pcn.R
* Fig 3: count_sv.R
* Fig 4: parse_ecdna.R
* Fig 5: abc/check_abc_sc_sim.R
* Fig 6: abc/check_abc_sc_real.R
