# Introduction
This program simulates cell divisions with a stochastic birth-death branching process.
It generates double strand breaks and repair them in each cell cycle.
It can output copy number alterations and structural variants (SVs) for all cells in the final population.

# Installation
## Required library
* [GSL](https://www.gnu.org/software/gsl/)
* [boost](https://www.boost.org)

## Compilation
Download the source code and then run `make` in folder "src" to generate the simulation program 'simsv'.


# Usage
Please run `./simsv -h` for all options.

Please run `bash script/run_sv.sh` for an example.

The results can be visualized by script/visualize_sv.R.

## Important parameters
`-n `: number of cells in the simulated population
`--div_break `: maximum ID of cell division when DSBs occurs
`--n_dsb`: number of double strand breaks introduced in G1
`--frac_unrepaired`: fraction of unrepaired double strand breaks in G1
`--n_local_frag`: mean number of double strand breaks introduced by local fragmentation during mitosis
`--chr_prob`: the types of assigning probability of double strand breaks
`--fchr_prob`: the file containing the probability of double strand breaks on each chromosome


## Output files
sumStats_ : files including summary statistics for each cell and the whole simulation 
CNData_*: files including copy numbers for each cell, format compatible with [ShatterSeek](https://github.com/parklab/ShatterSeek)
SVData_*: files including SVs for each cell, format compatible with [ShatterSeek](https://github.com/parklab/ShatterSeek)
c*/rck.acnt.tsv: files including copy numbers for each cell, format compatible with [RCK](https://github.com/aganezov/rck)
c*/rck.scnt.tsv: files including SVs for each cell, format compatible with [RCK](https://github.com/aganezov/rck)

## Visualize output
Please use functions in `visualize_sv.R`.
