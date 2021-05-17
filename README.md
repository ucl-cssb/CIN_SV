# Introduction
This program simulates cell divisions with a stochastic birth-death branching process.
It generates double strand breaks and repair them in each cell cycle.
It can output copy number alterations and structural variants for all cells in the final population.

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
