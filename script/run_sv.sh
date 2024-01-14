#!/usr/bin/env bash

if [[ ! -d example ]]; then
  mkdir example
fi

cd example

seed=$RANDOM

# change verbose level to get different output
verbose=0

track_all=0  # not keep track of all cells to save computation

chr_prob=2   # chr1 has 2/3 probablity of DSBs whereas the other chromosomes have the same probablity whose sum equals to 1
fchr=../data/hg38_size.tsv # size of the reference genome


num_cell=2  # simulated 2 cells, namely one cell division

n_dsb=50  # introduce 50 DSBs
frac_unrepaired=0  # all DSBs are repaired

n_local_frag=0  # no local fragmentation
circular_prob=0  # no additional formation of circular DNA

# output options
outdir=./
write_rck=1
write_shatterseek=1
write_genome=1
write_sumstats=1


../bin/simsv -n $num_cell --n_dsb $n_dsb --frac_unrepaired $frac_unrepaired --chr_prob $chr_prob --fchr $fchr --track_all $track_all -o $outdir --write_rck $write_rck --write_shatterseek $write_shatterseek  --write_genome $write_genome --write_sumstats $write_sumstats --n_local_frag $n_local_frag --seed $seed  --circular_prob $circular_prob  --verbose $verbose > std_example
