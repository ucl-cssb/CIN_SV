#!/bin/bash

if [[ ! -d example ]]; then
  mkdir example
fi

cd example

seed=$RANDOM

# change verbose level to get different output
verbose=0

track_all=0

chr_prob=2
fchr=../data/hg38_size.tsv # size of the reference genome


num_cell=2

n_dsb=50
frac_unrepaired=0

# output options
write_rck=1
write_shatterseek=1
write_genome=1
write_sumstats=1

n_local_frag=0

circular_prob=0.5

outdir=./

../bin/simsv -n $num_cell --n_dsb $n_dsb --frac_unrepaired $frac_unrepaired --chr_prob $chr_prob --fchr $fchr --track_all $track_all -o $outdir --write_rck $write_rck --write_shatterseek $write_shatterseek  --write_genome $write_genome --write_sumstats $write_sumstats --n_local_frag $n_local_frag --seed $seed  --circular_prob $circular_prob  --verbose $verbose > std_example
