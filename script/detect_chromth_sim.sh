#!/usr/bin/env bash

# This script is used to detect chromothripsis from batch simulations 

# need to change directory paths accordingly
bdir=./
sdir=${bdir}CIN_SV/script
bodir=${bdir}/data_svmodel/Fig3

model=0
chr=1

# num_cell=2
num_cell=3

for div_break in 0 1;
do
t="_ncell${num_cell}_break${div_break}"
for n_unrepaired in 0 0.1 0.3;
do
echo $n_unrepaired
for n_dsb in 10 30;
do
echo $n_dsb
for n_local_frag in 0 10 30;
do
echo $n_local_frag
for i in {1..50..1};
do
echo $i

if [[ $model -eq 0 ]]; then
  odir=${bodir}/test${t}/nDSB${n_dsb}_Un${n_unrepaired}_frag${n_local_frag}/r${i}
else
  odir=${bodir}/test${t}/nDSB${n_dsb}_Un${n_unrepaired}_frag${n_local_frag}_model${model}/r${i}
fi

echo $odir
cd $odir

# used when there are two cells
# Rscript  ${sdir}/detect_chromothripsis.R  -d $odir -n 1 -c 2 -r $chr
# Rscript  ${sdir}/detect_chromothripsis.R  -d $odir -n 1 -c 3 -r $chr

Rscript  ${sdir}/detect_chromothripsis.R  -d $odir -n 2 -c 4 -r $chr
Rscript  ${sdir}/detect_chromothripsis.R  -d $odir -n 2 -c 5 -r $chr

done
done
done
done
done
