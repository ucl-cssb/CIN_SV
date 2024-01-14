#!/usr/bin/env bash

# This script is used to call the Julia script to run ABC on simulated data

export JULIA_NUM_THREADS=2


# change paths accordingly
bdir="./fitsc/"
pdir="./CIN_SV/"

ddir=${bdir}/data/
sdir=${bdir}/script/

# simulating "true" data
num_cell=10
dataset="2340225"
epsilon=0.2


r_dsb_target=$1
frac_unrepaired_target=$2
prob_wgd_target=$3

max_n_dsb=$4

# for inference
max_frac_unrepaired=1
max_wgd=1
nparam=3

echo $num_cell    
echo $dataset
echo $max_n_dsb


/usr/bin/time julia $sdir/abc_smc_sim.jl $num_cell "$dataset" $r_dsb_target $frac_unrepaired_target $prob_wgd_target $epsilon $bdir $pdir $ddir $max_n_dsb 