#!/usr/bin/env bash

# This script is used to call the Julia script to run ABC on real data

export JULIA_NUM_THREADS=2

# change paths accordingly
bdir="./fitsc/"
pdir="./CIN_SV/"

ddir=${bdir}/data/
sdir=${bdir}/script/

epsilon=0.5
max_frac_unrepaired=1
max_wgd=1
max_selection_strength=10
max_local_frag=50

while read line; do
nparams=3
# echo $line 
dataset=`echo $line | cut -d' ' -f1`
if [[ $dataset == "sample" ]]; then
continue
fi
num_cell=`echo $line | cut -d' ' -f2`

num_wgd=`echo $line | cut -d' ' -f10`
if [[ $num_wgd == 0 ]]; then 
nparams=$((nparams - 1))
# continue
fi
echo $nparams


rdsb=`echo $line | cut -d' ' -f11`

if [[ $rdsb -lt 50 ]]; then
  max_n_dsb=50
  # continue
else if [[ $rdsb -lt 100 ]]; then
  max_n_dsb=100
  # continue
else if [[ $rdsb -lt 150 ]]; then
  max_n_dsb=150
  # continue
else if [[ $rdsb -lt 200 ]]; then
  max_n_dsb=200
  # continue
else if [[ $rdsb -lt 250 ]]; then
  max_n_dsb=250
  # continue
fi
fi
fi 
fi 
fi 

echo $num_cell
echo $dataset
echo $max_n_dsb
div_break=$num_cell

fout=`ls $bdir/smc/ | grep  "smc_.*${dataset}"`
echo $fout
if [[ -f $fout ]]; then 
echo "file exist"
continue
fi

/usr/bin/time julia $sdir/abc_smc_real.jl $epsilon $max_n_dsb $max_frac_unrepaired $max_wgd $div_break $num_cell "$dataset" $nparams $bdir $pdir $ddir

done < $ddir/patient_nclone_nbp.tsv

# Below is the content for patient_nclone_nbp.tsv

# sample	n	nSV_unique	nSV_clonal	nSV_subclonal	n_sequenced_cells	n_high_quality_cells	cell_type	signature_type	nWGD_clone	rbp	ntc	ntc_all	nsv	frac_tc
# DG1134	7	40	0	40	455	133	HGSC	FBI	7	7	2	2	40	0.025
# DG1197	3	37	9	25	345	115	HGSC	FBI	3	13	0	0	37	0
# SA039	8	40	0	40	2061	878	184-hTert	BASE	2	6	9	9	40	0.1125
# SA1049	8	573	19	544	3564	1283	HGSC	FBI	8	78	136	138	573	0.12041884816753928
# SA1054	5	142	15	132	636	382	184-hTert	TP53-BRCA	4	33	20	20	142	0.07042253521126761
# SA1055	10	18	0	18	582	391	184-hTert	TP53-BRCA	6	2	0	0	18	0
# SA1056	2	75	68	7	758	496	184-hTert	TP53-BRCA	0	7	3	13	75	0.08666666666666667
# SA1091	5	366	68	272	1629	506	HGSC	FBI	5	68	63	83	366	0.1133879781420765
# SA1096	7	326	11	317	3944	802	HGSC	FBI	5	53	72	76	326	0.1165644171779141
# SA1162	8	144	0	144	3176	254	HGSC	FBI	8	21	15	15	144	0.05208333333333333
# SA1180	13	320	13	317	1228	774	HGSC	FBI	13	27	38	38	320	0.059375
# SA1182	2	211	153	58	487	214	HGSC	FBI	2	58	4	14	211	0.03317535545023697
# SA1188	11	37	0	37	1245	472	184-hTert	TP53-BRCA	0	4	7	7	37	0.0945945945945946
# SA1292	4	28	1	27	1289	377	184-hTert	TP53-BRCA	0	9	8	8	28	0.14285714285714285
# SA530	3	123	27	75	743	324	TNBC	FBI	0	38	28	48	123	0.1951219512195122
# SA604	14	385	24	381	6630	2139	TNBC	FBI	14	30	68	73	385	0.0948051948051948
# SA609	16	558	5	542	14008	6033	TNBC	FBI	0	37	57	57	558	0.0510752688172043
# SA610	2	114	62	33	1250	268	TNBC	FBI	2	33	16	39	114	0.17105263157894737
# SA906a	13	38	1	38	1441	650	184-hTert	TP53	2	4	9	9	38	0.11842105263157894
# SA906b	10	39	0	38	1355	984	184-hTert	TP53	0	5	2	2	39	0.02564102564102564
