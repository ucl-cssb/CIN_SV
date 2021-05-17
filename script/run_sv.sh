if [[ ! -d test ]]; then
  mkdir test
  cd test

  n_cycle=1
  type_dsb=2
  n_unrepaired=5
  dsb_rate=20
  fchr=../data/hg38_size.tsv
  verbose=2
  ../bin/simsv -n $n_cycle --type_dsb $type_dsb --fchr $fchr --verbose $verbose > std
fi
