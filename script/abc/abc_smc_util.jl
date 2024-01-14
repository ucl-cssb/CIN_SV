using ApproxBayes
using Distributions
using DelimitedFiles
using Distances
using DataFrames

# constants and functions for running ABC SMC

model=1
n_local_frag=0

dir = "./"
fchr_prob = ""
n_dsb=0 # overload by r_dsb

selection_type=1  # arm level
# growth_type=0   # only birth
growth_type=3  # either birth or death

pair_type=1
prob_correct_repaired=0
only_repair_new=0

track_all=0
chr_prob=0  # random chr
circular_prob=0

verbose=0
bin_level_sumstat=1
write_bin=0
write_rck=0
write_shatterseek=0
write_genome=0
write_sumstats=0
write_selection=0

min_n_dsb = 0
min_frac_unrepaired = 0
min_break = 0
min_selection_strength = 1 
min_local_frag = 0
min_wgd = 0

selection_strength = 1

sstart = 1


function get_simdata(r_dsb, frac_unrepaired, selection_strength, n_local_frag, prob_wgd, seed, model, div_break, num_cell, fbp, fbp_common)
  cmdSim = `$progSim -o $dir -n $num_cell --bin_level_sumstat $bin_level_sumstat --write_bin $write_bin --selection_type $selection_type --growth_type $growth_type --pair_type $pair_type --prob_correct_repaired $prob_correct_repaired --model $model --selection_strength $selection_strength --only_repair_new $only_repair_new --write_selection $write_selection --dsb_rate $r_dsb --n_dsb $n_dsb --div_break $div_break --frac_unrepaired $frac_unrepaired --chr_prob $chr_prob --fchr_prob $fchr_prob --fchr "$fchr" --fbin "$fbin" --fbp "$fbp"  --fbp_common "$fbp_common"  --track_all $track_all --write_rck $write_rck --write_shatterseek $write_shatterseek  --write_genome $write_genome --write_sumstats $write_sumstats --seed $seed --circular_prob $circular_prob --n_local_frag $n_local_frag --prob_wgd $prob_wgd --verbose $verbose`
  res = chomp(read(cmdSim, String))
  arr = readdlm(IOBuffer(res))
  # simdata = arr[:]
  simdata = convert(Array{Float64}, arr[sstart:send])
  # print(simdata)
  return simdata
end


# Get the distance between simulated and real data (read from files) from copy numbers
function getDistance(params, constants, targetdata)
  # call simulation program to generate data
  # print(params)
  r_dsb = params[1]
  if length(params) == 2
    frac_unrepaired = params[2]
    selection_strength = 1
    n_local_frag = 0
    prob_wgd = 0
    div_break = constants[2]
  elseif length(params) == 3
    frac_unrepaired = params[2]
    # n_local_frag = params[3]
    prob_wgd = params[3]
    selection_strength = 1
    n_local_frag = 0
    div_break = constants[2]
  elseif length(params) == 4
    frac_unrepaired = params[2]
    # n_local_frag = params[3]
    prob_wgd = params[3]
    div_break = convert(UInt32, round(params[4]))  
    selection_strength = 1
    n_local_frag = 0      
  else
    frac_unrepaired = 0
    selection_strength = 1
    n_local_frag = 0
    prob_wgd = 0
    div_break = constants[2]
  end
  seed = rand(UInt64, 1)[1]
  model = constants[1]
  num_cell = constants[2]
  fbp = constants[3]
  fbp_common = constants[4]

  simdata = get_simdata(r_dsb, frac_unrepaired, selection_strength, n_local_frag, prob_wgd, seed, model, div_break, num_cell, fbp, fbp_common)

  r = euclidean(targetdata, simdata)

  return r, 1
end
