include("abc_smc_util.jl")


# export JULIA_NUM_THREADS=8

# used for simulating target data
num_cell = parse(UInt32, ARGS[1])
dataset = ARGS[2] # "2340225"
r_dsb_target = parse(UInt32, ARGS[3])
frac_unrepaired_target = parse(Float64, ARGS[4])
prob_wgd_target = parse(Float64, ARGS[5])
# for inference
epsilon = parse(Float64, ARGS[6])
bdir = ARGS[7]
pdir = ARGS[8]
ddir = ARGS[9]

max_n_dsb=parse(UInt32, ARGS[10])

# num_cell=10
# dataset="2340225"
# r_dsb_target=10
# frac_unrepaired_target=0.1
# prob_wgd_target=0.1
# seed_target=3824843719262025999

if size(ARGS, 1) > 10
  seed_target=ARGS[11]
else
  seed_target=rand(UInt64, 1)[1]
end
print(seed_target) 

dbreak_target=num_cell 

# for inference
max_frac_unrepaired=1
max_wgd=1
nparam=3

max_break = num_cell - 2
send = 9 + num_cell + 5


dir_smc = bdir * "smc_sim/seed" * string(seed_target) * "/"
if !isdir(dir_smc)
  mkpath(dir_smc)
end

progSim = pdir * "/bin/simsv"
fbin = pdir * "/data/bin_hg19_500K_5776.tsv"
fchr = pdir * "/data/hg19_size.tsv"

suffix =  "_ncell" * string(num_cell) * "_dsb" * string(r_dsb_target) * "_frac" * string(frac_unrepaired_target) * "_wgd" * string(prob_wgd_target)

fbp = ddir * dataset
fbp_common = ""
write_sumstats=1
dir=dir_smc * "target/"
targetdata = get_simdata(r_dsb_target, frac_unrepaired_target, selection_strength, n_local_frag, prob_wgd_target, seed_target, model, dbreak_target, num_cell, fbp, fbp_common)
println(targetdata)
# write out the target data for reference
ftarget = dir_smc * "target" * suffix
writedlm(ftarget, targetdata)
write_sumstats=0
dir="./"



# stuck when using DiscreteUniform(min_n_dsb, max_n_dsb)
setup1 = ABCSMC(getDistance, #simulation function
    1,  # number of parameters
    epsilon,   # target ϵ
    # parse(Float64, epsilon),
    Prior([Uniform(min_n_dsb, max_n_dsb)]), #Prior for each of the parameters
    maxiterations = 1*10^6, #Maximum number of simulations before the algorithm terminates
    nparticles = 500,
    α = 0.3, # used to determine sets of epsilon
    ϵ1 = 10000.0,
    convergence = 0.05,
    constants = [model, num_cell, fbp, fbp_common],
  )


setup2 = ABCSMC(getDistance, #simulation function
  2,  # number of parameters
  epsilon,   # target ϵ
  # parse(Float64, epsilon),
  Prior([Uniform(min_n_dsb, max_n_dsb), Uniform(min_frac_unrepaired, max_frac_unrepaired)]), #Prior for each of the parameters
  maxiterations = 1*10^6, #Maximum number of simulations before the algorithm terminates
  nparticles = 500,
  α = 0.3, # used to determine sets of epsilon
  ϵ1 = 10000.0,
  convergence = 0.05,
  constants = [model, num_cell, fbp, fbp_common],
)


# , Uniform(min_selection_strength, max_selection_strength), Uniform(min_local_frag, max_local_frag)
setup3 = ABCSMC(getDistance, #simulation function
  3,  # number of parameters
  epsilon,   # target ϵ
  # parse(Float64, epsilon),
  Prior([Uniform(min_n_dsb, max_n_dsb), Uniform(min_frac_unrepaired, max_frac_unrepaired), Uniform(min_wgd, max_wgd)]), #Prior for each of the parameters
  maxiterations = 1*10^6, #Maximum number of simulations before the algorithm terminates
  nparticles = 500,
  α = 0.3, # used to determine sets of epsilon
  ϵ1 = 10000.0,
  convergence = 0.05,
  constants = [model, num_cell, fbp, fbp_common],
)

setup32 = ABCSMC(getDistance, #simulation function
  3,  # number of parameters
  epsilon,   # target ϵ
  # parse(Float64, epsilon),
  Prior([Uniform(min_n_dsb, max_n_dsb), Uniform(min_frac_unrepaired, max_frac_unrepaired), Uniform(min_break, max_break)]), #Prior for each of the parameters
  maxiterations = 1*10^6, #Maximum number of simulations before the algorithm terminates
  nparticles = 500,
  α = 0.3, # used to determine sets of epsilon
  ϵ1 = 10000.0,
  convergence = 0.05,
  constants = [model, num_cell, fbp, fbp_common],
)

setup4 = ABCSMC(getDistance, #simulation function
  4,  # number of parameters
  epsilon,   # target ϵ
  # parse(Float64, epsilon),
  Prior([Uniform(min_n_dsb, max_n_dsb), Uniform(min_frac_unrepaired, max_frac_unrepaired), Uniform(min_wgd, max_wgd), Uniform(min_break, max_break)]), #Prior for each of the parameters
  maxiterations = 1*10^6, #Maximum number of simulations before the algorithm terminates
  nparticles = 500,
  α = 0.3, # used to determine sets of epsilon
  ϵ1 = 10000.0,
  convergence = 0.05,
  constants = [model, num_cell, fbp, fbp_common],
)

fout = "smc_sim" * suffix * "_nparam" * string(nparam) * "_epsilon" * string(epsilon) * "_maxDSB" * string(max_n_dsb)

if nparam == 1
  println("inferring 1 parameter")
  smc = runabc(setup2, targetdata, verbose = true, progress = true, parallel = true)
elseif nparam == 2
  println("inferring 2 parameters")
  smc = runabc(setup2, targetdata, verbose = true, progress = true, parallel = true)
elseif nparam == 3
  println("inferring 3 parameters")
  if max_wgd <= 1
    println("inferring WGD")
    fout = fout * "_maxWGD" * string(max_wgd) 
    smc = runabc(setup3, targetdata, verbose = true, progress = true, parallel = true)
  else
    println("inferring cell cycle ID")
    fout = fout * "_maxBreak" * string(max_break) 
    smc = runabc(setup32, targetdata, verbose = true, progress = true, parallel = true)
  end
else
  println("inferring 4 parameters")
  fout = fout * "_maxWGD" * string(max_wgd) * "_maxBreak" * string(max_break) 
  smc = runabc(setup4, targetdata, verbose = true, progress = true, parallel = true)
end 


fout = fout * "_" * dataset

println(dir_smc)
println(fout)
writeoutput(smc, dir=dir_smc, file=fout)