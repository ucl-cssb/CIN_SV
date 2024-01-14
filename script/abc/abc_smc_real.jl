include("abc_smc_util.jl")

# export JULIA_NUM_THREADS=8

epsilon = parse(Float64, ARGS[1])
max_n_dsb = parse(UInt32, ARGS[2])
max_frac_unrepaired = parse(Float64, ARGS[3])
max_wgd = parse(Float64, ARGS[4])
div_break = parse(UInt32, ARGS[5])  # used for predicting mode of evolution, gradual vs punctuated
num_cell = parse(UInt32, ARGS[6])
dataset = ARGS[7]
nparam = parse(UInt32, ARGS[8])
bdir = ARGS[9]
pdir = ARGS[10]
ddir = ARGS[11]

send = 9 + num_cell + 5
if size(ARGS, 1) > 11
  send = parse(UInt32, ARGS[12])
end

max_break = num_cell - 1

# n_local_frag=0  # simple break

# println(send)
# epsilon = 0.15
# max_n_dsb=100
# max_frac_unrepaired=1
# max_wgd = 0.5

# nparam = 3
# num_cell = 7
# dataset = "SA1096"
# num_cell = 11
# dataset = "SA1188"
# num_cell=3
# dataset="SA530"

progSim = pdir * "/bin/simsv"
fbin = pdir * "/data/bin_hg19_500K_5776.tsv"    # excluding chrX
fchr = pdir * "/data/hg19_size.tsv"

ftarget = ddir * "stat/" * dataset * "_sstat.tsv"

fbp = ddir * "bp_subclonal/" * "breakpoints_" * dataset * "_clones.tsv"
fbp_common = ddir * "bp_clonal/" * "breakpoints_" * dataset * "_clones.tsv"


fout = "smc" * "_nparam" * string(nparam) * "_epsilon" * string(epsilon) * "_ncell" * string(num_cell)  * "_maxDSB" * string(max_n_dsb) * "_stat" * string(send)

# targetdata = DataFrame(CSV.File(ftarget, comment="#", header=0))
targetdata = readdlm(ftarget)[:,1][sstart:send]
println(dataset)
println(targetdata)


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
dir_smc = bdir * "smc/"
if !isdir(dir_smc)
  mkpath(dir_smc)
end
println(dir_smc)
println(fout)
writeoutput(smc, dir=dir_smc, file=fout)
