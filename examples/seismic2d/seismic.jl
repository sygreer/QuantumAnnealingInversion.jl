import SparseArrays
import HydroThreeQ
import Random
import LinearAlgebra
using Printf
import FindQUBO
using PyPlot
using Statistics
using SeismicTools


########################
### HELPER FUNCTIONS ###
########################

function arr2mat(arr)
    return reshape(arr, sz2, sz)
end

function mat2arr(mat)
    return reshape(mat, sz*sz2)
end

function getWaveSeismicTools(k)
    # get velocity model
    modl = arr2mat(k)
    velocity = HydroThreeQ.twodinterp(modl, intval)

    # make synthetic truth model
    params = SeismicTools.makeModel(velocity, rho, dz, dx, dt, tmax)

    # get recording
    rec = SeismicTools.getShotRecord(src, irz, irx, params)
end

function f(k)
    rec = getWaveSeismicTools(k)
    return rec.p[:]
end

function forwardModelq(q)
    kjGuess = HydroThreeQ.b2f(q, kl, kh)    # convert to acutal perm
    mod = f(kjGuess)
    return sum(abs.(mod.^2 .- hhat.^2))
end

function plotIter(iternum)
    PyPlot.close("all")
    PyPlot.figure(figsize=(15,10))
    interpsavesoln = HydroThreeQ.twodinterp(arr2mat(soln),intval)
    PyPlot.subplot(2, 2, 2)
    PyPlot.imshow(interpsavesoln, cmap="summer")
    PyPlot.scatter([obsI.-1], [obsJ.-1], s=100, marker="v", c="r", clip_on=false)
    PyPlot.title("soln")
    PyPlot.xlim([size(interpsavesoln)[2]-0.5, 0-0.5])
    PyPlot.ylim([size(interpsavesoln)[1]-0.5, 0-0.5])
    PyPlot.xticks([])
    PyPlot.yticks([])

    interpkBin = HydroThreeQ.twodinterp(arr2mat(kBin), intval)
    PyPlot.subplot(2, 2, 1)
    PyPlot.imshow(interpkBin, cmap="summer")
    PyPlot.scatter([obsI.-1], [obsJ.-1], s=100, marker="v", c="r", clip_on=false)
    PyPlot.title("truth")
    PyPlot.xlim([size(interpsavesoln)[2]-0.5, 0-0.5])
    PyPlot.ylim([size(interpsavesoln)[1]-0.5, 0-0.5])
    PyPlot.xticks([])
    PyPlot.yticks([])

    PyPlot.suptitle(@sprintf("iteration # %i",iternum))

    PyPlot.subplot(2,2,(3,4))

    PyPlot.xticks([])
    PyPlot.yticks([])
    nwrong = sum(abs.(kBin.-savesoln))
    percentc = (100*(sizeModel - nwrong)/sizeModel)
    PyPlot.text(0.5, 0.85, @sprintf("# sensors = %i",numobs ), fontsize=12, horizontalalignment="center")
    PyPlot.text(0.5, 0.80, @sprintf("# qubits = %i",length(kBin) ), fontsize=12, horizontalalignment="center")
    PyPlot.text(0.5, 0.75, @sprintf("# possible solutions = %i",2^length(kBin) ), fontsize=12, horizontalalignment="center")
    PyPlot.text(0.5, 0.70, @sprintf("# iters = %i",iternum), fontsize=12, horizontalalignment="center")
    
    PyPlot.text(0.5, 0.60, @sprintf("non-linear objective function = %f",nobj ), fontsize=12, horizontalalignment="center")
    PyPlot.text(0.5, 0.55, @sprintf("percent correct = %.1f%%",percentc ), fontsize=12, horizontalalignment="center")
    PyPlot.text(0.5, 0.50, @sprintf("number wrong = %i",nwrong ), fontsize=12, horizontalalignment="center")
    
    PyPlot.text(0.5, 0.40, @sprintf("random seed = %i",randseed), fontsize=12, horizontalalignment="center")

    PyPlot.text(0.5, 0.25, @sprintf("u_l = %.2f, u_h = %.2f",kl, kh), fontsize=12, horizontalalignment="center")
    
    PyPlot.savefig(@sprintf("%02d.png",iternum))
end




################################
### INITIALIZE PROBLEM SETUP ###
################################

# seed for model
randseed=32             # seed for model

#### velocity values
#kl = 4250       # klow
#dk = 500        # difference between klow and khigh

kl = 2000       # klow
dk = 4000        # difference between klow and khigh
kh = kl + dk    # khigh

# size of domain
sz = 10 
sz2 = 5 

#sz = 4 
#sz2 = 3 

intval = 10

# param chain value
param_chain=2.4

# get receivers
numobs = 7
hlfobs = Int((numobs-1)/2)
szdomain = sz*intval
center = Int(floor(szdomain/2))
difobs = Int(floor(szdomain/numobs))
obsI = Array(range(center-hlfobs*difobs, step=difobs, length=numobs))
obsJ = (1)*Int.(ones(size(obsI)))

# make source 
kstart = zeros(sz * sz2)
kv = HydroThreeQ.b2f(kstart, kl, kh)      # array of known k
modl = arr2mat(kv)
velocity = HydroThreeQ.twodinterp(modl, intval)

# get spacing in vertical and horizontal direction, in meters
dz = 10 # m
dx = 10 # m

# get time sampling
dt = 0.001      # sampling interval
tmax = 0.401    # time of recording

# size of domain 
nz = size(velocity, 1)
nx = size(velocity, 2)

# get density model (constant density)
rho = 2.0 * ones(nz, nx);  # kg/m^3

# make synthetic truth model for params...
params = SeismicTools.makeModel(velocity, rho, dz, dx, dt, tmax)

# location of source
pz = 2                  # z direction
px = floor(Int,nx/2)    # x direction centered

# get source
fr = 25                 # frequency of source wavelet
src = SeismicTools.makeRicker(pz, px, fr, params, startTime=0.2) 

# rename receivers
irx = obsI
irz = obsJ

########################
### FORWARD MODELING ###
########################

klow = kl
khigh = kh



########################
### SET UP INVERSION ###
########################

# get truth model
Random.seed!(randseed)
kBin = map(x->x > 0 ? 0 : 1, randn(sz * sz2))
sizeModel = length(kBin)
ktrue = HydroThreeQ.b2f(kBin, kl, kh)      # array of known k

# get initial model
initGuess = floor.(Int, zeros(length(ktrue)))


# visualize velocity model
PyPlot.close("all")
PyPlot.figure(figsize=(5,5))
interpsavesoln = HydroThreeQ.twodinterp(arr2mat(ktrue),intval)
PyPlot.imshow(interpsavesoln, cmap="summer")
PyPlot.scatter([obsI.-1], [obsJ.-1], s=100, marker="v", c="r", clip_on=false)
PyPlot.title("soln")
PyPlot.xlim([size(interpsavesoln)[2]-0.5, 0-0.5])
PyPlot.ylim([size(interpsavesoln)[1]-0.5, 0-0.5])
PyPlot.xticks([])
PyPlot.yticks([])



@printf("Generating truth data...\n")
hhat = f(ktrue)
## add sigma = 0.1 gaussian noise
#hhat = hhat .+ 0.1 .*randn(size(hhat)) .* std(hhat)
@printf("Done.\n\n")

inum = 1       # size of j index in i index
@printf("Getting initial guess...\n")
h0 = HydroThreeQ.geth0(f, initGuess, kl, kh, inum, 1)
@printf("Done.\n\n")

# get Q
Qsave = HydroThreeQ.getQ(f, initGuess, kl, kh, inum, hhat)
Q = copy(Qsave)

# display parameters
@printf("# sensors = %i, # qubits = %i, # possible solutions = %i",numobs, length(kBin), 2^length(kBin))
@printf("\n\n")

# display solution
@printf("True value    = "); HydroThreeQ.printMat(arr2mat(kBin))
objSoln = kBin'*Qsave*kBin
nobj = HydroThreeQ.checkObj(hhat, hhat)
@printf("\n      # wrong = %i, lobj = %f, nobj = %f\n",sum(abs.(kBin.-kBin)), objSoln[1], nobj)


@printf("\nInitial guess = "); HydroThreeQ.printMat(floor.(Int, arr2mat(initGuess)))
objGuess = initGuess'*Qsave*initGuess
nobj = HydroThreeQ.checkObj(h0, hhat)
@printf("\n      # wrong = %i, lobj = %f, nobj = %f\n",sum(abs.(kBin.-initGuess)), objGuess[1], nobj)

@printf("\n\nBegin iterations to solve for permeability:\n\n")

savesoln = initGuess
soln = initGuess
nobjB = nobj
#plotIter(0)

########################
### START ITERATIONS ###
########################

global nwrongarr
nwrongarr = NaN.*zeros(20+1)
nwrongarr[1] = sum(abs.(kBin.-initGuess))
for o = 1:20
    global iternum
    iternum = o
    global soln
    global nobjB
    global savesoln
    global objVals
    global allsolns
    global lobj
    if o == 1
        soln2 = copy(initGuess)
        @printf("init     = ")
        HydroThreeQ.printArray(initGuess)
        @printf("\n\n")
        Q = HydroThreeQ.getQ(f, initGuess, kl, kh, inum, hhat)
        soln, lobj, nobj = HydroThreeQ.oneIterDwaveLinearnl(f, Q, initGuess, kl, kh, inum, 1, hhat, 1, kBin, mytoken, flip=true, param_chain=param_chain, retSolns=false)
    else 
        soln2 = copy(soln)
        Q = HydroThreeQ.getQ(f, soln, kl, kh, inum, hhat)
        soln, lobj, nobj = HydroThreeQ.oneIterDwaveLinearnl(f, Q, soln, kl, kh, inum, 1, hhat, o, kBin, mytoken, flip=true, param_chain=param_chain, retSolns=false)
    end
    ##########################################
    ### VISUALIZATIONS FOR EVERY ITERATION ###
    ##########################################

    # stopping / saving criteria
    global nobj
    epsilon = 0
    nwrongarr[o] = sum(abs.(kBin.-savesoln))

    if nobj <= epsilon # == 0.0
        @printf("Stopping b/c nobjB == 0\n")
        nobjB = nobj
        savesoln = copy(soln)
        nwrongarr[o+1] = sum(abs.(kBin.-savesoln))
        break
    end
    #@show nobj, nobjB
    if nobj < nobjB
        @printf("NEW MINIMUM\n\n")
        nobjB = nobj
        savesoln = copy(soln)
    else    # if it's not better than the last, keep the last one.
        soln = savesoln
    end

    nwrongarr[o+1] = sum(abs.(kBin.-savesoln))
    #plotIter(o)
    if soln2 == soln
        @printf("Stopping b/c of repeated solutions\n")
        break
    end
end

#########################
### DISPLAY SOLUTIONS ###
#########################

@printf("True value = ");   HydroThreeQ.printMat(arr2mat(kBin))
@printf("\nSolution   = "); HydroThreeQ.printMat(arr2mat(savesoln))
sizeModel = length(kBin)
nwrong = sum(abs.(kBin.-savesoln))
percentc = (100*(sizeModel - nwrong)/sizeModel)
@printf("\n          # wrong = %i, nobj = %.8E, %% correct = %.1f%% \n", nwrong, nobjB, percentc)
@printf("\n")
#include("visualize.jl")     # visualize plots
#savefig("visualize.png")
@printf("")

