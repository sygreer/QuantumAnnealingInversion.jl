import SparseArrays
import HydroThreeQ
import JLD2
import Random
using Printf

#solve 1d groundwater equation assuming the left and right boundary have 0 head, with a spacing of 1 between nodes
#inputs: permeability, recharge
#output: heads
function forwardModel(ks::Array{Float64, 1})
    rs = zeros(length(ks)-1)
    I = Array{Int}(undef, 3 * (length(ks) - 1) + 2)
    J = Array{Int}(undef, 3 * (length(ks) - 1) + 2)
    V = Array{Float64}(undef, 3 * (length(ks) - 1) + 2)
    b = Array{Float64}(undef, length(ks) + 1)
    I[1] = 1; J[1] = 1; V[1] = 1.0; b[1] = length(ks)
    I[2] = length(ks) + 1; J[2] = length(ks) + 1; V[2] = 1.0; b[2] = 0.0
    push!(I, length(ks) + 1); push!(J, length(ks) + 1); push!(V, 1.); b[end] = 0.0
    l = 3
    dx = 1
    for i = 2:length(ks)
        b[i] = -rs[i - 1]
        I[l] = i
        J[l] = i - 1
        V[l] = ks[i - 1] / dx^2
        I[l + 1] = i
        J[l + 1] = i + 1
        V[l + 1] = ks[i] / dx^2
        I[l + 2] = i
        J[l + 2] = i
        V[l + 2] = -(ks[i - 1] + ks[i]) / dx^2
        l += 3
    end
    A = SparseArrays.sparse(I, J, V)
    return A \ b
end

function chooseQ(choice, percent, difference, epsilon, percent2, scaleVal, forwardModel, initGuess, kl, kh, inum, lvals, hhat)
    F = HydroThreeQ.getF(forwardModel, lvals, kjGuess, kl, kh, inum)
    if choice == 1
        Q = HydroThreeQ.getQ(forwardModel, initGuess, kl, kh, inum, lvals, hhat)
        #@printf("getQ\n")
    elseif choice == 2
        Q = HydroThreeQ.thresholdQubo(hhat, F, h0, percent)
        #@printf("thresholdQubo, %f \n", percent)
    elseif choice == 3
        Q = HydroThreeQ.thresholdQubo2(hhat, F, h0, difference)
        #@printf("thresholdQubo2, %f \n", difference)
    elseif choice == 4
        Q = HydroThreeQ.tikhonovQubo(hhat, F, h0, epsilon)
        #@printf("tikhonovQubo, %f \n", epsilon)
    elseif choice == 5
        Q = HydroThreeQ.tikhonovQubo2(hhat, F, h0, percent2)
        #@printf("tikhonovQubo2, %f \n", percent2)
    else
        Q = HydroThreeQ.residualQCorrectionQubo(h0, hhat, lvals, inum, scaleVal, F)
        #@printf("residualQCorrectionQubo, %f \n", scaleVal)
    end
    return Q
end


## Initialize problem
dk = 0.001      # difference between klow and khigh
kl = 1.0        # klow
kh = kl + dk    # khigh
inum = 10       # size of j index in i index

# chose inversion variables
global percent, difference, epsilon, percent2, scaleVal, choice
percent, difference, epsilon, percent2, scaleVal  = 0.9, 10, 1000, 0.9, 1
choice = 1

global dictQ
Random.seed!(87923847289374)

# get initial model
sizeModel = 18
kBin = rand(0:1, sizeModel) # true solution
#kBin = [1 0 1 0 0 1 1 0 1 0 0 0 0 0 0 1]'
#kBin = rand(0:1, sizeModel).*rand(0:1, sizeModel) # true solution ( ~ low pass filter)
kj = HydroThreeQ.b2f(kBin, kl, kh)      # array of known k
ki = HydroThreeQ.kinn(kj, inum)         # interpolate to grid i
htrue = forwardModel(ki)    # get true head values

# sensors 
nSens = floor(Int, length(kBin)*0.25)    # number of sensors
#nSens = 8
lvals = zeros(Int, nSens)
spacing = floor(Int, length(htrue)/nSens)
start = floor(Int, spacing/2)
for i=1:nSens
    lvals[i] = start + (i-1)*spacing
end
#lvals = sort(rand(1:length(htrue), nSens))    # index in i where our sensors are 
hhat = htrue[lvals]                     # measured values of h located at lvals

#initGuess = copy(kBin)
#initGuess = Array(reshape([zeros(Int(length(kBin)/2)) ones(Int(length(kBin)/2))]', length(kBin))) #
#initGuess = rand(0:1, length(kBin)) # init guess random
initGuess = zeros(length(kBin))
h0 = HydroThreeQ.geth0(forwardModel, initGuess, kl, kh, inum, lvals)
kjGuess = HydroThreeQ.b2f(initGuess, kl, kh)    # convert to acutal perm 

# get Q
Q = chooseQ(choice, percent, difference, epsilon, percent2, scaleVal, forwardModel, initGuess, kl, kh, inum, lvals, hhat)
dictQ = Dict("0" => Q, "truth" => kBin)

# check all possible solutions
    
# display parameters
@printf("Number of sensors = %i, I = %i",nSens, length(htrue))
@printf("\nIndex of sensors = "); HydroThreeQ.printArray(lvals)
@printf("\n\n")

# display solution
objSoln = kBin'*Q*kBin
nobj = HydroThreeQ.checkObj(hhat, hhat)
@printf("True value    = "); HydroThreeQ.printArray(kBin)
@printf("\n      # wrong = %i, lobj = %f, nobj = %f\n",sum(abs.(kBin.-kBin)), objSoln[1], nobj)

objGuess = initGuess'*Q*initGuess
nobj = HydroThreeQ.checkObj(h0, hhat)
@printf("\nInitial guess = "); HydroThreeQ.printArray(initGuess)
@printf("\n      # wrong = %i, lobj = %f, nobj = %f\n",sum(abs.(kBin.-initGuess)), objGuess[1], nobj)

@printf("\n\nBegin iterations to solve for permeability:\n\n")

savesoln = initGuess
nobjB = nobj

for o = 1:20
    global soln
    global nobjB
    global savesoln
    if o == 1
        soln2 = copy(initGuess)
        Q = chooseQ(choice, percent, difference, epsilon, percent2, scaleVal, forwardModel, initGuess, kl, kh, inum, lvals, hhat)
        soln, lobj, nobj = HydroThreeQ.oneIterDwave(forwardModel, Q, initGuess, kl, kh, inum, lvals, hhat, 1, kBin, mytoken)
        dictQ[string(o)] = Q
    else 
        soln2 = copy(soln)
        Q = chooseQ(choice, percent, difference, epsilon, percent2, scaleVal, forwardModel, soln, kl, kh, inum, lvals, hhat)
        soln, lobj, nobj = HydroThreeQ.oneIterDwave(forwardModel, Q, soln, kl, kh, inum, lvals, hhat, o, kBin, mytoken)
        dictQ[string(o)] = Q
    end

    # stopping / saving criteria
    if nobj == 0.0
        @printf("Stopping b/c nobjB == 0\n")
        nobjB = nobj
        savesoln = copy(soln)
        break
    end
    if nobj < nobjB
        nobjB = nobj
        savesoln = copy(soln)
    end
    if soln2 == soln
        @printf("Stopping b/c of repeated solutions\n")
        break
    end
end

@printf("True value = ");   HydroThreeQ.printArray(kBin)
@printf("\nSolution   = "); HydroThreeQ.printArray(savesoln)
nwrong = sum(abs.(kBin.-savesoln))
percentc = (100*(sizeModel - nwrong)/sizeModel)
@printf("\n          # wrong = %i, nobj = %.8E, %% correct = %.1f%% \n", nwrong, nobjB, percentc)
@printf("\n")
dictQ["soln"] = savesoln
include("visualize.jl")     # visualize plots
@printf("")
@JLD2.save @sprintf("Q-%i-%i-%i.jld2",sizeModel,inum, dk*10) dictQ
