import SparseArrays
import HydroThreeQ
import Random
import LinearAlgebra
using Printf
import FiniteDifference2D
import FindQUBO
using Statistics

################################
### INITIALIZE PROBLEM SETUP ###
################################

dk = 1.0     # difference between klow and khigh
kl = 1.0        # klow
kh = kl + dk    # khigh

#kl = 10e-8        # klow
#kh = 10e-7
#kh = 10e-2
#dk = kh - kl


# good seeds to look at include: 
#   5: concerges almost immedaitely. 
#   8: converges in 1-2 iters
#   9: converges almost immediately

#randseed1=0             # seed for model
randseed1=10
randseed2=randseed1     # seed for receiver locations (only necessary for randomly placed receivers)
sqrtnump = 7 
intval  = 10
numobs = floor(Int,sqrtnump^2/2)

param_chain=2.4

nfd = sqrtnump*intval

###########################
### INITIALIZE SENSORS ####
###########################

## RANDOM SENSORS
#Random.seed!(randseed1)
#obsI = rand(1:nfd, numobs)
#obsJ = rand(1:nfd, numobs)

## GRID SENSOR PATTERN
#obsI = Int64[]
#obsJ = Int64[]
#n = intval*2
#for i=1:nfd
#    for j=1:nfd
#        if (i%n +j%n) == 0
#            push!(obsI, floor(Int, i-0.5*intval+1))
#            push!(obsJ, floor(Int, j-0.5*intval+1))
#        end
#    end
#end
#numobs = length(obsJ)


# CHECKERBOARD SENSOR PATTERN
obsI = Int64[]
obsJ = Int64[]
n = intval*2
for i=1:sqrtnump
    for j=1:sqrtnump
        if (i%2) == 0 && (j%2-1) == 0
            push!(obsI, floor(Int, i*intval-0.5*intval))
            push!(obsJ, floor(Int, j*intval-0.5*intval))
        elseif (i%2-1) == 0 && (j%2) == 0
            push!(obsI, floor(Int, i*intval-0.5*intval))
            push!(obsJ, floor(Int, j*intval-0.5*intval))
        end
    end
end
numobs = length(obsJ)

#obsI = [1, 3, 5, 7, 2, 4, 6, 1, 3, 5, 7, 2, 4, 6].*10 .-5
#obsJ = [1, 1, 1, 1, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7].*10 .-5

#obsI = [1, 3, 5, 7, 9, 2, 4, 6, 8, 1, 3, 5, 7, 9, 2, 4, 6, 8, 1, 3, 5, 7, 9].*10 .-5
#obsJ = [1, 1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 5, 7, 7, 7, 7, 9, 9, 9, 9, 9].*10 .-5


numobs = length(obsJ)

########################
### FORWARD MODELING ###
########################

weights = ones(numobs)
u_n(x, y) = 0
u_d(x, y) = y
f(x, y) = 0
a, b, c, d = 0, 1, 0, 1
function fin_diff(x)
    global ky, kx, kyi, kxi
    kx, ky = FiniteDifference2D.x2ks(x, nfd, nfd)
    u, jacobian = FiniteDifference2D.solvedsresiduals(obsI, obsJ, weights, zeros(numobs), kx, ky, u_d, u_n, f, nfd, nfd, a, b, c, d)
    retval = Array{Float64}(undef, numobs)
    for i = 1:numobs
        retval[i] = u[obsI[i], obsJ[i]]
    end
    return retval
end

klow = kl
khigh = kh
@generated function ks(p, style::Type{Val{T}}=Val{1}) where {T}
    if T == :J
        kdef = :(ks = zeros(Float64, 2 * nfd * (nfd + 1), length(p)))
        setkscode = quote
            ks[i + (j - 1) * (nfd + 1), pix +  (pjx - 1) * sqrtnump] = 1.
            ks[(nfd + 1) * nfd + j + (i - 1) * nfd, piy + (pjy - 1) * sqrtnump] = 1.
        end
    else
        kdef = :(ks = Array{Float64}(undef, 2 * nfd * (nfd + 1)))
        setkscode = quote
            ks[i + (j - 1) * (nfd + 1)] = newp[pix, pjx]
            ks[(nfd + 1) * nfd + j + (i - 1) * nfd] = newp[piy, pjy]
        end
    end
    q = quote
        #ks = Array(Float64, 2 * nfd * (nfd + 1))
        $kdef
        newp = reshape(p, sqrtnump, sqrtnump)
        for i = 1:nfd + 1
            for j = 1:nfd
                pix = ceil(Int, i * sqrtnump / (nfd + 1))
                pjx = ceil(Int, j * sqrtnump / nfd)
                pjy = ceil(Int, i * sqrtnump / (nfd + 1))
                piy = ceil(Int, j * sqrtnump / nfd)
                $setkscode
            end
        end
        return ks
    end
end

function f(k)
    p = arr2mat(k)
    if typeof(k[1]) == Int
        @printf("ERROR: f(k) doesn't take binary input\n")
    end
    return fin_diff(ks(p))
end

function forwardModelq(q)
    kjGuess = HydroThreeQ.b2f(q, kl, kh)    # convert to acutal perm
    mod = f(kjGuess)
    return sum(abs.(mod.^2 .- hhat.^2))
end

function arr2mat(arr)
    return reshape(arr, sz, sz)
end

function mat2arr(mat)
    return reshape(mat, sz*sz)
end


########################
### SET UP INVERSION ###
########################

sz = sqrtnump 
Random.seed!(randseed2)
kBin = map(x->x > 0 ? 0 : 1, randn(sqrtnump * sqrtnump))
ktrue = HydroThreeQ.b2f(kBin, kl, kh)      # array of known k
initGuess = floor.(Int, zeros(length(ktrue)))

hhat = f(ktrue)
## add sigma = 0.1 gaussian noise
#hhat = hhat .+ 0.05 .*randn(size(hhat)) .* std(hhat)

inum = 1       # size of j index in i index
h0 = HydroThreeQ.geth0(f, initGuess, kl, kh, inum, 1)

# check all possible solutions
    
# display parameters
@printf("# sensors = %i, # qubits = %i, # possible solutions = %i",numobs, length(kBin), 2^length(kBin))
@printf("\n\n")

## get Q
@show initGuess
Qsave = HydroThreeQ.getQ(f, initGuess, kl, kh, inum, hhat)
Q = copy(Qsave) #HydroThreeQ.getQ(f, initGuess, kl, kh, inum, hhat)
#Qsave, A, qubocoeffs, outputs = FindQUBO.getqubo(forwardModelq, Int.(initGuess))
#include("visualizeqsoln.jl")

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
nobjB = nobj

########################
### START ITERATIONS ###
########################


global nwrongarr
nwrongarr = NaN.*zeros(20+1)
nwrongarr[1] = sum(abs.(kBin.-initGuess))
nwrongarr[1] = sum(abs.(kBin.-savesoln))
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

        @show initGuess
        Q = HydroThreeQ.getQ(f, initGuess, kl, kh, inum, hhat)
        #Q = HydroThreeQ.tikhonovQ(Q, 0.1)
        #Q, A, qubocoeffs, outputs = FindQUBO.getqubo(forwardModelq, Int.(initGuess))
        #soln, lobj, nobj, allsolns, objVals = HydroThreeQ.oneIterDwaveLinearnl(f, Q, initGuess, kl, kh, inum, 1, hhat, 1, kBin, mytoken, flip=true, param_chain=param_chain, retSolns=true)
        soln, lobj, nobj = HydroThreeQ.oneIterDwaveLinearnl(f, Q, initGuess, kl, kh, inum, 1, hhat, 1, kBin, mytoken, flip=true, param_chain=param_chain, retSolns=false)

    else 
        soln2 = copy(soln)
        Q = HydroThreeQ.getQ(f, soln, kl, kh, inum, hhat)
        #Q = HydroThreeQ.tikhonovQ(Q, 0.1)
        #Q, A, qubocoeffs, outputs = FindQUBO.getqubo(forwardModelq, Int.(soln))
        #soln, lobj, nobj, allsolns, objVals = HydroThreeQ.oneIterDwaveLinearnl(f, Q, soln, kl, kh, inum, 3, hhat, o, kBin, mytoken, flip=true, param_chain=param_chain, retSolns=true)
        #soln, lobj, nobj = HydroThreeQ.oneIterDwaveLinearnl(f, Q, soln, kl, kh, inum, 3, hhat, o, kBin, mytoken, flip=true, param_chain=param_chain, retSolns=true)
        soln, lobj, nobj = HydroThreeQ.oneIterDwaveLinearnl(f, Q, soln, kl, kh, inum, 1, hhat, o, kBin, mytoken, flip=true, param_chain=param_chain, retSolns=false)

    end
    ##########################################
    ### VISUALIZATIONS FOR EVERY ITERATION ###
    ##########################################
    
    #include("visualizesolnsiter.jl")            # only shows solutions for each iteration, similar to visualize.jl
    #include("visualizetheorsolnsiter.jl")       # shows all theoretical solutions (only use for small sqrtnump)
    #include("visualizealliter.jl")              # all visualizations (only use for small sqrtnump)
    #include("visualizelargeniter.jl")           # only returns d-wave solns


    # stopping / saving criteria
    #
    #
    epsilon = 0#4e-3
    nwrongarr[o+1] = sum(abs.(kBin.-soln))
    if nobj  <= epsilon
        @printf("Stopping b/c nobjB == 0\n")
        nobjB = nobj
        savesoln = copy(soln)
        break
    end
    #@show nobj, nobjB
    if nobj < nobjB
        #@printf("NEW MINIMUM\n\n\n\n")
        nobjB = nobj
        savesoln = copy(soln)
    end

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
@printf("")

