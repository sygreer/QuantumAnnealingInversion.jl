import SparseArrays
import HydroThreeQ
#import JLD2
import Random
import LinearAlgebra
using Printf
import PyPlot
import FiniteDifference2D

## Initialize problem

## UNCOMMENT THIS IF RUNNING INDIVIDUALLY
#dk = 1.0     # difference between klow and khigh
#kl = 1.0        # klow
#kh = kl + dk    # khigh
#
#randseed1=0
#randseed2=0
#sqrtnump = 4 
#intval  = 5
#numobs = floor(Int,sqrtnump^2/2)
#
#param_chain = 0



nfd = sqrtnump*intval

Random.seed!(randseed1)
obsI = rand(1:nfd, numobs)
obsJ = rand(1:nfd, numobs)

#obsI = Int64[]
#obsJ = Int64[]
#n = intval
#for i=1:nfd
#    for j=1:nfd
#        if (i%n +j%n) == 0
#            push!(obsI, floor(Int, i-0.5*intval+1))
#            push!(obsJ, floor(Int, j-0.5*intval+1))
#        end
#    end
#end
numobs = length(obsJ)

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
#kBin2 = rand(0:1, 2*(nfd+1)*nfd) # true solution
#x = HydroThreeQ.b2f(kBin, kl, kh)      # array of known k
#y = fin_diff(x)

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
    #    @printf("ERROR IN f(kBin)\n")
    #    @show(unique(kBin))
    ##    return fin_diff(ks(p))
    #else
    #    @printf("f(kBin) is correct\n")
    #end
    #    return fin_diff(ks(HydroThreeQ.b2f(p, kl, kh)))
    return fin_diff(ks(p))
end

function fnbin(kBin)
    #p = HydroThreeQ.b2f(kBin, kl, kh)
    #return fin_diff(ks(p))
    p = arr2mat(kBin)
    return fin_diff(ks(p))
end

#function g(p)
#    x = ks(p)
#    Jk = fin_diff(x, Val{:J})
#    dkdp = ks(p, Val{:J})
#    return Jk * dkdp
#end


function arr2mat(arr)
    return reshape(arr, sz, sz)
end

function mat2arr(mat)
    return reshape(mat, sz*sz)
end


# chose inversion variables
#global dictQ

#%kBin = rand(0:1, 2*(nfd+1)*nfd) # true solution
sz = sqrtnump 
Random.seed!(randseed2)
kBin = map(x->x > 0 ? 0 : 1, randn(sqrtnump * sqrtnump))
#@show kBinMat
#kBinMat = rand([0,1], sz, sz)
#kBin = mat2arr(kBinMat)
ktrue = HydroThreeQ.b2f(kBin, kl, kh)      # array of known k
initGuess = floor.(Int, zeros(length(ktrue)))

hhat = f(ktrue)

inum = 1       # size of j index in i index
h0 = HydroThreeQ.geth0(f, initGuess, kl, kh, inum, 1)

# check all possible solutions
    
# display parameters
@printf("# sensors = %i, # qubits = %i, # possible solutions = %i",numobs, length(kBin), 2^length(kBin))
#@printf("\nIndex of sensors = "); HydroThreeQ.printArray(1)
@printf("\n\n")

## get Q
@show initGuess
Q = HydroThreeQ.getQ(f, initGuess, kl, kh, inum, hhat)
#dictQ = Dict("0" => Q, "truth" => kBin)

# display solution
@printf("True value    = "); HydroThreeQ.printMat(arr2mat(kBin))
objSoln = kBin'*Q*kBin
nobj = HydroThreeQ.checkObj(hhat, hhat)
@printf("\n      # wrong = %i, lobj = %f, nobj = %f\n",sum(abs.(kBin.-kBin)), objSoln[1], nobj)

@printf("\nInitial guess = "); HydroThreeQ.printMat(floor.(Int, arr2mat(initGuess)))
objGuess = initGuess'*Q*initGuess
nobj = HydroThreeQ.checkObj(h0, hhat)
@printf("\n      # wrong = %i, lobj = %f, nobj = %f\n",sum(abs.(kBin.-initGuess)), objGuess[1], nobj)

@printf("\n\nBegin iterations to solve for permeability:\n\n")

savesoln = initGuess
nobjB = nobj

for o = 1:1
    global iternum
    iternum = o
    global soln
    global nobjB
    global savesoln
    if o == 1
        soln2 = copy(initGuess)

        @show initGuess
        Q = HydroThreeQ.getQ(f, initGuess, kl, kh, inum, hhat)
        soln, lobj, nobj = HydroThreeQ.oneIterDwave(f, Q, initGuess, kl, kh, inum, 1, hhat, 1, kBin, mytoken, flip=true, param_chain=param_chain)
        #dictQ[string(o)] = Q
    else 
        soln2 = copy(soln)
        Q = HydroThreeQ.getQ(f, soln, kl, kh, inum, hhat)
        soln, lobj, nobj = HydroThreeQ.oneIterDwave(f, Q, soln, kl, kh, inum, 1, hhat, o, kBin, mytoken, flip=true, param_chain=param_chain)
        #dictQ[string(o)] = Q
    end

    # stopping / saving criteria
    if nobj == 0.0
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

@printf("True value = ");   HydroThreeQ.printMat(arr2mat(kBin))
@printf("\nSolution   = "); HydroThreeQ.printMat(arr2mat(savesoln))
sizeModel = length(kBin)
nwrong = sum(abs.(kBin.-savesoln))
percentc = (100*(sizeModel - nwrong)/sizeModel)
@printf("\n          # wrong = %i, nobj = %.8E, %% correct = %.1f%% \n", nwrong, nobjB, percentc)
@printf("\n")
#dictQ["soln"] = savesoln
#include("visualize.jl")     # visualize plots
@printf("")
#@JLD2.save @sprintf("Q-%i-%i-%i.jld2",sizeModel,inum, dk*10) dictQ


#kx, ky = FiniteDifference2D.x2ks(kBin .- savesoln, nfd, nfd)
#k = arr2mat(kBin .- savesoln)

#PyPlot.close("all")
#PyPlot.subplot(1, 2, 1)
#PyPlot.imshow(arr2mat(kBin))
#PyPlot.title("truth")
#PyPlot.subplot(1, 2, 2)
#PyPlot.imshow(arr2mat(savesoln))
#PyPlot.title("soln")
#PyPlot.savefig("k.pdf")
