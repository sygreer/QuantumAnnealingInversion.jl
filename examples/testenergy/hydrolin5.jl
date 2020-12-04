import SparseArrays
import HydroThreeQ
#import JLD2
import Random
import LinearAlgebra
using Printf
import FiniteDifference2D
import PyPlot
## Initialize problem

# UNCOMMENT THIS IF RUNNING INDIVIDUALLY
dk = 1.0     # difference between klow and khigh
kl = 1.0        # klow
kh = kl + dk    # khigh

randseed1 = 9
randseed2=randseed1
sqrtnump = 4 
intval  = 10
numobs = floor(Int,sqrtnump^2/2)

param_chain = 0


nfd = sqrtnump*intval

#Random.seed!(randseed1)
#obsI = rand(1:nfd, numobs)
#obsJ = rand(1:nfd, numobs)
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
    return fin_diff(ks(p))
end

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
ktrue = HydroThreeQ.b2f(kBin, kl, kh)      # array of known k
initGuess = floor.(Int, zeros(length(ktrue)))
#initGuess = map(x->x > 0 ? 0 : 1, randn(sqrtnump * sqrtnump))

hhat = f(ktrue)

inum = 1       # size of j index in i index
h0 = HydroThreeQ.geth0(f, initGuess, kl, kh, inum, 1)

Q = HydroThreeQ.getQ(f, initGuess, kl, kh, inum, hhat)
tw = 15
#Q = HydroThreeQ.thresholdQ2(Q, tw)

include("visualizeqsoln.jl")
