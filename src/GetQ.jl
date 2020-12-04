# Functions to get Q matrix
import DSP
using Distributed
using SharedArrays

# get the qubo
# inputs: measured values of head, F matrix, initial head
# outputs: qubo matrix
function qubo(hhat, F, h0)
    #@show F
    Q = F'*F + LinearAlgebra.Diagonal(2*F'*(h0-hhat))
    Q2 = zeros(size(Q))

    for i = 1:size(Q)[1]
        Q2[i,i] = Q[i,i]
        for j = 1:i-1
            Q2[i,j] = Q[i,j] + Q[j,i]
        end
    end
    return Q2'
end

# get the F matrix LxJ matrix
# inputs: forward modeling operator, i index of known l values, initial head, low head value, high head value, number of times i grid is within j grid
# outputs: F matrix used in qubo function
function getF(forwardModel, lvals, k0, klow, khigh, inum)
    # outputs the known values of head
    function h(k0, lvals)
        kinnn = kinn(k0,inum)
        hi = forwardModel(kinnn)
        return hi[lvals]
    end
    hl = h(k0, lvals)
    F = zeros(length(hl), length(k0))
    #@printf("First getF")
    for j = 1:length(k0)
        k1 = copy(k0)
        if k0[j] == klow
            k1[j] = khigh
        else
            k1[j] = klow
        end
        F[:,j] = h(k1, lvals) - hl
    end
    return F
end

# get the F matrix LxJ matrix
# inputs: forward modeling operator, initial head, low head value, high head value, number of times i grid is within j grid
# outputs: F matrix used in qubo function
function getF(forwardModel, k0, klow, khigh, inum)
    # outputs the known values of head
    function h(k0)
        kinnn = kinn(k0,inum)
        hi = forwardModel(kinnn)
        return hi
    end
    hl = h(k0)
    #F = zeros(length(hl), length(k0))
    F = SharedArray{Float64}(length(hl), length(k0))
    # @distributed for j = 1:length(k0)
    #    @sprintf("\n %i \n", j)
    #for j = 1:length(k0)
    #@printf("Second getF")
    #@printf("Using multiple threads")
    #Threads.@threads 
    @sync @distributed for j = 1:length(k0)
        #@show j
        #@printf("\n %i %i \n", j, Threads.threadid())
        k1 = copy(k0)
        if k0[j] == klow
            k1[j] = khigh
        else
            k1[j] = klow
        end
        F[:,j] = h(k1) - hl
    end
    return F
end

# get standard Q value
# inputs: forward modeling operator, guess for permeability, low permeability, high permeability, number of times i grid is within j grid, location of head measurements, measured head values
# output: Q matrix
function getQ(forwardModel, initGuess, kl, kh, inum, lvals, hhat)
    @printf("Calculating Q...\n")
    h0 = geth0(forwardModel, initGuess, kl, kh, inum, lvals)
    kjGuess = b2f(initGuess, kl, kh)    # convert to acutal perm 
    F = getF(forwardModel, lvals, kjGuess, kl, kh, inum)
    Q = qubo(hhat, F, h0)               # get Q matrix
    @printf("Done.\n\n")
    return Q
end


# get standard Q value
# inputs: forward modeling operator, guess for permeability, low permeability, high permeability, number of times i grid is within j grid, location of head measurements, measured head values
# output: Q matrix
function getQ(forwardModel, initGuess, kl, kh, inum, hhat)
    @printf("Calculating Q...\n")
    if typeof(initGuess[1]) != Int
        @printf("Wrong input: getQ \n")
    end
    h0 = geth0(forwardModel, initGuess, kl, kh, inum)
    kjGuess = b2f(initGuess, kl, kh)    # convert to acutal perm 
    F = getF(forwardModel, kjGuess, kl, kh, inum)
    Q = qubo(hhat, F, h0)               # get Q matrix
    @printf("Done.\n\n")
    return Q
end


# zeros out values smaller than a percentage in a matrix
# inputs: the Q matrix, the percent which smaller than that = 0
# output: the threshold Q matrix
function thresholdQ(Q, percent)
    idx = floor(Int, length(Q)*percent)
    threshold = sort(abs.(Q[:]))[idx]
    mask = abs.(Q) .> threshold
    return Q.*mask
end

# zeros out values smaller than a percentage in a matrix
# inputs: head measurements, F matrix, forward model head measurements, the percent which smaller than that = 0
# output: the threshold Q matrix
function thresholdQubo(hhat, F, h0, percent)
    Q1 = qubo(hhat, F, h0)
    Q = thresholdQ(Q1, percent)
    return Q
end

# zeros out values for two qs farther than diff apart
# input: Q matrix, distance between q values
# output: the threshold Q matrix
function thresholdQ2(Q, diff)
    Q2 = copy(Q)
    for i = 1:size(Q, 1)
        for j = 1:size(Q, 2)
            if abs(i-j) > diff
                Q2[i, j] = 0
            end
        end
    end
    return Q2
end

# zeros out values for two qs farther than diff apart
# inputs: head measurements, F matrix, forward model head measurements, distance between q values
# output: the threshold Q matrix
function thresholdQubo2(hhat, F, h0, diff)
    Q1 = qubo(hhat, F, h0)
    Q = thresholdQ2(Q1, diff)
    return Q
end

# add epsilon*identity matrix to Q matrix
# inputs: Q matrix, epsilon value
# output: updated Q matrix
function tikhonovQ(Q, epsilon)
    return Q + epsilon.* LinearAlgebra.I
end

# add epsilon*identity matrix to Q matrix
# inputs: head measurements, F matrix, forward model head measurements, epsilon value
# output: updated Q matrix
function tikhonovQubo(hhat, F, h0, epsilon)
    Q1 = qubo(hhat, F, h0)
    Q = tikhonovQ(Q1, epsilon)
    return Q
end

# add epsilon*identity matrix to Q matrix, where epsilon is a percent of the maximum value in Q
# inputs: Q matrix, percent
# output: updated Q matrix
function tikhonovQ2(Q, percent)
    maxVal = sort(abs.(Q[:]))[1]
    epsilon = percent*maxVal
    return Q + epsilon.* LinearAlgebra.I
end

# add epsilon*identity matrix to Q matrix, where epsilon is a percent of the maximum value in Q
# inputs: head measurements, F matrix, forward model head measurements, epsilon value
# output: updated Q matrix
function tikhonovQubo2(hhat, F, h0, epsilon)
    Q1 = qubo(hhat, F, h0)
    Q = tikhonovQ2(Q1, epsilon)
    return Q
end


# Updating Q value based on how incorrect a solution is
# inputs: Q matrix, forward modeled head, measured head, locations of measured head values, number of times i grid is within j grid, scaling constant
# output: updated Q matrix
function residualQCorrection(Q, h0, hhat, lvals, inum, scaleVal)
    maxVal = maximum(Q)
    diff = abs.(h0 - hhat)
    diff = scaleVal*maxVal ./(diff /maximum(diff))
    interp = weightinterp(size(Q,2), inum, diff, lvals)
    return Q + LinearAlgebra.Diagonal(interp)
end

# create array with interpolated weight values
# inputs: size of returning array, number of times i grid is within j grid, weight values, index of weight values
# output: interpolated array
function weightinterp(sizediag, scale, vals, lidx)
    arr = zeros(sizediag*scale)
    j = 1
    for i =1:sizediag*scale
        if i in lidx
            arr[i] = vals[j]
            j += 1
        end
    end
    #sizeFilt = minimum(diff(lidx))
    filt = DSP.conv(ones(10),ones(10))
    filt = DSP.conv(filt,ones(10))
    filt = filt./maximum(filt)
    foldConvVal = foldConv(arr, filt)
    decimated = foldConvVal[1:scale:end]
    return decimated
end

# convolve with folding edges
# inputs: two arrays to be convolved
# output: convolved array
function foldConv(a, b)
    foldSize = floor(Int, length(b) / 2)
    convArray = DSP.conv(a, b)
    foldArr = convArray[(foldSize+1):(length(convArray)-foldSize)]
    #for i = 1:foldSize
    #    foldArr[i] = foldArr[i] + convArray[i]
    #    foldArr[length(foldArr)+1 - i] = foldArr[length(foldArr)+1 - i] + convArray[length(foldArr)+1-i]
    #end
    return foldArr
end

# Updating Q value based on how incorrect a solution is
# inputs: forward modeled head, measured head, locations of measured head values, number of times i grid is within j grid, scaling constant, F matrix
# output: updated Q matrix
function residualQCorrectionQubo(h0, hhat, lvals, inum, scaleVal, F)
    Q1 =  qubo(hhat, F, h0)
    Q = residualQCorrection(Q1, h0, hhat, lvals, inum, scaleVal)
    return Q
end
