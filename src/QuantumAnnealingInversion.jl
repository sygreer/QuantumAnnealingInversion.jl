module HydroThreeQ

import LinearAlgebra
using ThreeQ
using Printf
using Debugger

include("GetQ.jl")
include("HelperFunctions.jl")

# get head measurements
# inputs: forward modeling operator, binary value of permeability, low permeability, high permeability, number of times i grid is within j grid, index of head measurements
# output: measured head values
function geth0(forwardModel, soln, kl, kh, inum, lvals)
    kjGuess = b2f(soln, kl, kh)    # convert to acutal perm 
    kiGuess = kinn(kjGuess, inum)       # interpolate
    hGuess = forwardModel(kiGuess)      # forward model to get head
    #if !isa(hGuess,Array)
    #    hGuess = hGuess.p[:]
    #end
    h0 = hGuess[lvals]                  # get initial head value
    return hGuess
end

# get head measurements
# inputs: forward modeling operator, binary value of permeability, low permeability, high permeability, number of times i grid is within j grid, index of head measurements
# output: measured head values
function geth0(forwardModel, soln, kl, kh, inum)
    if typeof(soln[1])!= Int
        @printf("Wrong input: get h0\n")
    end
    kjGuess = b2f(soln, kl, kh)    # convert to acutal perm 
    kiGuess = kinn(kjGuess, inum)       # interpolate
    hGuess = forwardModel(kiGuess)      # forward model to get head
    return hGuess
end

# inputs: binary representation of permeability, low and high permeability
# output: permeability
function b2f(kb, klow, khigh)
    ks = Array{Float64}(undef, length(kb))
    for i = 1:length(kb)
        if kb[i] == 1
            ks[i] = khigh
        else
            ks[i] = klow
        end
    end
    return ks
end

# nearest neighbor interpolation
# inputs: k values on j grid, number of times i grid is within j grid
# output: nn interpolated k values on i grid
function kinn(kj, inum)
    ki = ones(inum*length(kj))  # kj nn interp to grid i
    for i=1:length(ki)          # convert kj to ki
        val = ceil(Int,i/inum)
        ki[i] = kj[val]
    end
    return ki
end

# interpolate in 2d
function twodinterp(arr, intval)
    retarr = []
    for i = 1:size(arr)[1]
        arri = kinn(arr[i,: ],intval)'
        for j = 1:intval-1
            arri = [arri; kinn(arr[i,: ],intval)']
        end
        if i==1
            retarr = arri
        else
            retarr = [retarr; arri]
        end
    end
    return retarr
end


# get all possible solutions and their value
# input: qubo matrix
# output: all permutations, the objective function value of that permutation, the indices of permutations that minimize the objective function
function getPermutations(Q)
    N = size(Q, 1)
    perms = 2^N-1
    padTo = length(digits(perms, base=2))
    permsValue = zeros(padTo, perms+1)
    objVal = zeros(perms+1)
    for i = 0:perms
        val = digits(i, base=2)
        padVal =  [val; zeros(padTo - length(val) )]
        permsValue[:,i+1] = padVal
        objVal[i+1] = padVal'*Q*padVal
    end
    minVal = minimum(objVal)
    idxMinVal = findall(x->x==minVal, objVal)
    return permsValue, objVal, idxMinVal
end

# get the answers from the D-Wave
# input: Q matrix, D-Wave token
# output: kbs solutions
function getDwaveAnswer(Q, mytoken; param_chain = 0)
    num_reads = 10000   # maximum number 10,000
    kbs = Array{Int}(undef, size(Q,2) , num_reads)
    function finishsolve_helper(m, embans, emb)
        global embansg
        global embg
        embansg = embans
        embg = emb
        #answer = 0.5 * (ThreeQ.DWQMI.unembedanswer(embans["solutions"], emb, broken_chains="weighted_random")' .+ 1)#convert from ising to qubo
        
        #h = LinearAlgebra.diag(Q)
        #Q2 = Q - LinearAlgebra.Diagonal(h)
        #j = ThreeQ.matrix2Q(ThreeQ.DWQMI.qubo2ising(Q2))

        dict = ThreeQ.matrix2Q(ThreeQ.DWQMI.qubo2ising(Matrix(Q)))
        h = zeros(size(Q,1))
        j = Dict()
        for ((n1, n2), v) in dict
            if n1 == n2
                h[n1+1] = v
            else
                j[(n1, n2)] = v
            end
        end
        #h =  
        #j = 
        answer = 0.5 * (ThreeQ.DWQMI.unembedanswer(embans["solutions"], emb, broken_chains="minimize_energy", h=h, j=j)' .+ 1)#convert from ising to qubo
        i = 1
        for j = 1:size(answer, 2)
            occurrences = embans["num_occurrences"][j]
            for l = 1:occurrences
                kbs[:, i] = answer[:, j]
                i += 1
            end
        end
    end
    #dwsolver2 = ThreeQ.DWQMI.getlownoise(mytoken2)
    #x2 = ThreeQ.solvesapi!(Q; token=mytoken, auto_scale=true, solver=dwsolver2, num_reads=num_reads, finishsolve_helper=finishsolve_helper, reuse_embedding=true, param_chain_factor=1/6)
    #x2 = ThreeQ.solvesapi!(Q; token=mytoken, auto_scale=true, solver=dwsolver2, num_reads=num_reads, finishsolve_helper=finishsolve_helper, reuse_embedding=true, postprocessing="optimization", param_chain=15)
    #kbs2 = copy(kbs)
    
    #dwsolver = ThreeQ.DWQMI.getvfyc(mytoken)
    dwsolver = ThreeQ.DWQMI.getdw2q(mytoken)

    if param_chain==0
        param_chain = maximum(abs.(Q))/20
    else
        param_chain = maximum(abs.(Q))/param_chain
    end

    x = ThreeQ.solvesapi!(Q; token=mytoken, auto_scale=true, solver=dwsolver, num_reads=num_reads, finishsolve_helper=finishsolve_helper, reuse_embedding=true, param_chain=param_chain)

    #dwsolver = ThreeQ.DWQMI.defaultsolver
    #x = ThreeQ.solvesapi!(Q; token=mytoken, auto_scale=false, solver=dwsolver, num_reads=num_reads, finishsolve_helper=finishsolve_helper, reuse_embedding=true)
    return kbs#, kbs2
end

# get the answers from the D-Wave
# input: Q matrix, D-Wave token
# output: kbs solutions
function getDwaveAnswerqbsolv(Q, mytoken)
    num_reads = 1000
    kbs = Array{Int}(undef, size(Q,2) , num_reads)
    #function finishsolve_helper(m, embans, emb)
    #    global embansg
    #    global embg
    #    embansg = embans
    #    embg = emb
    #    answer = 0.5 * (ThreeQ.DWQMI.unembedanswer(embans["solutions"], emb)' .+ 1)#convert from ising to qubo
    #    i = 1
    #    for j = 1:size(answer, 2)
    #        occurrences = embans["num_occurrences"][j]
    #        for l = 1:occurrences
    #            kbs[:, i] = answer[:, j]
    #            i += 1
    #        end
    #    end
    #end

    modelargs=("lin_model", "laptop", "c4-sw_sample", "workingdir", "c4")
	model = ThreeQ.Model(modelargs...)
    @defvar model q[1:size(Q)[1]]

	for i = 1:size(Q, 1)
		@addterm model Q[i, i] * q[i]
		for j = 1:i - 1
			@addterm model (Q[i, j] + Q[j, i]) * q[i] * q[j]
		end
	end

	ThreeQ.qbsolv!(model; timeout=60)

    kbs = q.value
    #dwsolver = ThreeQ.DWQMI.getvfyc(mytoken)

    #x = ThreeQ.solvesapi!(Q; token=mytoken, auto_scale=true, solver=dwsolver, num_reads=num_reads, finishsolve_helper=finishsolve_helper, reuse_embedding=true, param_chain=maximum(Q))

    return kbs
end

# get solutions from the D-Wave 
# input: forward modeling operator, qubo matrix, h measurements, initial guess of permeability, low and high permeability, number of times i grid is within j grid, location of measured head values, D-Wave token
# output: all valid permutations, the linear and non-linear objective function value of that permutation, the indices of permutations that minimize the non-linear objective function
function getPermutationsDwaveLinear(forwardModel, Q, hhat, initGuess, kl, kh, inum, lvals, mytoken; flip=true, param_chain=0)
    ##################################################
    ### CHANGE BETWEEN QBSOLV AND DWAVE HERE ###
    @printf("Getting D-Wave solutions...\n")
    kbs = getDwaveAnswer(Q, mytoken, param_chain=param_chain)  
    #kbs = getDwaveAnswerqbsolv(Q, mytoken)      # works without internet
    @printf("Done.\n\n")
    ##################################################


    #kbs2 = getDwaveAnswerqbsolv(Q, mytoken)  
    kbs = unique(kbs, dims=2)      
    norm = zeros(size(kbs,2))
    #nlObj = zeros(size(kbs,2))

    #norm2 = zeros(size(kbs2,2))
    #nlObj2 = zeros(size(kbs2,2))

    for i = 1:size(kbs,2)
        norm[i] = kbs[:,i]'*Q*kbs[:,i]
        #@show sum(kbs[:,i])

        ## get non-linear objective function value
        #if flip
        #    h0 = geth0(forwardModel, flipBits(initGuess, kbs[:,i]), kl, kh, inum, lvals)
        #else
        #    h0 = geth0(forwardModel, kbs[:,i], kl, kh, inum, lvals)
        #end
        #nlObj[i] = checkObj(h0, hhat)
    end

    #for i = 1:size(kbs2,2)
    #    norm2[i] = kbs2[:,i]'*Q*kbs2[:,i]
    #    #@show sum(kbs[:,i])

    #    # get non-linear objective function value
    #    if flip
    #        h02 = geth0(forwardModel, flipBits(initGuess, kbs2[:,i]), kl, kh, inum, lvals)
    #    else
    #        h02 = geth0(forwardModel, kbs2[:,i], kl, kh, inum, lvals)
    #    end
    #    nlObj2[i] = checkObj(h02, hhat)
    #end

    #minVal = minimum(nlObj)
    #idxMinVal = findall(x->x==minVal, nlObj)[1]
    minVallin = minimum(norm)
    idxMinVallin = findall(x->x==minVallin, norm)[1]
    #@printf("qQq non-linear = %f\n", kbs[:,idxMinVal]'*Q*kbs[:,idxMinVal])
    #@printf("qQq linear = %f\n", minVallin)

    #minVal2 = minimum(nlObj2)
    #idxMinVal2 = findall(x->x==minVal2, nlObj2)[1]
    #minVallin2 = minimum(norm2)
    #idxMinVallin2 = findall(x->x==minVallin2, norm2)[1]

    #@printf("qQq qbsolv = %f\n", minVallin2)
    
    #if idxMinVal != idxMinVallin
    #    @printf("WARNING: non-linear and linear minimum NOT the same\n")
    #    @show idxMinVal, idxMinVallin
    #    @show kbs[:,idxMinVal] - kbs[:,idxMinVallin]
    #end
    #@printf("nlmin w/ corrections: %f; nlmin w/o corrections: %f\n", minVal2, minVal)
    return Array(Float64.(kbs)), norm, idxMinVallin
end

# get solutions from the D-Wave 
# input: forward modeling operator, qubo matrix, h measurements, initial guess of permeability, low and high permeability, number of times i grid is within j grid, location of measured head values, D-Wave token
# output: all valid permutations, the linear and non-linear objective function value of that permutation, the indices of permutations that minimize the non-linear objective function
function getPermutationsDwaveLinearnl(forwardModel, Q, hhat, initGuess, kl, kh, inum, lvals, mytoken; flip=true, param_chain=0)
    ##################################################
    ### CHANGE BETWEEN QBSOLV AND DWAVE HERE ###
    @printf("Getting D-Wave solutions...\n")
    kbs = getDwaveAnswer(Q, mytoken, param_chain=param_chain)  
    #kbs = getDwaveAnswerqbsolv(Q, mytoken)      # works without internet
    @printf("Done.\n\n")
    ##################################################


    #kbs2 = getDwaveAnswerqbsolv(Q, mytoken)  
    kbs = unique(kbs, dims=2)      
    norm = zeros(size(kbs,2))

    #norm2 = zeros(size(kbs2,2))
    #nlObj2 = zeros(size(kbs2,2))
    #

    for i = 1:size(kbs,2)
        norm[i] = kbs[:,i]'*Q*kbs[:,i]
        #@show sum(kbs[:,i])

        ## get non-linear objective function value
        #if flip
        #    h0 = geth0(forwardModel, flipBits(initGuess, kbs[:,i]), kl, kh, inum, lvals)
        #else
        #    h0 = geth0(forwardModel, kbs[:,i], kl, kh, inum, lvals)
        #end
        #nlObj[i] = checkObj(h0, hhat)
    end

    # get minimum numnl solutions
    rt = copy(norm)
    #numnl = 8
    #numnl = 24
    numnl = 48
    nlvals = floor.(Int, zeros(numnl))
    maxval = maximum(rt)
    for i = 1:numnl
        minidx = argmin(rt)
        nlvals[i] = minidx
        rt[minidx] = maxval
    end

    # get min nl val
    #nlObj = zeros(numnl)
    #for i = 1:numnl #size(kbs2,2)
    nlObj = SharedArray{Float64}(numnl)
    @sync @distributed for i = 1:numnl
        # get non-linear objective function value
        if flip
            h02 = geth0(forwardModel, flipBits(initGuess, kbs[:,nlvals[i]]), kl, kh, inum, lvals)
        else
            h02 = geth0(forwardModel, kbs[:,nlvals[i]], kl, kh, inum, lvals)
        end
        nlObj[i] = checkObj(h02, hhat)
    end

    minVal = minimum(nlObj)
    idxMinVal = nlvals[findall(x->x==minVal, nlObj)[1]]

    minVallin = minimum(norm)
    idxMinVallin = findall(x->x==minVallin, norm)[1]

    #@printf("qQq non-linear = %f\n", kbs[:,idxMinVal]'*Q*kbs[:,idxMinVal])
    #@printf("qQq linear = %f\n", minVallin)

    #minVal2 = minimum(nlObj2)
    #idxMinVal2 = findall(x->x==minVal2, nlObj2)[1]
    #minVallin2 = minimum(norm2)
    #idxMinVallin2 = findall(x->x==minVallin2, norm2)[1]

    #@printf("qQq qbsolv = %f\n", minVallin2)
    
    #if idxMinVal != idxMinVallin
    #    @printf("WARNING: non-linear and linear minimum NOT the same\n")
    #    @show idxMinVal, idxMinVallin
    #    @show kbs[:,idxMinVal] - kbs[:,idxMinVallin]
    #end
    #@printf("nlmin w/ corrections: %f; nlmin w/o corrections: %f\n", minVal2, minVal)
    return Array(Float64.(kbs))[:,idxMinVal], norm[idxMinVal], minVal
    #return Array(Float64.(kbs)), norm, idxMinVallin
end

# get solutions from the D-Wave 
# input: forward modeling operator, qubo matrix, h measurements, initial guess of permeability, low and high permeability, number of times i grid is within j grid, location of measured head values, D-Wave token
# output: all valid permutations, the linear and non-linear objective function value of that permutation, the indices of permutations that minimize the non-linear objective function
function getPermutationsDwave(forwardModel, Q, hhat, initGuess, kl, kh, inum, lvals, mytoken; flip=true, param_chain=0)
    ##################################################
    ### CHANGE BETWEEN QBSOLV AND DWAVE HERE ###
    kbs = getDwaveAnswer(Q, mytoken, param_chain=param_chain)  
    #kbs = getDwaveAnswerqbsolv(Q, mytoken)      # works without internet
    ##################################################


    #kbs2 = getDwaveAnswerqbsolv(Q, mytoken)  
    kbs = unique(kbs, dims=2)      
    norm = zeros(size(kbs,2))
    nlObj = zeros(size(kbs,2))

    #norm2 = zeros(size(kbs2,2))
    #nlObj2 = zeros(size(kbs2,2))

    for i = 1:size(kbs,2)
        norm[i] = kbs[:,i]'*Q*kbs[:,i]
        #@show sum(kbs[:,i])

        # get non-linear objective function value
        if flip
            h0 = geth0(forwardModel, flipBits(initGuess, kbs[:,i]), kl, kh, inum, lvals)
        else
            h0 = geth0(forwardModel, kbs[:,i], kl, kh, inum, lvals)
        end
        nlObj[i] = checkObj(h0, hhat)
    end

    #for i = 1:size(kbs2,2)
    #    norm2[i] = kbs2[:,i]'*Q*kbs2[:,i]
    #    #@show sum(kbs[:,i])

    #    # get non-linear objective function value
    #    if flip
    #        h02 = geth0(forwardModel, flipBits(initGuess, kbs2[:,i]), kl, kh, inum, lvals)
    #    else
    #        h02 = geth0(forwardModel, kbs2[:,i], kl, kh, inum, lvals)
    #    end
    #    nlObj2[i] = checkObj(h02, hhat)
    #end

    minVal = minimum(nlObj)
    idxMinVal = findall(x->x==minVal, nlObj)[1]
    minVallin = minimum(norm)
    idxMinVallin = findall(x->x==minVallin, norm)[1]
    #@printf("qQq non-linear = %f\n", kbs[:,idxMinVal]'*Q*kbs[:,idxMinVal])
    #@printf("qQq linear = %f\n", minVallin)

    #minVal2 = minimum(nlObj2)
    #idxMinVal2 = findall(x->x==minVal2, nlObj2)[1]
    #minVallin2 = minimum(norm2)
    #idxMinVallin2 = findall(x->x==minVallin2, norm2)[1]

    #@printf("qQq qbsolv = %f\n", minVallin2)
    
    if idxMinVal != idxMinVallin
        @printf("WARNING: non-linear and linear minimum NOT the same\n")
        @show idxMinVal, idxMinVallin
        @show kbs[:,idxMinVal] - kbs[:,idxMinVallin]
    end
    #@printf("nlmin w/ corrections: %f; nlmin w/o corrections: %f\n", minVal2, minVal)
    return Array(Float64.(kbs)), norm, nlObj,idxMinVal
end

# get one iteration of the model update
# inputs: forward modeling operator, initial guess of solution, k_low, k_high, number of times i grid is within j grid, locations where sensors are, measured values of head, iteration number
# output: updated solution, qubo matrix
function oneIter(forwardModel, initGuess, kl, kh, inum, lvals, hhat, iterNum)
    h0 = geth0(forwardModel, initGuess, kl, kh, inum, lvals)
    kjGuess = b2f(initGuess, kl, kh)    # convert to acutal perm 
    F = getF(forwardModel, lvals, kjGuess, kl, kh, inum)
    Q = qubo(hhat, F, h0)               # get Q matrix
    
    # check all possible solutions
    permsValue, objVal, binVal = getPermutations(Q)
    
    # find index and objective function value of true answer
    objSoln = kBin'*Q*kBin
    
    # get all q values that minimize
    q = permsValue[:, binVal]
    
    # display parameters
    @printf("ITERATION NUMBER %i",iterNum)
    nlObj = Array{Float64}(undef, size(q,2))
    for i=1:length(binVal)
        q = permsValue[:, binVal[i]]
        soln = flipBits(initGuess,q) 
        @printf("    h%i = ",i); printArray(soln)
    
        h0 = geth0(forwardModel, soln, kl, kh, inum, lvals)
        obj = checkObj(h0, hhat)
        nlObj[i] = obj
        lObj = soln'*Q*soln
        @printf("\n          # wrong = %i, lobj = %f, nobj = %.8E\n",sum(abs.(kBin.-soln)), lObj, obj)
    end
    idxSoln = findall(x->x==minimum(nlObj), nlObj)
    @printf("\n    Solution that minimizes non-linear objective function (%i total)\n", length(idxSoln))
    i = idxSoln[1]
    q = permsValue[:, binVal[i]]
    soln = flipBits(initGuess,q) 
    @printf("    h%i = ",i)
    printArray(soln)
    @printf("\n          # wrong = %i, ",sum(abs.(kBin .- soln)))
    
    h0 = geth0(forwardModel, soln, kl, kh, inum, lvals)
    obj = checkObj(h0, hhat)
    nlObj[i] = obj
    lObj = soln'*Q*soln
    @printf("lobj = %f",lObj)
    @printf(", nobj = %.8E\n",obj)
    @printf("\n\n")

    return soln, Q
end

# update the solution
# inputs: initial guess, q values
# output: updated guess
function flipBits(k0, q)
    flipped = copy(k0)
    for i=1:length(k0)
        if q[i] == 1 && k0[i]==1
            flipped[i] = 0
        elseif q[i] == 1 && k0[i]==0
            flipped[i] = 1
        end
    end
    return flipped
end

# get one iteration of the model update
# inputs: forward modeling operator, Q matrix, initial guess of solution, k_low, k_high, number of times i grid is within j grid, locations where sensors are, measured values of head, iteration number, binary permeability, D-wave token
# output: updated solution, linear and non-linear objective function values
function oneIterDwaveLinear(forwardModel, Q, initGuess, kl, kh, inum, lvals, hhat, iterNum, kBin, mytoken; flip=true, param_chain=0, retSolns=false)
    ##Q = getQFunction(initGuess, kl, kh, inum, lvals)
    #Q2 = getQ(initGuess, kl, kh, inum, lvals)
    #@show Q2 == Q

    # check all possible solutions
    #permsValue, objVal, binVal = getPermutations(Q)
    #
    permsValue, objVal, binVal = getPermutationsDwaveLinear(forwardModel, Q, hhat, initGuess, kl, kh, inum, lvals, mytoken, flip=flip, param_chain=param_chain)
    
    # get the first value...
    binVal = binVal[1]
    
    # display parameters
    # get all q values that minimize
    q = permsValue[:, binVal]
    objSoln = q'*Q*q
    
    #@printf("ITERATION NUMBER %i, lobj_true = %f\n",iterNum, objSoln[1])
    @printf("ITERATION NUMBER %i\n",iterNum)
    q = permsValue[:, binVal[1]]
    if flip
        soln = floor.(Int,flipBits(initGuess,q) )
    else
        soln = floor.(Int,q)
    end

    # get non-linear objective function
    if flip
        h0 = geth0(forwardModel, flipBits(initGuess, q), kl, kh, inum, lvals)
    else
        h0 = geth0(forwardModel, q, kl, kh, inum, lvals)
    end
    nlObj = checkObj(h0, hhat)

    @printf("     h%i = ",iterNum)
    printArray(soln)
    @printf("\n     # wrong = %i, ",sum(abs.(kBin .- soln)))
    
    @printf("lobj = %f",objVal[binVal])
    @printf(", nobj = %.8E\n",nlObj)
    @printf("\n\n")

    if retSolns
        if flip
            allsolns = 0 .*permsValue
            for i = 1:size(allsolns,2)
                allsolns[:,i] = HydroThreeQ.flipBits(permsValue[:,i],initGuess)
            end
            allsolns = floor.(Int, allsolns)
        else
            allsolns = permsValue
        end
        linobjvals = zeros(size(allsolns,2))
        for i = 1:size(allsolns,2)
            linobjvals[i] = allsolns[:,i]'*Q*allsolns[:,i]
        end
        return soln, objVal[binVal], nlObj, allsolns, linobjvals
    else
        return soln, objVal[binVal], nlObj
    end
end

# get one iteration of the model update
# inputs: forward modeling operator, Q matrix, initial guess of solution, k_low, k_high, number of times i grid is within j grid, locations where sensors are, measured values of head, iteration number, binary permeability, D-wave token
# output: updated solution, linear and non-linear objective function values
function oneIterDwaveLinearnl(forwardModel, Q, initGuess, kl, kh, inum, lvals, hhat, iterNum, kBin, mytoken; flip=true, param_chain=0, retSolns=false)
    ##Q = getQFunction(initGuess, kl, kh, inum, lvals)
    #Q2 = getQ(initGuess, kl, kh, inum, lvals)
    #@show Q2 == Q

    # check all possible solutions
    #permsValue, objVal, binVal = getPermutations(Q)
    #
    q, objSoln, nlObj = getPermutationsDwaveLinearnl(forwardModel, Q, hhat, initGuess, kl, kh, inum, lvals, mytoken, flip=flip, param_chain=param_chain)
    
    
    #@printf("ITERATION NUMBER %i, lobj_true = %f\n",iterNum, objSoln[1])
    @printf("ITERATION NUMBER %i\n",iterNum)
    if flip
        soln = floor.(Int,flipBits(initGuess,q) )
    else
        soln = floor.(Int,q)
    end

    ## get non-linear objective function
    #if flip
    #    h0 = geth0(forwardModel, flipBits(initGuess, q), kl, kh, inum, lvals)
    #else
    #    h0 = geth0(forwardModel, q, kl, kh, inum, lvals)
    #end
    #nlObj = checkObj(h0, hhat)

    @printf("     h%i = ",iterNum)
    printArray(soln)
    @printf("\n     # wrong = %i, ",sum(abs.(kBin .- soln)))
    
    @printf("lobj = %f",objSoln)
    @printf(", nobj = %.8E\n",nlObj)
    @printf("\n\n")

    if retSolns # did not fix
        if flip
            allsolns = 0 .*permsValue
            for i = 1:size(allsolns,2)
                allsolns[:,i] = HydroThreeQ.flipBits(permsValue[:,i],initGuess)
            end
            allsolns = floor.(Int, allsolns)
        else
            allsolns = permsValue
        end
        linobjvals = zeros(size(allsolns,2))
        for i = 1:size(allsolns,2)
            linobjvals[i] = allsolns[:,i]'*Q*allsolns[:,i]
        end
        return soln, objSoln, nlObj, allsolns, linobjvals
    else
        return soln, objSoln, nlObj
    end
end

# get one iteration of the model update
# inputs: forward modeling operator, Q matrix, initial guess of solution, k_low, k_high, number of times i grid is within j grid, locations where sensors are, measured values of head, iteration number, binary permeability, D-wave token
# output: updated solution, linear and non-linear objective function values
function oneIterDwave(forwardModel, Q, initGuess, kl, kh, inum, lvals, hhat, iterNum, kBin, mytoken; flip=true, param_chain=0, retSolns=false)
    ##Q = getQFunction(initGuess, kl, kh, inum, lvals)
    #Q2 = getQ(initGuess, kl, kh, inum, lvals)
    #@show Q2 == Q

    # check all possible solutions
    #permsValue, objVal, binVal = getPermutations(Q)
    #
    permsValue, objVal, nlObjVal, binVal = getPermutationsDwave(forwardModel, Q, hhat, initGuess, kl, kh, inum, lvals, mytoken, flip=flip, param_chain=param_chain)
    
    # get the first value...
    binVal = binVal[1]
    
    # display parameters
    # get all q values that minimize
    q = permsValue[:, binVal]
    objSoln = q'*Q*q
    
    #@printf("ITERATION NUMBER %i, lobj_true = %f\n",iterNum, objSoln[1])
    @printf("ITERATION NUMBER %i\n",iterNum)
    q = permsValue[:, binVal[1]]
    if flip
        soln = floor.(Int,flipBits(initGuess,q) )
    else
        soln = floor.(Int,q)
    end

    @printf("     h%i = ",iterNum)
    printArray(soln)
    @printf("\n     # wrong = %i, ",sum(abs.(kBin .- soln)))
    
    @printf("lobj = %f",objVal[binVal])
    @printf(", nobj = %.8E\n",nlObjVal[binVal])
    @printf("\n\n")

    if retSolns
        if flip
            allsolns = 0 .*permsValue
            for i = 1:size(allsolns,2)
                allsolns[:,i] = HydroThreeQ.flipBits(permsValue[:,i],initGuess)
            end
            allsolns = floor.(Int, allsolns)
        else
            allsolns = permsValue
        end
        linobjvals = zeros(size(allsolns,2))
        for i = 1:size(allsolns,2)
            linobjvals[i] = allsolns[:,i]'*Q*allsolns[:,i]
        end
        return soln, objVal[binVal], nlObjVal[binVal], allsolns, linobjvals
    else
        return soln, objVal[binVal], nlObjVal[binVal]
    end
end

function oneIterDwaveNH(forwardModel, initGuess; mytoken=Main.mytoken, epsilon=0.0, groundstate=false)
    forwardModel0 = forwardModel(initGuess)
    F = getF(forwardModel, 1:length(forwardModel0), initGuess, 0, 1, 1)
    Q = qubo(zeros(length(forwardModel0)), F, forwardModel0)
    Q = Q + epsilon.*LinearAlgebra.I
    
    # check all possible solutions
    permsValue, objVal, nlObjVal, binVal = HydroThreeQ.getPermutationsDwave(forwardModel, Q, zeros(length(forwardModel0)), initGuess, 0, 1, 1, 1:length(forwardModel0), mytoken)
   
    # get the first value...
    binVal = binVal[1]
  
    # display parameters
    # get all q values that minimize
    q = permsValue[:, binVal]
    objSoln = q'*Q*q
 
    q = permsValue[:, binVal[1]]
    soln = HydroThreeQ.flipBits(initGuess,q)
    HydroThreeQ.printArray(soln)
    if groundstate != false
        @printf("\n     # wrong = %i, ",sum(abs.(groundstate .- soln)))
    end
    
    @printf("lobj = %f",objVal[binVal])
    @printf(", nobj = %.8E\n",nlObjVal[binVal])
    @printf("\n\n")
    
    return soln, objVal[binVal], nlObjVal[binVal]
end


end
