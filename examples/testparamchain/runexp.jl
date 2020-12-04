numtries = 75
#paramvals = floor.(Int, collect(range(5, length=20, stop=100)))
#paramvals = collect(5:5:100)
#paramvals = collect(60:5:100)
#paramvals = 10. .^(-3:6)
#paramvals = collect(0:0.1:2)
#paramvals = collect(2:0.1:4)
paramvals = collect(2.25:0.025:2.55)
storey=zeros(numtries, length(paramvals))
storex=zeros(numtries, length(paramvals))
for i = 4:4         # size of model space
    for j = 1:length(paramvals)    # number of experiments
        for k = 1:(numtries)    # random variables
            storex[k,j] = paramvals[j]
            global dk, kl, kh, randseed1, randseed2, intval, numobs, sqrtnump, trynum, param_chain
            param_chain = paramvals[j]
            trynum = 1
            sqrtnump = i
            dk = 0.1     # difference between klow and khigh
            kl = 1.0        # klow
            kh = kl + dk    # khigh
            randseed1=k + 80
            randseed2=k + 80
            intval  = 5
            numobs = floor(Int,sqrtnump^2/2)
            include("hydrolin5.jl")
            storey[k,j] = nwrong # number it got wrong
        end
    end
end

PyPlot.ion()
PyPlot.close("all")
fig, ax = PyPlot.subplots()
#PyPlot.xscale("log")
for i = 1:length(storex[:,1])
    ax.scatter(storex[i,:], storey[i,:], label=@sprintf("Seed = %i",i))
end
ax.scatter(storex[1,:], sum(storey, dims=1)./size(storey)[1], label="Average", marker="X", s=100)
#PyPlot.scatter(storex, storey)
ax.legend()
PyPlot.xlabel("param chain weight = maximum(sum.(Q))./#")
PyPlot.ylabel("number incorrect")
PyPlot.show()
