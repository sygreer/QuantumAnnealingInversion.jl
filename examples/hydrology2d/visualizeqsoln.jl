
solns, energies, idxminval = HydroThreeQ.getPermutations(Q)
st = sort(energies)

soln = floor.(Int, solns[:,idxminval])
#@show soln
#@show kBin
solnsint = floor.(Int, solns)

solnPos = [i for i in 1:size(solnsint,2) if solnsint[:,i] == kBin][1]
ensoln = energies[solnPos]

idxsoln = findall(x->x==ensoln, st)


endv = 30
PyPlot.close("all")
PyPlot.figure(figsize=(15,8))

PyPlot.subplot(2,2,1)
PyPlot.plot(st)
PyPlot.scatter(idxsoln, st[idxsoln], c="r")
PyPlot.xlabel("solutions")
PyPlot.ylabel("energies")
#PyPlot.title(@sprintf("%ix%i energies, findiff = %i, tikhonov wt = %f",sqrtnump, sqrtnump, intval, tw))
PyPlot.title(@sprintf("%ix%i energies, findiff = %i",sqrtnump, sqrtnump, intval))

PyPlot.subplot(2,2,2)
PyPlot.scatter([1:endv],st[1:endv])
#PyPlot.scatter(1:length(st),st)
PyPlot.scatter(idxsoln, st[idxsoln], marker="X", c="r", s=150)
PyPlot.ylim([st[1]-0.0001, st[endv]+0.0001])
PyPlot.xlabel("solutions")
PyPlot.xlabel("solutions")
PyPlot.ylabel("energies")
#PyPlot.title(@sprintf("%ix%i energies, findiff = %i, tikhonov wt = %f",sqrtnump, sqrtnump, intval, tw))
PyPlot.title(@sprintf("%ix%i energies, findiff = %i",sqrtnump, sqrtnump, intval))

PyPlot.subplot(2,2,(3,4))
PyPlot.xticks([])
PyPlot.yticks([])

function strarr(arr)
    str = ""
    for i = 1:length(arr)
        str = string(str, arr[i] , " ")
    end
    return str
end


PyPlot.text(0.3, 0.85, @sprintf("%s","kBin = "), fontsize=12, horizontalalignment="center")
PyPlot.text(0.3, 0.78, @sprintf("%s","soln = "), fontsize=12, horizontalalignment="center" )
PyPlot.text(0.5, 0.85, @sprintf("%s",strarr(kBin)), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.78, @sprintf("%s",strarr(soln)), fontsize=12, horizontalalignment="center" )

PyPlot.text(0.5, 0.58, @sprintf("grid = %ix%i",sqrtnump,sqrtnump), fontsize=12, horizontalalignment="center" )
PyPlot.text(0.5, 0.51, @sprintf("findiff = %i",intval), fontsize=12, horizontalalignment="center" )

PyPlot.text(0.5, 0.38, @sprintf("seed1 = %i",randseed1), fontsize=12, horizontalalignment="center" )
PyPlot.text(0.5, 0.31, @sprintf("seed2 = %i",randseed2), fontsize=12, horizontalalignment="center" )

PyPlot.text(0.5, 0.21, @sprintf("k = %f, dk = %f",kl, dk), fontsize=12, horizontalalignment="center" )


#PyPlot.savefig(@sprintf("%i-%i-t.png", sqrtnump,intval))
PyPlot.savefig(@sprintf("%i-%i-%i.png", sqrtnump,intval,randseed1))
