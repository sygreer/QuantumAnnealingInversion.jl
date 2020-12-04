import PyPlot
solns, energies, idxminval = HydroThreeQ.getPermutations(Q)
st = sort(energies)

minsoln = floor.(Int, solns[:,idxminval])

solnsint = floor.(Int, solns)
solnPos = [i for i in 1:size(solnsint,2) if solnsint[:,i] == kBin][1]
ensoln = energies[solnPos]
idxsoln = findall(x->x==ensoln, st)

pickedsoln = floor.(Int, soln)
psolnPos = [i for i in 1:size(solnsint,2) if solnsint[:,i] == soln][1]
pensoln = energies[psolnPos]
pidxsoln = findall(x->x==pensoln, st)




endv = 30
#PyPlot.close("all")
PyPlot.figure(figsize=(15,15))

PyPlot.subplot(3,4,(1,2))
PyPlot.plot(st)
PyPlot.scatter(idxsoln, st[idxsoln], c="m")
PyPlot.scatter(pidxsoln, st[pidxsoln], c="c")
PyPlot.xlabel("solutions")
PyPlot.ylabel("energies")
#PyPlot.title(@sprintf("%ix%i energies, findiff = %i, tikhonov wt = %f",sqrtnump, sqrtnump, intval, tw))
PyPlot.title(@sprintf("Q energies"))

PyPlot.subplot(3,4,(3,4))
PyPlot.scatter([1:endv],st[1:endv])
#PyPlot.scatter(1:length(st),st)
PyPlot.scatter(idxsoln, st[idxsoln], marker="P", c="m", s=150)
PyPlot.scatter(pidxsoln, st[pidxsoln], marker="*", c="c", s=150)
PyPlot.ylim([st[1]-0.0001, st[endv]+0.0001])
PyPlot.xlim([0, endv+1])
PyPlot.xlabel("solutions")
PyPlot.xlabel("solutions")
PyPlot.ylabel("energies")
#PyPlot.title(@sprintf("%ix%i energies, findiff = %i, tikhonov wt = %f",sqrtnump, sqrtnump, intval, tw))
PyPlot.title(@sprintf("Q energies"))

PyPlot.subplot(3,4,(5,8))
PyPlot.xticks([])
PyPlot.yticks([])

function strarr(arr)
    str = ""
    for i = 1:length(arr)
        str = string(str, arr[i] , " ")
    end
    return str
end


PyPlot.text(0.5, 0.93, @sprintf("iter = %i",iternum), fontsize=15, horizontalalignment="center")
PyPlot.text(0.5, 0.88, @sprintf("k = %f, dk = %f",kl, dk), fontsize=12, horizontalalignment="center" )

PyPlot.text(0.5, 0.80, @sprintf("# sensors = %i",numobs ), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.75, @sprintf("# qubits = %i",length(kBin) ), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.70, @sprintf("# possible solutions = %i",2^length(kBin) ), fontsize=12, horizontalalignment="center")

PyPlot.text(0.5, 0.65, @sprintf("grid = %ix%i, findiff = %i",sqrtnump,sqrtnump, intval), fontsize=12, horizontalalignment="center" )

PyPlot.text(0.5, 0.60, @sprintf("random seed = %i",randseed2), fontsize=12, horizontalalignment="center" )

PyPlot.text(0.35, 0.50, @sprintf("%s","kBin = "), fontsize=12, horizontalalignment="right")
PyPlot.text(0.5, 0.50, @sprintf("%s",strarr(kBin)), fontsize=12, horizontalalignment="center")
PyPlot.text(0.35, 0.450, @sprintf("%s","minsoln = "), fontsize=12, horizontalalignment="right" )
PyPlot.text(0.5, 0.450, @sprintf("%s",strarr(minsoln)), fontsize=12, horizontalalignment="center" )
PyPlot.text(0.35, 0.40, @sprintf("%s","pickedsoln = "), fontsize=12, horizontalalignment="right" )
PyPlot.text(0.5, 0.40, @sprintf("%s",strarr(soln)), fontsize=12, horizontalalignment="center" )


sizeModel = length(kBin)
nwrong = sum(abs.(kBin.-soln))
percentc = (100*(sizeModel - nwrong)/sizeModel)

PyPlot.text(0.5, 0.3, @sprintf("non-linear objective function = %f",nobj ), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.25, @sprintf("percent correct = %.1f%%",percentc ), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.20, @sprintf("number wrong = %i", nwrong ), fontsize=12, horizontalalignment="center")




PyPlot.text(0.5, 0.10, @sprintf("mageneta + = true soln, cyan * = picked soln"), fontsize=12, horizontalalignment="center" )
PyPlot.text(0.5, 0.05, @sprintf("green = kl, yellow = kh, red x = source location"), fontsize=12, horizontalalignment="center" )







interpkBin = HydroThreeQ.twodinterp(arr2mat(kBin), intval)
interpsoln = HydroThreeQ.twodinterp(arr2mat(soln),intval)

PyPlot.subplot(3, 4, 9)
PyPlot.imshow(interpkBin, cmap="summer")
PyPlot.scatter([obsI.-1], [obsJ.-1], s=150, marker="X", c="r")
PyPlot.title("truth")
PyPlot.xlim([0-0.5, size(interpkBin)[1]-0.5])
PyPlot.ylim([0-0.5, size(interpkBin)[2]-0.5])
PyPlot.xticks([])
PyPlot.yticks([])

PyPlot.subplot(3, 4, 10)
PyPlot.imshow(interpsoln, cmap="summer")
PyPlot.scatter([obsI.-1], [obsJ.-1], s=150, marker="X", c="r")
PyPlot.title("soln")
PyPlot.xlim([0-0.5, size(interpsoln)[1]-0.5])
PyPlot.ylim([0-0.5, size(interpsoln)[2]-0.5])
PyPlot.xticks([])
PyPlot.yticks([])

PyPlot.subplot(3, 4, 11)
PyPlot.imshow(interpsoln.-interpkBin, cmap="PRGn", vmin=-1, vmax=1)
PyPlot.scatter([obsI.-1], [obsJ.-1], s=150, marker="X", c="r")
PyPlot.title("difference")
PyPlot.xlim([0-0.5, size(interpsoln)[1]-0.5])
PyPlot.ylim([0-0.5, size(interpsoln)[2]-0.5])
PyPlot.xticks([])
PyPlot.yticks([])

PyPlot.subplot(3, 4, 12)
PyPlot.imshow(Qsave)
PyPlot.title("Q")
PyPlot.xticks([])
PyPlot.yticks([])

#PyPlot.savefig(@sprintf("%i-%i-t.png", sqrtnump,intval))
#PyPlot.savefig(@sprintf("%i-%i-%i.png", sqrtnump,intval,randseed1))
