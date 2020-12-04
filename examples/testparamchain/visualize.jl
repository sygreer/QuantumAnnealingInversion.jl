import HydroThreeQ
import PyPlot

PyPlot.ioff()

interpkBin = HydroThreeQ.twodinterp(arr2mat(kBin), intval)
interpsavesoln = HydroThreeQ.twodinterp(arr2mat(savesoln),intval)


PyPlot.close("all")
PyPlot.figure(figsize=(14,10))

PyPlot.subplot(2, 4, 1)
PyPlot.imshow(interpkBin, cmap="summer")
PyPlot.scatter([obsI.-1], [obsJ.-1], s=150, marker="X", c="r")
PyPlot.title("truth")
PyPlot.xticks([])
PyPlot.yticks([])

PyPlot.subplot(2, 4, 2)
PyPlot.imshow(interpsavesoln, cmap="summer")
PyPlot.scatter([obsI.-1], [obsJ.-1], s=150, marker="X", c="r")
PyPlot.title("soln")
PyPlot.xticks([])
PyPlot.yticks([])

PyPlot.subplot(2, 4, 3)
PyPlot.imshow(interpsavesoln.-interpkBin, cmap="PRGn", vmin=-1, vmax=1)
PyPlot.scatter([obsI.-1], [obsJ.-1], s=150, marker="X", c="r")
PyPlot.title("difference")
PyPlot.xticks([])
PyPlot.yticks([])

PyPlot.subplot(2, 4, 4)
PyPlot.imshow(Q)
PyPlot.title("Q")
PyPlot.xticks([])

PyPlot.subplot(2, 4, (5,8))
PyPlot.xticks([])
PyPlot.yticks([])
PyPlot.text(0.5, 0.85, @sprintf("# sensors = %i",numobs ), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.80, @sprintf("# qubits = %i",length(kBin) ), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.75, @sprintf("# possible solutions = %i",2^length(kBin) ), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.70, @sprintf("# iters = %i",iternum), fontsize=12, horizontalalignment="center")

PyPlot.text(0.5, 0.60, @sprintf("non-linear objective function = %f",nobjB ), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.55, @sprintf("percent correct = %.1f%%",percentc ), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.50, @sprintf("number wrong = %i",nwrong ), fontsize=12, horizontalalignment="center")

PyPlot.text(0.5, 0.40, @sprintf("random seed 1 = %i",randseed1), fontsize=12, horizontalalignment="center")
PyPlot.text(0.5, 0.35, @sprintf("random seed 2 = %i",randseed2), fontsize=12, horizontalalignment="center")

if @isdefined trynum
    @printf("Saving seed %i experiment %i\n", randseed1, trynum)
    PyPlot.text(0.5, 0.93, @sprintf("experiment number = %i",trynum ), fontsize=15, horizontalalignment="center")
    #PyPlot.savefig(@sprintf("size%i-%i-%i.png",sqrtnump,randseed2, trynum), bbox_inches="tight")
    PyPlot.savefig(@sprintf("param_chain%i-%i-%i.png",sqrtnump, randseed2, param_chain), bbox_inches="tight")
end


#PyPlot.show()
