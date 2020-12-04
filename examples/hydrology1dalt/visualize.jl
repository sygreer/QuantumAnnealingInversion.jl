import HydroThreeQ
import PyPlot

solni = HydroThreeQ.kinn(HydroThreeQ.b2f(savesoln, kl, kh), inum)

#fig, ax = PyPlot.subplots()
#
##scale = 3
##PyPlot.figure(figsize=(scale*length(ki)/inum,2*scale))
#
#sizeA = length(ki)/inum
##scale2=175
#scale2=100
#PyPlot.figure(figsize=(scale2,scale2/sizeA))
#
#ax = PyPlot.axes(frameon=false)
#ax.get_xaxis().set_visible(false)
#ax.get_yaxis().set_visible(false)
#
##PyPlot.imshow([ki solni]', cmap="summer", aspect="auto" )
#PyPlot.imshow([ki ki]', cmap="summer", aspect="auto" )
##x = range(0,sizeModel,length=length(ki))
##PyPlot.plot(1.5.*ones(length(ki)), linewidth=10, c="k")
#
#y = ones(length(lvals)).*0.5
#PyPlot.scatter(lvals,y, s=5000, marker="X", c="r")
#
#PyPlot.annotate("Truth", xy=[1;0], fontsize=50)
#PyPlot.annotate(@sprintf("Solution, %.1f%% correct", percentc), xy=[1;1], fontsize=50)
#PyPlot.savefig("model.pdf", bbox_inches="tight", pad_inches=0)

PyPlot.close("all")


function plotSoln(soln, name)
    fig, ax = PyPlot.subplots()
    sizeA = length(soln)/inum

    scale2=100
    PyPlot.figure(figsize=(scale2,scale2/sizeA))
    
    ax = PyPlot.axes(frameon=false)
    ax.get_xaxis().set_visible(false)
    ax.get_yaxis().set_visible(false)
    
    PyPlot.imshow([soln soln]', cmap="summer", aspect="auto" )
    
    y = ones(length(lvals)).*0.5
    PyPlot.scatter(lvals,y, s=5000, marker="X", c="r")
    
    #PyPlot.annotate("Truth", xy=[1;0], fontsize=50)
    #PyPlot.annotate(@sprintf("Solution, %.1f%% correct", percentc), xy=[1;1], fontsize=50)
    PyPlot.savefig(@sprintf("%s.pdf",name), bbox_inches="tight", pad_inches=0)
    
    PyPlot.close("all")

end

plotSoln(ki, "model")
plotSoln(solni, "soln4")
