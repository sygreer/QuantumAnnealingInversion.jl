import HydroThreeQ
import PyPlot

solni = HydroThreeQ.kinn(HydroThreeQ.b2f(savesoln, kl, kh), inum)

fig, ax = PyPlot.subplots()

#scale = 3
#PyPlot.figure(figsize=(scale*length(ki)/inum,2*scale))

sizeA = length(ki)/inum
#scale2=175
scale2=100
PyPlot.figure(figsize=(scale2,scale2/sizeA*2))

ax = PyPlot.axes(frameon=false)
ax.get_xaxis().set_visible(false)
ax.get_yaxis().set_visible(false)

PyPlot.imshow([ki solni]', cmap="summer", aspect="auto" )
x = range(0,16,length=length(ki))
y = ones(length(lvals)).*0.5
PyPlot.plot(0.5.*ones(length(ki)), linewidth=10, c="k")
PyPlot.scatter(lvals,y, s=3000, marker="X", c="k")

PyPlot.annotate("Truth", xy=[1;0], fontsize=50)
PyPlot.annotate(@sprintf("Solution, %.1f%% correct", percentc), xy=[1;1], fontsize=50)
PyPlot.savefig("model.pdf", bbox_inches="tight")

PyPlot.close("all")
