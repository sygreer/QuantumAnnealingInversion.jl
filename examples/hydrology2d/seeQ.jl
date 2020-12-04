PyPlot.subplot(1,2,1)
PyPlot.imshow(Qlarge)
PyPlot.title("High permeability contrast")
PyPlot.xticks([])
PyPlot.yticks([])
PyPlot.colorbar()

PyPlot.subplot(1,2,2)
PyPlot.imshow(Qsmall)
PyPlot.title("Low permeability contrast")
PyPlot.xticks([])
PyPlot.yticks([])
PyPlot.colorbar()

PyPlot.show()