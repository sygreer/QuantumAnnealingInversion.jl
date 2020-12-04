using PyPlot
#for i = 1:size(savenwrong,2)
#    scatter([0:20], savenwrong[:,i])
#end
#
#xlim([0,21])
#ylim([0, 26])

figure(figsize=[5,5])
x = Array(0:20)
for i = 1:size(savenwrong,2)
    plot(x, savenwrong[:,i], alpha=0.5, linewidth=1)
end

endptsx = zeros(size(savenwrong,2))
endptsy = zeros(size(savenwrong,2))
for i = 1:size(savenwrong,2)
    stoploc = findfirst(isnan, savenwrong[:,i])
    if stoploc isa Nothing
        #scatter(10, y[i,10]-1, alpha=0.5)
        endptsx[i] = 20
        endptsy[i] = savenwrong[20,i]-1
    else
        stoploc = stoploc-1
        yvalthere = savenwrong[stoploc,i]
        #scatter(stoploc-1, yvalthere, alpha=0.5)
        endptsx[i] = stoploc-1
        endptsy[i] = yvalthere
    end
end
for i = 1:size(savenwrong,2)
    scatter(endptsx[i], endptsy[i], alpha=0.2)
end

## savenwronglarge
#text(11,0,"2",  horizontalalignment="center", verticalalignment="center")
#text(8,0,"3",  horizontalalignment="center", verticalalignment="center")
#text(6,7,"1",  horizontalalignment="center", verticalalignment="center")
#text(11,3,"3",  horizontalalignment="center", verticalalignment="center")
#text(12,1,"1",  horizontalalignment="center", verticalalignment="center")

## savenwronglargehc
#text(8,4,"1",  horizontalalignment="center", verticalalignment="center")
#text(11,6,"1",  horizontalalignment="center", verticalalignment="center")
#text(6,10,"2",  horizontalalignment="center", verticalalignment="center")
#text(12,7,"1",  horizontalalignment="center", verticalalignment="center")
#text(12,8,"1",  horizontalalignment="center", verticalalignment="center")
#text(9,4,"1",  horizontalalignment="center", verticalalignment="center")
#text(8,5,"2",  horizontalalignment="center", verticalalignment="center")
#text(9,5,"1",  horizontalalignment="center", verticalalignment="center")

#text(9,0,"10",  horizontalalignment="center", verticalalignment="center")
text(2,17,"5",  horizontalalignment="center", verticalalignment="center")
text(4,17,"2",  horizontalalignment="center", verticalalignment="center")
text(4,15,"2",  horizontalalignment="center", verticalalignment="center")
text(7,15,"1",  horizontalalignment="center", verticalalignment="center")

xlim([-0.3, 12.3])
ylim([-0.3, 28])
#title("Convergence pattern, low contrast")
title("Convergence pattern, D-Wave ")
xlabel("Number of iterations (k)")
ylabel("Number of incorrect permeability values")
yticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28])
xticks([0, 2, 4, 6, 8, 10, 12])

show()

#savefig("conv.pdf", bbox_inches="tight")
