ioff()

savenwrong=zeros(21,10)
for iterv = 1:10
    @show iterv
    @time include("seismic.jl")
    savenwrong[:,iterv] = nwrongarr
end
