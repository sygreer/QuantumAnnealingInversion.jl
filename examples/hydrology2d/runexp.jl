
for i = 4:4         # size of model space
    for k = 2:9     # random variables
    for j = 1:10    # number of experiments
        global dk, kl, kh, randseed1, randseed2, intval, numobs, sqrtnump, trynum
        sqrtnump = i
        trynum = j
        dk = 0.1     # difference between klow and khigh
        kl = 1.0        # klow
        kh = kl + dk    # khigh
        randseed1=k
        randseed2=k
        intval  = 5
        numobs = floor(Int,sqrtnump^2/2)
        include("hydrolin5.jl")
    end
end
