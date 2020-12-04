# helper functions

# print an array with @printf
# input: array to print
# output: none (printed array)
function printArray(arr)
    if typeof(arr[1]) == Int64
        for i in arr
            @printf("%d ",i)
        end
    else
        for i in arr
            @printf("%f ",i)
        end
    end
end

# print an array with @printf
# input: array to print
# output: none (printed array)
function printMat(mat)
    I = size(mat)[1]
    J = size(mat)[2]
    @printf("\n")
    if typeof(mat[1]) == Int64
        for i = 1:I
            for j = 1:J
                @printf("%d ",mat[i,j])
            end
            @printf("\n")
        end
    else
        for i = 1:I
            for j = 1:J
                @printf("%f ",mat[i,j])
            end
            @printf("\n")
        end
    end
end


# checks the full objective function
# input: forward modeled head at all l, measured head at all l
# output: value of objective function for the input model
function checkObj(h0, hhat)
    sumval = 0
    for i=1:length(hhat)
        sumval += (h0[i] - hhat[i])^2
    end
    return sumval
end

