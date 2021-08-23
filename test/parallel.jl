module TestModule
using LinearAlgebra
using Distributed
using SharedArrays



export innerfunction
function innerfunction(x, a)
    return x*a
end

export batchfunction
function batchfunction(X, a)

    result = SharedArray{Float64}(size(X))

    @sync @distributed for i = 1:length(X)
        result[i] = innerfunction(X[i], a)
    end
    return result
end

end
