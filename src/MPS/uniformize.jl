"""
    uniformize!(obj::DenseInfiniteMPS{L, T})

Uniformize an iMPS/iMPO with canonical form.
"""
function uniformize!(obj::DenseInfiniteMPS{L, T}) where {L, T}

    iscanonical(obj) || throw(ArgumentError("Illegal behavior of uniformizing an iMPS/iMPO with uniform form"))

    AL, _, _, C = getAllCanonicalFormTensors(obj)
    A = map(l -> inv(sqrt(C[l])) * AL[l] * sqrt(C[mod(l-1, L)+1]), 1:L)

    obj.A[:] = A
    obj.Center = nothing

    return obj
end

"""
    uniformize(obj::DenseInfiniteMPS{L, T})

Pure function version of `uniformize!`
"""
function uniformize(obj::DenseInfiniteMPS{L, T}) where {L, T}
    obj′ = deepcopy(obj)
    return uniformize!(obj′)
end