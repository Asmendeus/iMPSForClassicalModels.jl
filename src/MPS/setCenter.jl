"""
    setCenter!(obj::DenseInfiniteMPS{L, T}, si::Int64; kwargs...)

Set `Center` of an iMPS/iMPO with canonical form to `si`.

# Arguments
`obj::DenseInfiniteMPS{L}`: iMPS or iMPO with cell length `L`
`si::Int64`: expected canonical center

# Keyword Arguments
`kwargs`: propagated to `leftorth` or `rightorth`
"""
function setCenter!(obj::DenseInfiniteMPS{L, T}, si::Int64; kwargs...) where {L, T}

    iscanonical(obj) || throw(ArgumentError("Illegal behavior of setting center for an iMPS/iMPO with uniform form"))

    if si == obj.Center
        nothing
    elseif si > obj.Center
        for l in obj.Center:si-1
            obj.A[l], X, _ = leftorth(obj.A[l]; kwargs...)
            obj.A[l+1] = X * obj.A[l+1]
        end
    else
        for l in obj.Center:-1:si+1
            X, obj.A[l], _ = rightorth(obj.A[l]; kwargs...)
            obj.A[l-1] = obj.A[l-1] * X
        end
    end
    obj.Center = si
    return obj
end

"""
    setCenter(obj::DenseInfiniteMPS{L}, si::Int64; kwargs...)

Pure function version of `setCenter!`
"""
function setCenter(obj::DenseInfiniteMPS{L, T}, si::Int64; kwargs...) where {L, T}
    obj′ = deepcopy(obj)
    return setCenter!(obj′, si; kwargs...)
end
