"""
    setCenter!(obj::DenseInfiniteMPS{L}, si::Int64;
        check::Bool=false, tol::Float64=Defaults.tol,
        alg::EigenAlgorithm=Defaults.alg_eig, kwargs...)

Set `obj.Center` to `si` according to subsidiary conditions.

# Arguments
`obj::DenseInfiniteMPS{L}`: iMPS or iMPO with cell length `L`
`si::Int64`: expected canonical center

# Keyword Arguments
`check::Bool=false`: whether check local tensors on the left side of `si` are all left-canonical and the ones on the right side are all right-canonical
`tol::Bool=Defaults.tol`: propagated to `isLeftIsometric` and `isRightIsometric`
`alg::EigenAlgorithm`: propagated to `canonicalize!`
`kwargs`: propagated to `leftorth` or `rightorth`
"""
function setCenter!(obj::DenseInfiniteMPS{L}, si::Int64;
            check::Bool=false, tol::Float64=Defaults.tol,
            alg::EigenAlgorithm=Defaults.alg_eig, kwargs...) where L

    if isnothing(obj.Center) || (check && !(all(x->isLeftIsometric(x, false; tol=tol), obj.A[1:si-1]) && all(x->isRightIsometric(x, false; tol=tol), obj.A[si+1:end])))
        c = obj.c
        canonicalize!(obj, si; alg=alg)
        obj.c = c
    else
        if si == obj.Center
            nothing
        elseif si < obj.Center
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
        return obj
    end
end
