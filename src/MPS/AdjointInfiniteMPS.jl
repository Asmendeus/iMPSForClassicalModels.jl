"""
    struct AdjointInfiniteMPS{L} <: AbstractInfiniteMPS{L}
        parent::DenseInfiniteMPS{L}
    end

Lazy wrapper type for adjoint of iMPS.

    adjoint(::DenseInfiniteMPS) -> ::AdjointInfiniteMPS
    adjoint(::AdjointInfiniteMPS) -> ::DenseInfiniteMPS

Functions to be directly propagated to the parent:

    lastindex, length, keys, norm, normalize!, Center, iterate, canonicalize!

Functions to be propagated to the parent with some adaptations:

    getindex, setindex!, coef
"""
struct AdjointInfiniteMPS{L} <: AbstractInfiniteMPS{L}
    parent::DenseInfiniteMPS{L}
end
adjoint(A::DenseInfiniteMPS) = AdjointInfiniteMPS(A)
adjoint(A::AdjointInfiniteMPS) = A.parent

function show(io::IO, obj::AdjointInfiniteMPS)
     print(io, "Adjoint of ")
     show(io, obj.parent)
end

# apply lazy adjoint when obtaining the local tensors
getindex(obj::AdjointInfiniteMPS, inds...) = getindex(obj.parent, inds...)'
setindex!(obj::AdjointInfiniteMPS, X, inds...) = setindex!(obj.parent, X', inds...)
"""
    coef(obj::AdjointInfiniteMPS) = coef(obj.parent)'
"""
coef(obj::AdjointInfiniteMPS) = coef(obj.parent)'

# some functions to be directly propagated to the parent
for func in (:lastindex, :length, :keys,:norm, :normalize!, :Center, :iterate, :canonicalize!, :scalartype)
     @eval $func(obj::AdjointInfiniteMPS, args...) = $func(obj.parent, args...)
end

