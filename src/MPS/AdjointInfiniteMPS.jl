"""
    struct AdjointInfiniteMPS{L, T} <: AbstractInfiniteMPS{L}
        parent::Union{DenseUniformMPS{L, T}, DenseCanonicalMPS{L, T}}
    end

Lazy wrapper type for adjoint of iMPS.

    adjoint(::DenseInfiniteMPS) -> ::AdjointInfiniteMPS
    adjoint(::AdjointInfiniteMPS) -> ::DenseInfiniteMPS

Functions to be directly propagated to the parent:

    normalize!

Functions to be propagated to the parent with some adaptations:

    getAL, getAR, getAC, getC, leftVirtualSpace, rightVirtualSpace, physicalSpace, extraPhysicalSpace
"""
struct AdjointInfiniteMPS{L, T} <: AbstractInfiniteMPS{L}
    parent::DenseInfiniteMPS
end
adjoint(A::DenseInfiniteMPS) = AdjointInfiniteMPS(A)
adjoint(A::AdjointInfiniteMPS) = A.parent

function show(io::IO, obj::AdjointInfiniteMPS)
     print(io, "Adjoint of ")
     show(io, obj.parent)
end

# some functions to be directly propagated to the parent
normalize!(obj::AdjointInfiniteMPS) = normalize!(obj.parent)

# apply lazy adjoint when obtaining the local tensors
for func in (:getAL, :getAR, :getAC, :getC, :leftVirtualSpace, :rightVirtualSpace, :physicalSpace, :extraPhysicalSpace)
    @eval $func(obj::AdjointInfiniteMPS, args...) = adjoint.($func(obj.parent, args...))
end

for func in (:EE,)
    @eval $func(obj::AdjointInfiniteMPS, args...) = ($func(obj.parent, args...))
end
