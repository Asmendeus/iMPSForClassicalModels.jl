"""
     abstract type AbstractInfiniteMPS{L}

Abstract type of all iMPS with length `L`, marking the size of cell.
"""
abstract type AbstractInfiniteMPS{L} end

length(::AbstractInfiniteMPS{L}) where L = L
for func in (:getindex, :lastindex, :setindex!, :iterate, :keys, :isassigned)
     @eval Base.$func(obj::AbstractInfiniteMPS, args...) = $func(obj.A, args...)
end

abstract type DenseInfiniteMPS{L} <: AbstractInfiniteMPS{L} end
