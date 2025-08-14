"""
     abstract type AbstractInfiniteMPO{L}

Abstract type of all iMPO with length `L`, marking the size of cell.
"""
abstract type AbstractInfiniteMPO{L} end
const AbstractInfiniteMPSOrMPO{L} = Union{AbstractInfiniteMPS{L}, AbstractInfiniteMPO{L}}

length(::AbstractInfiniteMPO{L}) where L = L
for func in (:getindex, :lastindex, :setindex!, :iterate, :keys, :isassigned)
     @eval Base.$func(obj::AbstractInfiniteMPO, args...) = $func(obj.A, args...)
end

abstract type DenseInfiniteMPO{L} <: AbstractInfiniteMPO{L} end
const DenseInfiniteMPSOrMPO{L} = Union{DenseInfiniteMPS{L}, DenseInfiniteMPO{L}}
