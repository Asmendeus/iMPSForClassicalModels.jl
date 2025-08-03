"""
     abstract type AbstractInfiniteMPO{W, L}

Abstract type of all iMPO with length `L` and width `W` of cell.
"""
abstract type AbstractInfiniteMPO{W, L} end

length(::AbstractInfiniteMPO{W, L}) where {W, L} = W * L
size(::AbstractInfiniteMPO{W, L}) where {W, L} = (W, L)
for func in (:getindex, :lastindex, :setindex!, :iterate, :keys, :isassigned)
     @eval Base.$func(obj::AbstractInfiniteMPO, args...) = $func(obj.A, args...)
end