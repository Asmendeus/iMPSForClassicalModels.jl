mutable struct InfiniteMPS{L, T <:Union{Float64, ComplexF64}} <: DenseMPS{L, T}
    const A::AbstractVector{AbstractMPSTensor}
    const Center::Vector{Int64}
    c::T
end
