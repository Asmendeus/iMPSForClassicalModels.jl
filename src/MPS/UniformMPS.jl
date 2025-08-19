"""
    mutable struct UniformMPS{L, T<:Union{FLoat64, ComplexF64}} <: DenseInfiniteMPS{L, T}
        const A::AbstractVector{<:AbstractMPSTensor}
    end

Concrete type of iMPS with uniform form, where `L` is the cell length, `T == Float64` or `ComplexF64` is the number type of local tensors.

# Fields
    const A::AbstractVector{<:AbstractMPSTensor}
Length `L` vector to store the local tensors. Note the vector `A` is immutable while the local tensors in it are mutable.

# Constructors
    UniformMPS{L, T}(A::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L))
    UniformMPS(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype)
"""
mutable struct UniformMPS{L, T<:Union{FLoat64, ComplexF64}} <: DenseInfiniteMPS{L, T}
    const A::AbstractVector{<:AbstractMPSTensor}
end
const UMPS = UniformMPS
