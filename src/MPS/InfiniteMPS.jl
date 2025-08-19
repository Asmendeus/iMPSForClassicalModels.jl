"""
    mutable struct InfiniteMPS{L, T<:Union{Float64, ComplexF64}} <: DenseInfiniteMPS{L, T}
        const AL::AbstractVector{<:AbstractMPSTensor}
        const AR::AbstractVector{<:AbstractMPSTensor}
        const AC::AbstractVector{<:AbstractMPSTensor}
        const C::AbstractVector{<:AbstractBondTensor}
    end

Concrete type of iMPS with canonical form, where `L` is the cell length, `T == Float64` or `ComplexF64` is the number type of local tensors.

# Fields
    const AL::AbstractVector{<:AbstractMPSTensor}
Length `L` vector to store the left-canonical tensors.

    const AR::AbstractVector{<:AbstractMPSTensor}
Length `L` vector to store the right-canonical tensors.

    const AC::AbstractVector{<:AbstractMPSTensor}
Length `L` vector to store the center tensors.

    const C::AbstractVector{<:AbstractBondTensor}
Length `L` vector to store the center bond tensors.

# Constructors
"""
mutable struct InfiniteMPS{L, T<:Union{Float64, ComplexF64}} <: DenseInfiniteMPS{L, T}
    const AL::AbstractVector{<:AbstractMPSTensor}
    const AR::AbstractVector{<:AbstractMPSTensor}
    const AC::AbstractVector{<:AbstractMPSTensor}
    const C::AbstractVector{<:AbstractBondTensor}
end
const IMPS = InfiniteMPS
