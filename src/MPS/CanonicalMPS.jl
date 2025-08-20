"""
    mutable struct CanonicalMPS{L, T<:Union{Float64, ComplexF64}} <: DenseInfiniteMPS{L, T}
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

# Notes
By convention, we appoint that:
- `AL[i] * C[i]` = `AC[i]` = `C[i-1] * AR[i]`
- `AL[i]' * AL[i] = I`
- `AR[i] * AR[i]' = I`
The first virtual space is between AL[1] and AL[2], and the one behind AL[L] same as the one ahead AL[1] is the last one.

# Constructors
    CanonicalMPS{L, T}(AL::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                       AR::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                       AC::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                       C::AbstractVector{<:AbstractBondTensor}=Vector{BondTensor}(undef, L))
    CanonicalMPS(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype)
    CanonicalMPS(AL::AbstractVector{<:AbstractMPSTensor}, C::AbstractVector{<:AbstractBondTensor})
    CanonicalMPS(obj::UniformMPS)
"""
mutable struct CanonicalMPS{L, T<:Union{Float64, ComplexF64}} <: DenseInfiniteMPS{L, T}
    const AL::AbstractVector{<:AbstractMPSTensor}
    const AR::AbstractVector{<:AbstractMPSTensor}
    const AC::AbstractVector{<:AbstractMPSTensor}
    const C::AbstractVector{<:AbstractBondTensor}
end
const CMPS = CanonicalMPS

"""
    randCMPS([::Type{T},] pspace::Vector{VectorSpace}, aspace::Vector{VectorSpace}) -> CMPS{L}

Generate a length `L` random CMPS with given length `L` vector `pspace` and `aspace`. `T = ComplexF64`(default) or `Float64` is the number type.

    randCMPS([::Type{T},] L::Int64, pspace::VectorSpace, apsace::VectorSpace) -> CMPS{L}

Assume the same `pspace` and `aspace`, except for the boundary bond, which is assumed to be trivial.

    randCMPS([::Type{T},] pdim::Vector{Int64}, adim::Vector{Int64}) -> CMPS{L}

Assume `pspace = [spacetype] .^ pdim` and `aspace = [spacetype] .^ adim`, where `spacetype = ℂ(T=ComplexF64) or ℝ(T=Float64)`

    randCMPS([::Type{T},] L::Int64, D::Int64, d::Int64) -> CMPS{L}

Assume the same `pdim` and `adim`.
"""
function randCMPS(::Type{T}, pspace::Vector{VectorSpace}, aspace::Vector{VectorSpace}) where T<:Union{Float64, ComplexF64}
end
function randCMPS(::Type{T}, L::Int64, pspace::VectorSpace, apsace::VectorSpace) where T<:Union{Float64, ComplexF64}
end
function randCMPS(::Type{T}, pdim::Vector{Int64}, adim::Vector{Int64}) where T<:Union{Float64, ComplexF64}
end
function randCMPS(::Type{T}, L::Int64, pdim::Int64, adim::Int64) where T<:Union{Float64, ComplexF64}
end
