"""
    mutable struct UniformMPS{L, T<:Union{FLoat64, ComplexF64}} <: DenseInfiniteMPS{L, T}
        const A::AbstractVector{<:AbstractMPSTensor}
    end

Concrete type of iMPS with uniform form, where `L` is the cell length, `T == Float64` or `ComplexF64` is the number type of local tensors.

# Fields
    const A::AbstractVector{<:AbstractMPSTensor}
Length `L` vector to store the local tensors. Note the vector `A` is immutable while the local tensors in it are mutable.

# Notes
By convention, we appoint that the first virtual space is between A[1] and A[2], and the one behind A[L] same as the one ahead A[1] is the last one.

# Constructors
    UniformMPS{L, T}(A::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L))
    UniformMPS(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype)
    UniformMPS(obj::CanonicalMPS)
"""
mutable struct UniformMPS{L, T<:Union{FLoat64, ComplexF64}} <: DenseInfiniteMPS{L, T}
    const A::AbstractVector{<:AbstractMPSTensor}

    function UniformMPS{L, T}(A::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L)) where {L, T<:Union{Float64, ComplexF64}}
        L == length(A) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(A))"))
    end
    UniformMPS(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype) = UniformMPS{L, T}()
    UniformMPS(obj::CanonicalMPS) = canonicalize(obj)
end
const UMPS = UniformMPS

"""
    randUMPS([::Type{T},] pspace::Vector{VectorSpace}, aspace::Vector{VectorSpace}) -> UMPS{L}

Generate a length `L` random UMPS with given length `L` vector `pspace` and `aspace`. `T = ComplexF64`(default) or `Float64` is the number type.

    randUMPS([::Type{T},] L::Int64, pspace::VectorSpace, apsace::VectorSpace) -> UMPS{L}

Assume the same `pspace` and `aspace`, except for the boundary bond, which is assumed to be trivial.

    randUMPS([::Type{T},] pdim::Vector{Int64}, adim::Vector{Int64}) -> UMPS{L}

Assume `pspace = [spacetype] .^ pdim` and `aspace = [spacetype] .^ adim`, where `spacetype = ℂ(T=ComplexF64) or ℝ(T=Float64)`

    randUMPS([::Type{T},] L::Int64, D::Int64, d::Int64) -> UMPS{L}

Assume the same `pdim` and `adim`.
"""
function randUMPS(::Type{T}, pspace::Vector{VectorSpace}, aspace::Vector{VectorSpace}) where T<:Union{Float64, ComplexF64}
end
function randUMPS(::Type{T}, L::Int64, pspace::VectorSpace, apsace::VectorSpace) where T<:Union{Float64, ComplexF64}
end
function randUMPS(::Type{T}, pdim::Vector{Int64}, adim::Vector{Int64}) where T<:Union{Float64, ComplexF64}
end
function randUMPS(::Type{T}, L::Int64, pdim::Int64, adim::Int64) where T<:Union{Float64, ComplexF64}
end
