"""
    mutable struct CanonicalMPS{L, T<:Union{Float64, ComplexF64}} <: DenseCanonicalMPS{L, T}
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
- AL[i] * C[i] = AC[i] = C[i-1] * AR[i]
- AL[i]' * AL[i] = I
- AR[i] * AR[i]' = I
The first virtual space and C[1] are between AL[1] and AL[2], and the ones behind AL[L] same as the ones ahead AL[1] are the last ones.

# Constructors
    CanonicalMPS{L, T}(AL::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                       AR::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                       AC::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                       C::AbstractVector{<:AbstractBondTensor}=Vector{BondTensor}(undef, L))
    CanonicalMPS(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype)
    CanonicalMPS(AL::AbstractVector{<:AbstractMPSTensor},
                 AR::AbstractVector{<:AbstractMPSTensor},
                 AC::AbstractVector{<:AbstractMPSTensor},
                 C::AbstractVector{<:AbstractBondTensor})
    CanonicalMPS(AL::AbstractVector{<:AbstractMPSTensor}, C::AbstractVector{<:AbstractBondTensor})
    CanonicalMPS(AL::AbstractMPSTensor, C::AbstractBondTensor)
"""
mutable struct CanonicalMPS{L, T<:Union{Float64, ComplexF64}} <: DenseCanonicalMPS{L, T}
    const AL::AbstractVector{<:AbstractMPSTensor}
    const AR::AbstractVector{<:AbstractMPSTensor}
    const AC::AbstractVector{<:AbstractMPSTensor}
    const C::AbstractVector{<:AbstractBondTensor}

    function CanonicalMPS{L, T}(AL::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                                AR::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                                AC::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                                C::AbstractVector{<:AbstractBondTensor}=Vector{BondTensor}(undef, L)) where {L, T<:Union{Float64, ComplexF64}}
        (L == length(AL) == length(AR) == length(AC) == length(C)) || throw(ArgumentError("Mismatched lengths: ($L, $(length(AL)), $(length(AR)), $(length(AC)), $(length(C)))"))

        if T == ComplexF64
            for i in 1:L
                (!isassigned(AL, i) && eltype(AL[i]) != T) && (AL[i] *= one(T))
                (!isassigned(AR, i) && eltype(AR[i]) != T) && (AR[i] *= one(T))
                (!isassigned(AC, i) && eltype(AC[i]) != T) && (AC[i] *= one(T))
                (!isassigned(C, i) && eltype(C[i]) != T) && (C[i] *= one(T))
            end
        end

        return new{L, T}(AL, AR, AC, C)
    end
    CanonicalMPS(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype) = CanonicalMPS{L, T}()

    function CanonicalMPS(AL::AbstractVector{<:AbstractMPSTensor},
                          AR::AbstractVector{<:AbstractMPSTensor},
                          AC::AbstractVector{<:AbstractMPSTensor},
                          C::AbstractVector{<:AbstractBondTensor})
        L = length(AL)
        T = promote_type(mapreduce(eltype, promote_type, AL),
                         mapreduce(eltype, promote_type, AR),
                         mapreduce(eltype, promote_type, AC),
                         mapreduce(eltype, promote_type, C))
        return CanonicalMPS{L, T}(AL, AR, AC, C)
    end
    function CanonicalMPS(AL::AbstractVector{<:AbstractMPSTensor}, C::AbstractVector{<:AbstractBondTensor})
        (L = length(AL)) == length(C) || throw(ArgumentError("Mismatched lengths: $(length(AL)) ≠ $(length(C))"))
        T = promote_type(mapreduce(eltype, promote_type, AL),
                         mapreduce(eltype, promote_type, C))
        AC = AL .* C
        AR = inv.(C)[[L, (1:L-1)...]] * AC
        return CanonicalMPS{L, T}(AL, AR, AC, C)
    end
    CanonicalMPS(AL::AbstractMPSTensor, C::AbstractBondTensor) = CanonicalMPS([AL,], [C,])
end
const CMPS = CanonicalMPS
const InfiniteMPS{L, T} = Union{UniformMPS{L, T}, CanonicalMPS{L, T}}

"""
    uniformize(::Type{CanonicalMPS{L, T}}) -> ::Type{UniformMPS{L, T}}

Interface of `uniformize` indicates `uniformize(::CanonicalMPS) -> ::UniformMPS`.
"""
uniformize(::Type{CanonicalMPS{L, T}}) where {L, T} = UniformMPS{L, T}

"""
    randCMPS([::Type{T},] pspace::Vector{VectorSpace}, aspace::Vector{VectorSpace}; kwargs...) -> CMPS{L, T}
    randCMPS([::Type{T},] L::Int64, pspace::VectorSpace, apsace::VectorSpace; kwargs...) -> CMPS{L, T}
    randCMPS([::Type{T},] pdim::Vector{Int64}, adim::Vector{Int64}; kwargs...) -> CMPS{L, T}
    randCMPS([::Type{T},] L::Int64, pdim::Int64, adim::Int64; kwargs...) -> CMPS{L, T}

Generate a length `L` random CMPS by canonicalize a random UMPS. The arguments are propagated to `randUMPS` and `kwargs` are propagated to `canonicalize`.
"""
function randCMPS(::Type{T}, pspace::Vector{VectorSpace}, aspace::Vector{VectorSpace}; kwargs...) where T<:Union{Float64, ComplexF64}
    umps = randUMPS(T, pspace, aspace)
    return canonicalize(umps; kwargs...)
end
function randCMPS(::Type{T}, L::Int64, pspace::VectorSpace, apsace::VectorSpace; kwargs...) where T<:Union{Float64, ComplexF64}
    umps = randUMPS(T, L, pspace, aspace)
    return canonicalize(umps; kwargs...)
end
function randCMPS(::Type{T}, pdim::Vector{Int64}, adim::Vector{Int64}; kwargs...) where T<:Union{Float64, ComplexF64}
    umps = randUMPS(T, pdim, adim)
    return canonicalize(umps; kwargs...)
end
function randCMPS(::Type{T}, L::Int64, pdim::Int64, adim::Int64; kwargs...) where T<:Union{Float64, ComplexF64}
    umps = randUMPS(T, L, pdim, adim)
    return canonicalize(umps; kwargs...)
end
