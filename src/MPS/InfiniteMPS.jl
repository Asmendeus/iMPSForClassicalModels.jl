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

# Notes
By convention, we appoint that:
- AL[i] * C[i] = AC[i] = C[i-1] * AR[i]
- AL[i]' * AL[i] = I
- AR[i] * AR[i]' = I
The first virtual space and C[1] are between AL[1] and AL[2], and the ones behind AL[L] same as the ones ahead AL[1] are the last ones.

# Constructors
    InfiniteMPS{L, T}(AL::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                       AR::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                       AC::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
                       C::AbstractVector{<:AbstractBondTensor}=Vector{BondTensor}(undef, L))
    InfiniteMPS(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype)

    InfiniteMPS(AL::AbstractVector{<:AbstractMPSTensor},
                 AR::AbstractVector{<:AbstractMPSTensor},
                 AC::AbstractVector{<:AbstractMPSTensor},
                 C::AbstractVector{<:AbstractBondTensor})

    InfiniteMPS(AL::AbstractVector{<:AbstractMPSTensor}, C::AbstractVector{<:AbstractBondTensor})
    InfiniteMPS(AL::AbstractMPSTensor, C::AbstractBondTensor)

    InfiniteMPS(A::AbstractVector{<:AbstractMPSTensor})
    InfiniteMPS(A::AbstractMPSTensor)
"""
mutable struct InfiniteMPS{L, T<:Union{Float64, ComplexF64}} <: DenseInfiniteMPS{L, T}
    const AL::AbstractVector{<:AbstractMPSTensor}
    const AR::AbstractVector{<:AbstractMPSTensor}
    const AC::AbstractVector{<:AbstractMPSTensor}
    const C::AbstractVector{<:AbstractBondTensor}

    function InfiniteMPS{L, T}(AL::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L),
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
    InfiniteMPS(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype) = InfiniteMPS{L, T}()

    function InfiniteMPS(AL::AbstractVector{<:AbstractMPSTensor},
                          AR::AbstractVector{<:AbstractMPSTensor},
                          AC::AbstractVector{<:AbstractMPSTensor},
                          C::AbstractVector{<:AbstractBondTensor})
        L = length(AL)
        T = promote_type(mapreduce(eltype, promote_type, AL),
                         mapreduce(eltype, promote_type, AR),
                         mapreduce(eltype, promote_type, AC),
                         mapreduce(eltype, promote_type, C))
        return InfiniteMPS{L, T}(AL, AR, AC, C)
    end
    function InfiniteMPS(AL::AbstractVector{<:AbstractMPSTensor}, C::AbstractVector{<:AbstractBondTensor})
        (L = length(AL)) == length(C) || throw(ArgumentError("Mismatched lengths: $(length(AL)) ≠ $(length(C))"))
        T = promote_type(mapreduce(eltype, promote_type, AL),
                         mapreduce(eltype, promote_type, C))
        AC = AL .* C
        AR = inv.(C)[[L, (1:L-1)...]] * AC
        return InfiniteMPS{L, T}(AL, AR, AC, C)
    end
    InfiniteMPS(AL::AbstractMPSTensor, C::AbstractBondTensor) = InfiniteMPS([AL,], [C,])

    InfiniteMPS(A::AbstractVector{<:AbstractMPSTensor}) = canonicalize(A)
    InfiniteMPS(A::AbstractMPSTensor) = canonicalize([A,])
end
const iMPS = InfiniteMPS

"""
    randInfiniteMPS([::Type{T}=Defaults.datatype,] pspace::Vector{VectorSpace}, aspace::Vector{VectorSpace}; kwargs...) -> iMPS{L, T}
    randInfiniteMPS([::Type{T}=Defaults.datatype,] L::Int64, pspace::VectorSpace, apsace::VectorSpace; kwargs...) -> iMPS{L, T}
    randInfiniteMPS([::Type{T}=Defaults.datatype,] pdim::Vector{Int64}, adim::Vector{Int64}; kwargs...) -> iMPS{L, T}
    randInfiniteMPS([::Type{T}=Defaults.datatype,] L::Int64, pdim::Int64, adim::Int64; kwargs...) -> iMPS{L, T}

Generate a length `L` random iMPS by canonicalize a random UMPS. `kwargs` are propagated to `canonicalize`.
"""
function randInfiniteMPS(::Type{T}, pspace::Vector{VectorSpace}, aspace::Vector{VectorSpace}; kwargs...) where T<:Union{Float64, ComplexF64}
    (L = length(pspace)) == length(aspace) || throw(ArgumentError("Mismatched lengths: $(length(pspace)) ≠ $(length(aspace))"))
    left_aspace = aspace[[L, (1:L-1)...]]
    right_aspace = aspace
    A = map(l->MPSTensor(TensorMap(rand, T, left_aspace[l]⊗pspace[l], right_aspace[l])), 1:L)
    return canonicalize(A; kwargs...)
end
randInfiniteMPS(pspace::Vector{VectorSpace}, aspace::Vector{VectorSpace}; kwargs...) = randInfiniteMPS(Defaults.datatype, pspace, aspace; kwargs...)

function randInfiniteMPS(::Type{T}, L::Int64, pspace::VectorSpace, apsace::VectorSpace; kwargs...) where T<:Union{Float64, ComplexF64}
    A = map(_->MPSTensor(TensorMap(rand, T, aspace⊗pspace, aspace)), 1:L)
    return canonicalize(A; kwargs...)
end
randInfiniteMPS(L::Int64, pspace::VectorSpace, apsace::VectorSpace; kwargs...) = randInfiniteMPS(Defaults.datatype, L, pspace, apsace; kwargs...)

function randInfiniteMPS(::Type{T}, pdim::Vector{Int64}, adim::Vector{Int64}; kwargs...) where T<:Union{Float64, ComplexF64}
    (L = length(pdim)) == length(adim) || throw(ArgumentError("Mismatched lengths: $(length(pdim)) ≠ $(length(adim))"))
    spacetype = T == ComplexF64 ? ℂ : ℝ
    pspace = map(x->spacetype ^ x, pdim)
    left_aspace = map(x->spacetype ^ x, adim[[L, (1:L-1)...]])
    right_aspace = map(x->spacetype ^ x, adim)
    A = map(l->MPSTensor(TensorMap(rand, T, left_aspace[l]⊗pspace[l], right_aspace[l])), 1:L)
    return canonicalize(A; kwargs...)
end
randInfiniteMPS(pdim::Vector{Int64}, adim::Vector{Int64}; kwargs...) = randInfiniteMPS(Defaults.datatype, pdim, adim; kwargs...)

function randInfiniteMPS(::Type{T}, L::Int64, pdim::Int64, adim::Int64; kwargs...) where T<:Union{Float64, ComplexF64}
    spacetype = T == ComplexF64 ? ℂ : ℝ
    pspace = spacetype ^ pdim
    aspace = spacetype ^ adim
    A = map(_->MPSTensor(TensorMap(rand, T, aspace⊗pspace, aspace)), 1:L)
    return canonicalize(A; kwargs...)
end
randInfiniteMPS(L::Int64, pdim::Int64, adim::Int64; kwargs...) = randInfiniteMPS(Defaults.datatype, L, pdim, adim; kwargs...)
