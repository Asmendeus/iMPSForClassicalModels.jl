"""
    mutable struct UniformMPS{L, T<:Union{Float64, ComplexF64}} <: DenseUniformMPS{L, T}
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
    UniformMPS(A::AbstractVector{<:AbstractMPSTensor})
"""
mutable struct UniformMPS{L, T<:Union{Float64, ComplexF64}} <: DenseUniformMPS{L, T}
    const A::AbstractVector{<:AbstractMPSTensor}

    function UniformMPS{L, T}(A::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L)) where {L, T<:Union{Float64, ComplexF64}}
        L == length(A) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(A))"))
        if T == ComplexF64
            for i in 1:L
                (!isassigned(A, i) && eltype(A[i]) != T) && (A[i] *= one(T))
            end
        end
        return new{L, T}(A)
    end
    UniformMPS(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype) = UniformMPS{L, T}()

    function UniformMPS(A::AbstractVector{<:AbstractMPSTensor})
        L = length(A)
        T = mapreduce(eltype, promote_type, A)
        return new{L, T}(A)
    end
end
const UMPS = UniformMPS

"""
    canonicalize(::Type{UniformMPS{L, T}}) -> ::Type{CanonicalMPS{L, T}}

Interface of `canonicalize` indicates `canonicalize(::UniformMPS) -> ::CanonicalMPS`.
"""
canonicalize(::Type{UniformMPS{L, T}}) where {L, T} = CanonicalMPS{L, T}

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
    (L = length(pspace)) == length(aspace) || throw(ArgumentError("Mismatched lengths: $(length(pspace)) ≠ $(length(aspace))"))
    left_aspace = aspace[[end, (1:end-1)...]]
    right_aspace = aspace
    A = map(l->TensorMap(rand, T, left_aspace[l]⊗pspace[l], right_aspace[l]), 1:L)
    return UniformMPS{L, T}(A)
end
function randUMPS(::Type{T}, L::Int64, pspace::VectorSpace, apsace::VectorSpace) where T<:Union{Float64, ComplexF64}
    A = map(_->TensorMap(rand, T, aspace⊗pspace, aspace), 1:L)
    return UniformMPS{L, T}(A)
end
function randUMPS(::Type{T}, pdim::Vector{Int64}, adim::Vector{Int64}) where T<:Union{Float64, ComplexF64}
    (L = length(pdim)) == length(adim) || throw(ArgumentError("Mismatched lengths: $(length(pdim)) ≠ $(length(adim))"))
    spacetype = T == ComplexF64 ? ℂ : ℝ
    pspace = map(x->spacetype ^ x, pdim)
    left_aspace = map(x->spacetype ^ x, adim[[end, (1:end-1)...]])
    right_aspace = map(x->spacetype ^ x, adim)
    A = map(l->TensorMap(rand, T, left_aspace[l]⊗pspace[l], right_aspace[l]), 1:L)
    return UniformMPS{L, T}(A)
end
function randUMPS(::Type{T}, L::Int64, pdim::Int64, adim::Int64) where T<:Union{Float64, ComplexF64}
    spacetype = T == ComplexF64 ? ℂ : ℝ
    pspace = spacetype ^ pdim
    aspace = spacetype ^ adim
    A = map(_->TensorMap(rand, T, aspace⊗pspace, aspace), 1:L)
    return UniformMPS{L, T}(A)
end
