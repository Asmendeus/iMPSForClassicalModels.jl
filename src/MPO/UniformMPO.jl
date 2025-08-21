"""
    mutable struct UniformMPO{L, T<:Union{Float64, ComplexF64}} <: DenseUniformMPO{L, T}
        const A::AbstractVector{<:AbstractMPOTensor}
    end

Concrete type of iMPO with uniform form, where `L` is the cell length, `T == Float64` or `ComplexF64` is the number type of local tensors.

All the fields and constructors are exactly the same to those of `UniformMPS`.
We redefine the type `UniformMPO` just for using multiple dispacth when implementing the algebra between `UniformMPS` and `UniformMPO`.
Details of constructors please see `UniformMPS`.
"""
mutable struct UniformMPO{L, T<:Union{Float64, ComplexF64}} <: DenseUniformMPO{L, T}
    const A::AbstractVector{<:AbstractMPOTensor}

    function UniformMPO{L, T}(A::AbstractVector{<:AbstractMPOTensor}=Vector{MPOTensor}(undef, L)) where {L, T<:Union{Float64, ComplexF64}}
        L == length(A) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(A))"))
        if T == ComplexF64
            for i in 1:L
                (!isassigned(A, i) && eltype(A[i]) != T) && (A[i] *= one(T))
            end
        end
        return new{L, T}(A)
    end
    UniformMPO(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype) = UniformMPO{L, T}()

    function UniformMPO(A::AbstractVector{<:AbstractMPOTensor})
        L = length(A)
        T = mapreduce(eltype, promote_type, A)
        return new{L, T}(A)
    end
end
const UMPO = UniformMPO

"""
    canonicalize(::Type{UniformMPO{L, T}}) -> ::Type{CanonicalMPO{L, T}}

Interface of `canonicalize` indicates `canonicalize(::UniformMPO) -> ::CanonicalMPO`.
"""
canonicalize(::Type{UniformMPO{L, T}}) where {L, T} = CanonicalMPO{L, T}

"""
    identityUMPO([::Type{T},] pspace::Vector{VectorSpace}, expspace::Vector{VectorSpace}=pspace) -> UMPO{L}

Generate a length `L` identity UMPO with given length `L` vector `pspace` and `expspace`.

    identityUMPO([::Type{T},] L::Int64, pspace::VectorSpace, expspace::VectorSpace=pspace) -> UMPO{L}

Assume the same `pspace` and `expspace`, except for the boundary bond, which is assumed to be trivial.

    identityUMPO([::Type{T},] pdim::Vector{Int64}, expdim::Vector{Int64}=pdim) -> UMPO{L}

Assume `pspace = [spacetype] .^ pdim` and `expspace = [spacetype] .^ expdim`, where `spacetype = ℂ(T=ComplexF64) or ℝ(T=Float64)`

    identityUMPO([::Type{T},] L::Int64, pdim::Int64, expdim::Int64=pdim) -> UMPO{L}

Assume the same `pdim` and `adim`.
"""
function identityUMPO(::Type{T}, pspace::Vector{VectorSpace}, expspace::Vector{VectorSpace}=pspace) where T<:Union{Float64, ComplexF64}
    (L = length(pspace)) == length(expspace) || throw(ArgumentError("Mismatched lengths: $(length(pspace)) ≠ $(length(expspace))"))

    aspace = trivial(pspace[1])
    va = TensorMap(ones, aspace, aspace)

    A = Vector{MPOTensor}(undef, L)
    for i in 1:L
        Id = TensorMap(Matrix(I, dim(pspace[i], dim(expspace[i]))) * one(T), pspace[i], expspace[i])
        @tensor Ai[-1 -2; -3 -4] := va[-1 -4] * Id[-2 -3]
        A[i] = Ai
    end

    return UniformMPO{L, T}(A)
end
function identityUMPO(::Type{T}, L::Int64, pspace::VectorSpace, expspace::VectorSpace=pspace) where T<:Union{Float64, ComplexF64}
    aspace = trivial(pspace)
    va = TensorMap(ones, aspace, aspace)

    Id = TensorMap(Matrix(I, dim(pspace), dim(expspace)) * one(T), pspace, expspace)
    @tensor A[-1 -2; -3 -4] := va[-1 -4] * Id[-2 -3]

    A = repeat([MPOTensor(A),], L)

    return UniformMPO{L, T}(A)
end
function identityUMPO(::Type{T}, pdim::Vector{Int64}, expdim::Vector{Int64}=pdim) where T<:Union{Float64, ComplexF64}
    (L = length(pdim)) == length(expdim) || throw(ArgumentError("Mismatched lengths: $(length(pdim)) ≠ $(length(expdim))"))

    spacetype = T == ComplexF64 ? ℂ : ℝ

    aspace = spacetype ^ 1
    va = TensorMap(ones, aspace, aspace)

    A = Vector{MPOTensor}(undef, L)
    pspace = map(x->spacetype ^ x, pdim)
    expspace = map(x->spacetype ^ x, expdim)
    for i in 1:L
        Id = TensorMap(Matrix(I, pdim[i], expdim[i]) * one(T), pspace[i], expspace[i])
        @tensor Ai[-1 -2; -3 -4] := va[-1 -4] * Id[-2 -3]
        A[i] = Ai
    end

    return UniformMPO{L, T}(A)
end
function identityUMPO(::Type{T}, L::Int64, pdim::Int64, expdim::Int64=pdim) where T<:Union{Float64, ComplexF64}
    spacetype = T == ComplexF64 ? ℂ : ℝ

    aspace = spacetype ^ 1
    va = TensorMap(ones, aspace, aspace)

    pspace = map(x->spacetype ^ x, pdim)
    expspace = map(x->spacetype ^ x, expdim)
    Id = TensorMap(Matrix(I, pdim, expdim) * one(T), pspace, expspace)
    @tensor A[-1 -2; -3 -4] := va[-1 -4] * Id[-2 -3]

    A = repeat([MPOTensor(A),], L)

    return UniformMPO{L, T}(A)
end
