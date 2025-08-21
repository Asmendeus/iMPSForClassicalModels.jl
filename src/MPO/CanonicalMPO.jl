"""
    mutable struct CanonicalMPO{L, T<:Union{Float64, ComplexF64}} <: DenseCanonicalMPO{L, T}
        const AL::AbstractVector{<:AbstractMPOTensor}
        const AR::AbstractVector{<:AbstractMPOTensor}
        const AC::AbstractVector{<:AbstractMPOTensor}
        const C::AbstractVector{<:AbstractBondTensor}
    end

Concrete type of iMPO with canonical form, where `L` is the cell length, `T == Float64` or `ComplexF64` is the number type of local tensors.

All the fields and constructors are exactly the same to those of `CanonicalMPS`.
We redefine the type `CanonicalMPO` just for using multiple dispacth when implementing the algebra between `UniformMPS` and `UniformMPO`.
Details of constructors please see `CanonicalMPS`.
"""
mutable struct CanonicalMPO{L, T<:Union{Float64, ComplexF64}} <: DenseCanonicalMPO{L, T}
    const AL::AbstractVector{<:AbstractMPOTensor}
    const AR::AbstractVector{<:AbstractMPOTensor}
    const AC::AbstractVector{<:AbstractMPOTensor}
    const C::AbstractVector{<:AbstractBondTensor}

    function CanonicalMPO{L, T}(AL::AbstractVector{<:AbstractMPOTensor}=Vector{MPOTensor}(undef, L),
                                AR::AbstractVector{<:AbstractMPOTensor}=Vector{MPOTensor}(undef, L),
                                AC::AbstractVector{<:AbstractMPOTensor}=Vector{MPOTensor}(undef, L),
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
    CanonicalMPO(L::Int64, T::Union{Float64, ComplexF64}=Defaults.datatype) = CanonicalMPO{L, T}()

    function CanonicalMPO(AL::AbstractVector{<:AbstractMPOTensor},
                          AR::AbstractVector{<:AbstractMPOTensor},
                          AC::AbstractVector{<:AbstractMPOTensor},
                          C::AbstractVector{<:AbstractBondTensor})
        L = length(AL)
        T = promote_type(mapreduce(eltype, promote_type, AL),
                         mapreduce(eltype, promote_type, AR),
                         mapreduce(eltype, promote_type, AC),
                         mapreduce(eltype, promote_type, C))
        return CanonicalMPO{L, T}(AL, AR, AC, C)
    end
    function CanonicalMPO(AL::AbstractVector{<:AbstractMPOTensor}, C::AbstractVector{<:AbstractBondTensor})
        (L = length(AL)) == length(C) || throw(ArgumentError("Mismatched lengths: $(length(AL)) ≠ $(length(C))"))
        T = promote_type(mapreduce(eltype, promote_type, AL),
                         mapreduce(eltype, promote_type, C))
        AC = AL .* C
        AR = inv.(C)[[end, (1:end-1)...]] * AC
        return CanonicalMPO{L, T}(AL, AR, AC, C)
    end
end
const CMPO = CanonicalMPO

"""
    uniformize(::Type{CanonicalMPO{L, T}}) -> ::Type{UniformMPO{L, T}}

Interface of `uniformize` indicates `uniformize(::CanonicalMPO) -> ::UniformMPO`.
"""
uniformize(::Type{CanonicalMPO{L, T}}) where {L, T} = UniformMPO{L, T}

"""
    identityCMPO([::Type{T},] pspace::Vector{VectorSpace}, expspace::Vector{VectorSpace}=pspace) -> CMPO{L}

Generate a length `L` identity CMPO with given length `L` vector `pspace` and `expspace`.

    identityCMPO([::Type{T},] L::Int64, pspace::VectorSpace, expspace::VectorSpace=pspace) -> CMPO{L}

Assume the same `pspace` and `expspace`, except for the boundary bond, which is assumed to be trivial.

    identityCMPO([::Type{T},] pdim::Vector{Int64}, expdim::Vector{Int64}=pdim) -> CMPO{L}

Assume `pspace = [spacetype] .^ pdim` and `expspace = [spacetype] .^ expdim`, where `spacetype = ℂ(T=ComplexF64) or ℝ(T=Float64)`

    identityCMPO([::Type{T},] L::Int64, pdim::Int64, expdim::Int64=pdim) -> CMPO{L}

Assume the same `pdim` and `adim`.
"""
function identityCMPO(::Type{T}, pspace::Vector{VectorSpace}, expspace::Vector{VectorSpace}=pspace) where T<:Union{Float64, ComplexF64}
    (L = length(pspace)) == length(expspace) || throw(ArgumentError("Mismatched lengths: $(length(pspace)) ≠ $(length(expspace))"))

    aspace = trivial(pspace[1])
    va = TensorMap(ones, aspace, aspace)

    A = Vector{MPOTensor}(undef, L)
    for i in 1:L
        Id = TensorMap(Matrix(I, dim(pspace[i], dim(expspace[i]))) * one(T), pspace[i], expspace[i])
        @tensor Ai[-1 -2; -3 -4] := va[-1 -4] * Id[-2 -3]
        A[i] = Ai
    end

    return CanonicalMPO{L, T}(A, A, A, repeat([BondTensor(va),], L))
end
function identityCMPO(::Type{T}, L::Int64, pspace::VectorSpace, expspace::VectorSpace=pspace) where T<:Union{Float64, ComplexF64}
    aspace = trivial(pspace)
    va = TensorMap(ones, aspace, aspace)

    Id = TensorMap(Matrix(I, dim(pspace), dim(expspace)) * one(T), pspace, expspace)
    @tensor A[-1 -2; -3 -4] := va[-1 -4] * Id[-2 -3]

    A = repeat([MPOTensor(A),], L)

    return CanonicalMPO{L, T}(A, A, A, repeat([BondTensor(va),], L))
end
function identityCMPO(::Type{T}, pdim::Vector{Int64}, expdim::Vector{Int64}=pdim) where T<:Union{Float64, ComplexF64}
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

    return CanonicalMPO{L, T}(A, A, A, repeat([BondTensor(va),], L))
end
function identityCMPO(::Type{T}, L::Int64, pdim::Int64, expdim::Int64=pdim) where T<:Union{Float64, ComplexF64}
    spacetype = T == ComplexF64 ? ℂ : ℝ

    aspace = spacetype ^ 1
    va = TensorMap(ones, aspace, aspace)

    pspace = map(x->spacetype ^ x, pdim)
    expspace = map(x->spacetype ^ x, expdim)
    Id = TensorMap(Matrix(I, pdim, expdim) * one(T), pspace, expspace)
    @tensor A[-1 -2; -3 -4] := va[-1 -4] * Id[-2 -3]

    A = repeat([MPOTensor(A),], L)

    return CanonicalMPO{L, T}(A, A, A, repeat([BondTensor(va),], L))
end
