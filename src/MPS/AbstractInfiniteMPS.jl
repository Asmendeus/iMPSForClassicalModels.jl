# ===== Infinite =====
"""
    abstract type AbstractInfiniteMPS{L}

Abstract type of all iMPS/iMPO of cell length `L`.
"""
abstract type AbstractInfiniteMPS{L} end

length(::AbstractInfiniteMPS{L}) where L = L

# ===== Uniform <: Infinite =====
"""
    abstract type AbstractUniformMPS{L} <: AbstractInfiniteMPS{L}

Abstract type of all iMPS/iMPO with uniform form of cell length `L`.

Note each concrete subtype must have a field:
    `A::AbstractVector{<:AbstractLocalTensor{R}}`: local tensors.
"""
abstract type AbstractUniformMPS{L} <: AbstractInfiniteMPS{L} end

for func in (:getindex, :lastindex, :setindex!, :iterate, :keys, :isassigned)
    @eval Base.$func(obj::AbstractUniformMPS, args...) = $func(obj.A, args...)
end

"""
    getA(obj::AbstractUniformMPS)

Interface of `AbstractUniformMPS`, return the local tensors
"""
getA(obj::AbstractUniformMPS) = obj.A

"""
    norm(obj::AbstractUniformMPS;
        XL₀::AbstractVector{<:BondTensor}=_default_X₀_leftFixedPoint(obj.A),
        alg::EigenAlgorithm=Defaults.alg_eig) -> ::Vector{Float64}

Return the largest eigenvalues in canonicalization, i.e., `λ` solved by

    -- X[l] -- A[l] --  =  λ[l] * -- AL[l] -- X[l+1] --
                |                     |
"""
function norm(obj::AbstractUniformMPS;
            XL₀::AbstractVector{<:BondTensor}=_default_X₀_leftFixedPoint(obj.A),
            alg::EigenAlgorithm=Defaults.alg_eig)
    λ, _, _ = leftFixedPoint(obj.A, XL₀, alg)
    @assert all(x->abs(imag(x)) < Defaults.tol, λ) "Incorrect result: norms of non-real number"
    return abs.(λ)
end

"""
    norm(obj::AbstractUniformMPS, si::Int64;
        XL₀::AbstractVector{<:BondTensor}=_default_X₀_leftFixedPoint(obj.A),
        alg::EigenAlgorithm=Defaults.alg_eig) -> ::Float64

Return the `si`-th largest eigenvalue in canonicalization.
"""
function norm(obj::AbstractUniformMPS, si::Int64;
            XL₀::AbstractVector{<:BondTensor}=_default_X₀_leftFixedPoint(obj.A),
            alg::EigenAlgorithm=Defaults.alg_eig)
    return norm(obj; XL₀=XL₀, alg=alg)[si]
end

"""
    normalize!(obj::AbstractUniformMPS{L};
            XL₀::AbstractVector{<:BondTensor}=_default_X₀_leftFixedPoint(obj.A),
            alg::EigenAlgorithm=Defaults.alg_eig)

Normalize a given uniform iMPS according to inner-induced norm.
"""
function normalize!(obj::AbstractUniformMPS{L};
            XL₀::AbstractVector{<:BondTensor}=_default_X₀_leftFixedPoint(obj.A),
            alg::EigenAlgorithm=Defaults.alg_eig) where L
    λ = norm(obj; XL₀=XL₀, alg=alg)
    obj.A ./= λ
    return obj
end

"""
    normalize(obj::AbstractUniformMPS{L};
            XL₀::AbstractVector{<:BondTensor}=_default_X₀_leftFixedPoint(obj.A),
            alg::EigenAlgorithm=Defaults.alg_eig)

Pure function of `normalize!`
"""
function normalize(obj::AbstractUniformMPS{L};
        XL₀::AbstractVector{<:BondTensor}=_default_X₀_leftFixedPoint(obj.A),
        alg::EigenAlgorithm=Defaults.alg_eig) where L
    obj′ = deepcopy(obj)
    return normalize!(obj′; XL₀=XL₀, alg=alg)
end

# space
leftVirtualSpace(obj::AbstractUniformMPS) = map(x->leftVirtualSpace(x), obj.A)
rightVirtualSpace(obj::AbstractUniformMPS) = map(x->rightVirtualSpace(x), obj.A)
physicalSpace(obj::AbstractUniformMPS) = map(x->physicalSpace(x), obj.A)
extraPhysicalSpace(obj::AbstractUniformMPS) = map(x->extraPhysicalSpace(x), obj.A)

# ===== Canonical <: Infinite =====
"""
    abstract type AbstractCanonicalMPS{L} <: AbstractInfiniteMPS{L}

Abstract type of all iMPS/iMPO with canonical form of cell length `L` and data type `T`.

Note each concrete subtype must have fields:
    `AL::AbstractVector{<:AbstractLocalTensor{R}}`: left-canonical tensors,
    `AR::AbstractVector{<:AbstractLocalTensor{R}}`: right-canonical tensors,
    `AC::AbstractVector{<:AbstractLocalTensor{R}}`: center tensors,
    `C::AbstractVector{<:AbstractBondTensor}`: center bond tensors.
"""
abstract type AbstractCanonicalMPS{L} <: AbstractInfiniteMPS{L} end

"""
    getAL(obj::AbstractCanonicalMPS)

Interface of `AbstractCanonicalMPS`, return the left-canonical tensors
"""
getAL(obj::AbstractCanonicalMPS) = obj.AL

"""
    getAR(obj::AbstractCanonicalMPS)

Interface of `AbstractCanonicalMPS`, return the right-canonical tensors
"""
getAR(obj::AbstractCanonicalMPS) = obj.AR

"""
    getAC(obj::AbstractCanonicalMPS)

Interface of `AbstractCanonicalMPS`, return the center tensors
"""
getAC(obj::AbstractCanonicalMPS) = obj.AC

"""
    getC(obj::AbstractCanonicalMPS)

Interface of `AbstractCanonicalMPS`, return the center bond tensors
"""
getC(obj::AbstractCanonicalMPS) = obj.C

"""
    norm(obj::AbstractCanonicalMPS) -> ::Vector{Float64}

Return the inner-induced norms.
"""
norm(obj::AbstractCanonicalMPS) = norm.(obj.C)

"""
    norm(obj::AbstractCanonicalMPS, si::Int64) -> ::Float64

Return the `si`-th inner-induced norm.
"""
norm(obj::AbstractCanonicalMPS, si::Int64) = norm(obj.C[si])

"""
    normalize!(obj::AbstractCanonicalMPS)

Normalize a given canonical iMPS according to inner-induced norm.
"""
function normalize!(obj::AbstractCanonicalMPS)
    nc = norm(obj)
    obj.C ./= nc
    obj.AC ./= nc
    return obj
end

"""
    normalize(obj::AbstractCanonicalMPS)

Pure function of `normalize!`
"""
function normalize(obj::AbstractCanonicalMPS)
    obj′ = deepcopy(obj)
    return normalize!(obj′)
end

# space
leftVirtualSpace(obj::AbstractCanonicalMPS) = map(x->leftVirtualSpace(x), obj.AL)
rightVirtualSpace(obj::AbstractCanonicalMPS) = map(x->rightVirtualSpace(x), obj.AL)
physicalSpace(obj::AbstractCanonicalMPS) = map(x->physicalSpace(x), obj.AL)
extraPhysicalSpace(obj::AbstractCanonicalMPS) = map(x->extraPhysicalSpace(x), obj.AL)

# ===== Dense Uniform <: Uniform =====
"""
    abstract type DenseUniformMPS{L, T<:Union{Float64, ComplexF64}} <: AbstractUniformMPS{L}

Abstract type of all dense iMPS/iMPO with uniform form of cell length `L`.
"""
abstract type DenseUniformMPS{L, T<:Union{Float64, ComplexF64}} <: AbstractUniformMPS{L} end

scalartype(::DenseUniformMPS{L, T}) where {L, T} = T

"""
    canonicalize(::Type{DenseUniformMPS}) -> ::Type{DenseCanonicalMPS}

Interface of `canonicalize(::DenseUniformMPS)` indicates the data type after transformation.
"""
function canonicalize(obj::Type{DenseUniformMPS})
    throw(ArgumentError("Not define `canonicalize` for `$(typeof(obj))`"))
end

function Base.show(io::IO, obj::DenseUniformMPS{L}) where L
    any(i -> !isassigned(obj.A, i), 1:L) && return println(io, "$(typeof(obj)): L = $L, to be initialized !")

    # avoid to show bond info
    memory = obj |> Base.summarysize |> Base.format_bytes
    println(io, "$(typeof(obj)): total memory = $memory")
    # bond dimenson
    lsi = ceil(Int64, log10(L)) # length of si to be printed
    for si = 1:L
        local A = obj.A[si]
        D, DD = dim(A, 1)
        println(io, "Bond ", lpad(si-1, lsi), "->", lpad(si, lsi), ": $(codomain(A).spaces[1]), dim = $(D) -> $(DD)")
    end
end

# ===== Dense Canonical <: Canonical =====
"""
    abstract type DenseCanonicalMPS{L, T<:Union{Float64, ComplexF64}} <: AbstractCanonicalMPS{L}

Abstract type of all dense iMPS/iMPO with canonical form of cell length `L`.
"""
abstract type DenseCanonicalMPS{L, T<:Union{Float64, ComplexF64}} <: AbstractCanonicalMPS{L} end

scalartype(::DenseCanonicalMPS{L, T}) where {L, T} = T

"""
    uniformize(::Type{DenseCanonicalMPS}) -> ::Type{DenseUniformMPS}

Interface of `uniformize(::DenseCanonicalMPS)` indicates the data type after transformation.
"""
function uniformize(obj::Type{DenseCanonicalMPS})
    throw(ArgumentError("Not define `uniformize` for `$(typeof(obj))`"))
end

function Base.show(io::IO, obj::DenseCanonicalMPS{L}) where L
    any(i -> !isassigned(obj.AL, i), 1:L) && return println(io, "$(typeof(obj)): L = $L, to be initialized !")
    any(i -> !isassigned(obj.AR, i), 1:L) && return println(io, "$(typeof(obj)): L = $L, to be initialized !")
    any(i -> !isassigned(obj.AC, i), 1:L) && return println(io, "$(typeof(obj)): L = $L, to be initialized !")
    any(i -> !isassigned(obj.C, i), 1:L) && return println(io, "$(typeof(obj)): L = $L, to be initialized !")

    # avoid to show bond info
    memory = obj |> Base.summarysize |> Base.format_bytes
    println(io, "$(typeof(obj)): total memory = $memory")
    # bond dimenson
    lsi = ceil(Int64, log10(L)) # length of si to be printed
    for si = 1:L
        local A = obj.AC[si]
        D, DD = dim(A, 1)
        println(io, "Bond ", lpad(si-1, lsi), "->", lpad(si, lsi), ": $(codomain(A).spaces[1]), dim = $(D) -> $(DD)")
    end
end

const DenseInfiniteMPS{L, T} = Union{DenseUniformMPS{L, T}, DenseCanonicalMPS{L, T}}

# Interface functions for `iterate`
lambda(obj::DenseInfiniteMPS) = norm.(obj.C) .* sign.(obj.C)
function division(obj::DenseInfiniteMPS, n::AbstractVector{<:Number})
    obj′ = deepcopy(obj)
    obj′.C ./= n
    obj′.AC ./= n
    return obj′
end
minus(obj1::T, _::T) where T <: DenseInfiniteMPS = obj1
norm_max(obj::DenseInfiniteMPS) = maximum(norm.(obj.AL .* obj.C .- obj.AC))
