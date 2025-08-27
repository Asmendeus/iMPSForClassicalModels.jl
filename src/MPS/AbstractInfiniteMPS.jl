
# ========== General iMPS/iMPO ==========
"""
    abstract type AbstractGeneralInfiniteMPS{W, L}

Abstract type of all general iMPS/iMPO of cell width `W` and length `L`.
"""
abstract type AbstractGeneralInfiniteMPS{W, L} end

size(::AbstractGeneralInfiniteMPS{W, L}) where {W, L} = (W, L)

# ========== iMPS/iMPO <: General iMPS/iMPO ==========
"""
    abstract type AbstractInfiniteMPS{L} <: AbstractGeneralInfiniteMPS{1, L}

Abstract type of all iMPS/iMPO with canonical form of cell length `L`.
"""
abstract type AbstractInfiniteMPS{L} <: AbstractGeneralInfiniteMPS{1, L} end

length(::AbstractInfiniteMPS{L}) where L = L

# ========== Dense iMPS/iMPO <: iMPS/iMPO ==========
"""
    abstract type DenseInfiniteMPS{L, T<:Union{Float64, ComplexF64}} <: AbstractInfiniteMPS{L}

Abstract type of all dense iMPS/iMPO with canonical form of cell length `L`.
"""
abstract type DenseInfiniteMPS{L, T<:Union{Float64, ComplexF64}} <: AbstractInfiniteMPS{L} end

scalartype(::DenseInfiniteMPS{L, T}) where {L, T} = T

"""
    getAL(obj::DenseInfiniteMPS)

Interface of `DenseInfiniteMPS`, return the left-canonical tensors
"""
getAL(obj::DenseInfiniteMPS) = obj.AL

"""
    getAR(obj::DenseInfiniteMPS)

Interface of `DenseInfiniteMPS`, return the right-canonical tensors
"""
getAR(obj::DenseInfiniteMPS) = obj.AR

"""
    getAC(obj::DenseInfiniteMPS)

Interface of `DenseInfiniteMPS`, return the center tensors
"""
getAC(obj::DenseInfiniteMPS) = obj.AC

"""
    getC(obj::DenseInfiniteMPS)

Interface of `DenseInfiniteMPS`, return the center bond tensors
"""
getC(obj::DenseInfiniteMPS) = obj.C

"""
    norm(obj::DenseInfiniteMPS) -> ::Vector{Float64}

Return the inner-induced norms.
"""
norm(obj::DenseInfiniteMPS) = norm.(obj.C)

"""
    norm(obj::DenseInfiniteMPS, si::Int64) -> ::Float64

Return the `si`-th inner-induced norm.
"""
norm(obj::DenseInfiniteMPS, si::Int64) = norm(obj.C[si])

"""
    normalize!(obj::DenseInfiniteMPS)

Normalize a given canonical iMPS according to inner-induced norm.
"""
function normalize!(obj::DenseInfiniteMPS)
    nc = norm(obj)
    obj.C ./= nc
    obj.AC ./= nc
    return obj
end

"""
    normalize(obj::DenseInfiniteMPS)

Pure function of `normalize!`
"""
function normalize(obj::DenseInfiniteMPS)
    obj′ = deepcopy(obj)
    return normalize!(obj′)
end

# space
leftVirtualSpace(obj::DenseInfiniteMPS) = map(x->leftVirtualSpace(x), obj.AL)
rightVirtualSpace(obj::DenseInfiniteMPS) = map(x->rightVirtualSpace(x), obj.AL)
physicalSpace(obj::DenseInfiniteMPS) = map(x->physicalSpace(x), obj.AL)
extraPhysicalSpace(obj::DenseInfiniteMPS) = map(x->extraPhysicalSpace(x), obj.AL)

function Base.show(io::IO, obj::DenseInfiniteMPS{L}) where L
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

# Interface functions for `iterate`
function lambda(obj::DenseInfiniteMPS; which::Symbol)
    if which == :approximate || which == :ViTEBD
        λ_C = norm.(obj.C) .* sign.(obj.C)
        AC_data = map(x->x.A, obj.AC)
        λ_AC = sqrt.(tr.(adjoint.(AC_data) .* AC_data)) .* sign.(obj.C)
        return [λ_C, λ_AC]
    else
        throw(ArgumentError("Undefined behavior"))
    end
end
function division(obj::DenseInfiniteMPS, n::AbstractVector{<:Vector{<:Number}}; which::Symbol)
    if which == :approximate || which == :ViTEBD
        obj′ = deepcopy(obj)
        obj′.C ./= n[1]
        obj′.AC ./= n[2]
        return obj′
    else
        throw(ArgumentError("Undefined behavior"))
    end
end
function minus(obj1::T, obj2::T; which::Symbol) where T <: DenseInfiniteMPS
    if which == :approximate
        return obj1
    elseif which == :ViTEBD
        T <: DenseUniformMPS && return T(obj1.A-obj2.A)
        T <: DenseInfiniteMPS && return T(obj1.AL-obj2.AL, obj1.AR-obj2.AR, obj1.AC-obj2.AC, obj1.C-obj2.C)
    else
        throw(ArgumentError("Undefined behavior"))
    end
end
function norm_max(obj::DenseInfiniteMPS; which::Symbol)
    if which == :approximate
        return maximum(norm.(obj.AL .* obj.C .- obj.AC))
    elseif which == :ViTEBD
        return maximum(norm.(obj.AC))
    else
        throw(ArgumentError("Undefined behavior"))
    end
end
