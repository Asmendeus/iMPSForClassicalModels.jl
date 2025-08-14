"""
    abstract type AbstractInfiniteMPS{L}

Abstract type of all iMPS/iMPO with unit cell length `L`, marking the size of cell.
"""
abstract type AbstractInfiniteMPS{L} end

length(::AbstractInfiniteMPS{L}) where L = L
for func in (:getindex, :lastindex, :setindex!, :iterate, :keys, :isassigned)
    @eval Base.$func(obj::AbstractInfiniteMPS, args...) = $func(obj.A, args...)
end

"""
    abstract type DenseInfiniteMPS{L, T<:Union{Float64, ComplexF64}} <: AbstractInfiniteMPS{L} end

Abstract type of dense iMPS/iMPO with unit cell length `L`.
"""
abstract type DenseInfiniteMPS{L, T<:Union{Float64, ComplexF64}} <: AbstractInfiniteMPS{L} end

# promote local tensors
scalartype(::DenseInfiniteMPS{L, T}) where {L, T} = T
function Base.setindex!(obj::DenseInfiniteMPS{L, ComplexF64}, A::T, si::Int64) where {L, T <:Union{AbstractTensorMap, AbstractLocalTensor}}
    if scalartype(A) != ComplexF64
        return setindex!(obj.A, A*one(ComplexF64), si)
    else
        return setindex!(obj.A, A, si)
    end
end

"""
     normalize!(obj::DenseInfiniteMPS) -> obj

Normalize a given iMPS/iMPO according to inner-induced norm.

Note we assume the iMPS/iMPO satisfies a canonical form and the center tensor is normalized, hence we only normalize `c`.
"""
function normalize!(obj::DenseInfiniteMPS)
    obj.c /= norm(obj)
    return obj
end

"""
     norm(obj::DenseInfiniteMPS) -> ::Float64

Return the inner-induced norm. Note we assume the iMPS/iMPO satisfies a canonical form and the center tensor is normalized, hence the norm is just `abs(c)`.
"""
function norm(obj::DenseInfiniteMPS)
    return abs(obj.c)
end

"""
     coef(obj::DenseInfiniteMPS) -> ::F

Interface of `DenseInfiniteMPS`, return the global coefficient, where `F` is the number type of given MPS.
"""
coef(obj::DenseInfiniteMPS) = obj.c

"""
     Center(obj::DenseInfiniteMPS) -> Vector (length 2)

Interface of `DenseInfiniteMPS`, return the info of canonical center. `[a, b]` means left-canonical from `1` to `a-1` and right-canonical from `b+1` to `L`.
"""
Center(obj::DenseInfiniteMPS) = obj.Center

"""
     complex(obj::DenseInfiniteMPS{L}) -> ::DenseInfiniteMPS{L, ComplexF64}

Return a copy of given iMPS/iMPO but with `ComplexF64` as basic field.
"""
function complex(obj::DenseInfiniteMPS{L, Float64}) where L
    obj_c = similar(ComplexF64, obj)
    obj_c.c = obj.c
    obj_c.Center[:] = obj.Center[:]
    for i = 1:L
        obj_c[i] = obj[i]
    end
    return obj_c
end
complex(obj::DenseInfiniteMPS{L, ComplexF64}) where L = deepcopy(obj)

function Base.show(io::IO, obj::DenseInfiniteMPS{L}) where L
    any(i -> !isassigned(obj.A, i), 1:L) && return println(io, "$(typeof(obj)): L = $L, to be initialized !")

    # avoid to show bond info
    memory = obj |> Base.summarysize |> Base.format_bytes
    println(io, "$(typeof(obj)): Center = $(Center(obj)), Norm = $(norm(obj)), total memory = $memory")
    # bond dimenson
    lsi = ceil(Int64, log10(L)) # length of si to be printed
    for si = 1:L
        local A = obj[si]
        D, DD = dim(A, 1)
        println(io, "Bond ", lpad(si-1, lsi), "->", lpad(si, lsi), ": $(codomain(A).spaces[1]), dim = $(D) -> $(DD)")
    end
end
