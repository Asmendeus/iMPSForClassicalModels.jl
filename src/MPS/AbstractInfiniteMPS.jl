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

Abstract type of dense iMPS/iMPO with unit cell length `L` and data type `T`.
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
    Center(obj::DenseInfiniteMPS) -> ::Union{Nothing, Int64}

Interface of `DenseInfiniteMPS`, return the info of canonical center.
"""
Center(obj::DenseInfiniteMPS) = obj.Center

"""
    iscanonical(obj::DenseInfiniteMPS) -> ::Bool

Return !isnothing(obj.Center)
"""
iscanonical(obj::DenseInfiniteMPS) = !isnothing(obj.Center)

"""
    isuniform(obj::DenseInfiniteMPS) -> ::Bool

Return isnothing(obj.Center)
"""
isuniform(obj::DenseInfiniteMPS) = isnothing(obj.Center)

"""
    complex(obj::DenseInfiniteMPS{L, Float64}) -> ::DenseInfiniteMPS{L, ComplexF64}

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

"""
    getAllCanonicalFormTensors(obj::DenseInfiniteMPS{L, T}; kwargs...)
        -> AL::Vector{LocalTensor{R}}, AR::Vector{LocalTensor{R}}, AC::Vector{LocalTensor{R}}, C::Vector{BondTensor}

Get all tensors of an iMPS/iMPO with canonical form.

`kwargs` is propagated to `setCenter`, `leftorth` and `rightorth`
"""
function getAllCanonicalFormTensors(obj::DenseInfiniteMPS{L, T}; kwargs...) where {L, T}
    A = setCenter(obj, 1; kwargs...).A

    AL = typeof(A)(undef, L)
    AR = typeof(A)(undef, L)
    AC = typeof(A)(undef, L)
    C = Vector{BondTensor}(undef, L)

    AC[1] = A[1]
    for l in 1:L-1
        A[l], C[l+1], x = leftorth(A[l]; ismerge=false, kwargs...)
        A[l+1] = C[l+1] * x * A[l+1]
        AL[l] = A[l]
        AC[l+1] = A[l+1]
    end
    AL[L], C[1], x = leftorth(A[L]; ismerge=false, kwargs...)

    AC = map(l -> AL[l] * C[mod(l, L)+1], 1:L)
    AR = map(l -> inv(C[l]) * AC[l], 1:L)

    @assert norm(C[1] * x * AR[1] - AC[1]) < tol
    # for l in 1:L
    #     l₊ = mod(l, L) + 1
    #     @assert norm(AL[l] * C[l₊] - AC[l]) < Defaults.tol "Mismatched results of left- and right-canonicalization"
    #     @assert norm(C[l] * AR[l] - AC[l]) < Defaults.tol "Mismatched results of left- and right-canonicalization"
    # end

    return AL, AR, AC, C
end

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
