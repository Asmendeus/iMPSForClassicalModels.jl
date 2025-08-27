"""
    struct SparseGeneralMPO{W, L} <: AbstractGeneralInfiniteMPS{W, L}
        O::AbstractMatrix{<:AbstractMPOTensor}
    end

An iMPO type stores local tensors of the classical system's partition function.

# Note
Although it is written as "SparseGeneralMPO", it is not always sparse for a classical system, which is just a naming
convention inherited from `FiniteMPS.jl`. For a classical model with discrete degrees of freedom (like Ising, Clock),
the local tensor corresponding to the partition function is usually dense, with dense degree equal to about 1.
However, for a XY-like model with continuous degrees of freedom, the local tensors are usually sparse after discretization.
For example, dense degree of classical XY model with truncation dimension 5 (physical dimension d = 11) on square lattice
approximates to 0.06, which is typically sparse.

# Constructors
    SparseGeneralMPO{W, L}(A::AbstractMatrix{<:AbstractMPOTensor})
    SparseGeneralMPO(A::AbstractMatrix{<:AbstractMPOTensor})
    SparseGeneralMPO(A::AbstractMPOTensor)
"""
struct SparseGeneralMPO{W, L} <: AbstractGeneralInfiniteMPS{W, L}
    A::AbstractMatrix{<:AbstractMPOTensor}

    function SparseGeneralMPO{W, L}(A::AbstractMatrix{<:AbstractMPOTensor}) where {W, L}
        (W, L) == size(A) || throw(ArgumentError("Mismatched size: ($W, $L) ≠ $(size(A))"))
        return new{W, L}(A)
    end
    function SparseGeneralMPO(A::AbstractMatrix{<:AbstractMPOTensor})
        W, L = size(A)
        return new{W, L}(A)
    end
    function SparseGeneralMPO(A::AbstractMPOTensor)
        return new{1, 1}([A;;])
    end
end

convert(::Type{<:SparseGeneralMPO}, A::AbstractMatrix{<:AbstractMPOTensor}) = SparseGeneralMPO(A)
for func in (:getindex, :lastindex, :setindex!, :iterate, :keys, :isassigned)
    @eval Base.$func(obj::SparseGeneralMPO, args...) = $func(obj.A, args...)
end

"""
    getA(obj::SparseGeneralMPO)

Interface of `SparseGeneralMPO`, return the local tensors
"""
getA(obj::SparseGeneralMPO) = obj.A

function ishermitian(obj::SparseGeneralMPO{W, L}; tol::Float64=Defaults.tol_low) where {W, L}
    for w in 1:ceil(Int, W/2), l in 1:L
        norm(convert(Array, obj[w, l]) - convert(Array, permute(obj[end+1-w, l], (1, 3), (2, 4)))) < tol || return false
    end
    return true
end

function show(io::IO, obj::SparseGeneralMPO{W, L}) where {W, L}

    memory = obj |> Base.summarysize |> Base.format_bytes
    println(io, "$(typeof(obj)): total memory = $memory")
    # bond dimenson
    for w in 1:W
        println("line $w:")
        lsi = ceil(Int64, log10(L)) # length of si to be printed
        for si in 1:L
            D, DD = dim(obj[w, si], 1)
            print(io, "  Bond ", lpad(si-1, lsi), "->", lpad(si, lsi), ": ")
            println(io, "  $(sum(D)) -> $(sum(DD))")
        end
    end
end
