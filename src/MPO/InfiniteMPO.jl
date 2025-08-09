"""
mutable struct InfiniteMPO{L, T} <: DenseInfiniteMPO{L}
    const A::AbstractVector{MPOTensor}
end

Concrete type of iMPO, where `L` is the cell size, `T == Float64` or `ComplexF64` is the number type of local tensors.

Graphic presentation:

           |            |
    ... -- A1 -- ... -- AL -- ...
           |            |

# Constructors
    InfiniteMPO{L, T}(::AbstractVector{<:MPOTensor})
    InfiniteMPO(::AbstractVector{<:MPOTensor})
    InfiniteMPO{L, T}()
    InfiniteMPO(L, T=Float64)
"""
mutable struct InfiniteMPO{L, T} <: DenseInfiniteMPO{L}
    const A::AbstractVector{MPOTensor}

    function InfiniteMPO{L, T}(A::AbstractVector{<:MPOTensor}) where {L, T}
        length(A) == L || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(A))"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types $T"))
        return new{L, T}(A)
    end
    function InfiniteMPO(A::AbstractVector{<:MPOTensor})
        L = length(A)
        T = mapreduce(eltype, promote_type, A)
        return new{L, T}(A)
    end

    function InfiniteMPO{L, T}() where {L, T}
        A = Vector{MPOTensor}(undef, L)
        return new{L, T}(A)
    end
    InfiniteMPO(L::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=InfiniteMPO{L, T}()
end

const iMPO = InfiniteMPO
