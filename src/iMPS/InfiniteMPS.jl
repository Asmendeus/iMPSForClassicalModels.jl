"""
    mutable struct InfiniteMPS{L} <: AbstractInfiniteMPS{L}
        const A::AbstractVector{AbstractTensorMap}
        const Centor::Vector{Int64}
     end

Concrete type of iMPS, where `L` is the cell size, `T == Float64` or `ComplexF64` is the number type of local tensors.
"""

mutable struct InfiniteMPS{L, T} <: AbstractInfiniteMPS{L}
    const A::AbstractVector{AbstractTensorMap}
    const Centor::Vector{Int64}

    function InfiniteMPS{L, T}() where {L, T}
        A = Vector{AbstractTensorMap}(undef, L)
        return new{L, T}(A, [1, L])
    end
    InfiniteMPS(L::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=InfiniteMPS{L, T}()

    function InfiniteMPS{L, T}(A::AbstractVector{<:AbstractTensorMap}, Centor::Vector{Int64}=[1, L]) where {L, T}
        length(A) == L || throw(ArgumentError("The InfiniteMPS length `L` does not match the length of the local tensor vector"))
        (length(Center) == 2 && Center[1] ≥ 1 && Center[2] ≤ L) || throw(ArgumentError("Illegal `Centor` vector"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types"))
        return new{L, T}(A, Centor)
    end
    function InfiniteMPS(A::AbstractVector{<:AbstractTensorMap}, Centor::Vector{Int64}=[1, length(A)])
        L = length(A)
        T = mapreduce(eltype, promote_type, A)
        return new{L, T}(A, Centor)
    end
end

const iMPS = InfiniteMPS

function Base.show(io::IO, mps::InfiniteMPS{L, T}) where {L, T}
    # data type
    println(io, "InfiniteMPS{$L, $T}:")
    println(io)

    # schematic diagram: line 1
    print(io, repeat(" ", 3) * "|")
    for l in 1:L-1
        print(io, repeat(" ", 3 + length(string(l))) * "|")
    end
    println(io)

    # schematic diagram: line 2
    for l in 1:L
        print(io, " — A$l")
    end
    println(io, " —")

    # local tensors
    println(io)
    println(io, "local tensors:")
    show(io, mps.A)
end
