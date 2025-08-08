"""
    mutable struct InfiniteMPS{L} <: AbstractInfiniteMPS{L}
        const A::AbstractVector{MPSTensor}
     end

Concrete type of iMPS, where `L` is the cell size, `T == Float64` or `ComplexF64` is the number type of local tensors.
"""

mutable struct InfiniteMPS{L, T} <: DenseInfiniteMPS{L}
    const A::AbstractVector{MPSTensor}

    function InfiniteMPS{L, T}() where {L, T}
        A = Vector{MPSTensor}(undef, L)
        return new{L, T}(A)
    end
    InfiniteMPS(L::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=InfiniteMPS{L, T}()

    function InfiniteMPS{L, T}(A::AbstractVector{<:MPSTensor}) where {L, T}
        length(A) == L || throw(ArgumentError("The InfiniteMPS length `L` does not match the number of `MPSTensor` in `A`"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types $T"))
        return new{L, T}(A)
    end
    function InfiniteMPS(A::AbstractVector{<:MPSTensor})
        L = length(A)
        T = mapreduce(eltype, promote_type, A)
        return new{L, T}(A)
    end
end

const iMPS = InfiniteMPS

function Base.show(io::IO, mps::InfiniteMPS{L, T}) where {L, T}
    # data type
    println(io, "InfiniteMPS{$L, $T}:")
    println(io)

    # graphic presentation: line 1
    print(io, repeat(" ", 3) * "|")
    for l in 1:L-1
        print(io, repeat(" ", 3 + length(string(l))) * "|")
    end
    println(io)

    # graphic presentation: line 2
    for l in 1:L
        print(io, " — A$l")
    end
    println(io, " —")
    println(io)

    # local tensors
    println(io, "MPS tensors:")
    show(io, mps.A)
    println(io)
end
