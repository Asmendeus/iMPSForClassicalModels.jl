"""
    mutable struct InfiniteMPS{L} <: AbstractInfiniteMPS{L}
        const A::AbstractVector{MPSTEnsor}
        const Λ::AbstractVector{MPSBondTensor}
     end

Concrete type of iMPS, where `L` is the cell size, `T == Float64` or `ComplexF64` is the number type of local tensors.
"""

mutable struct InfiniteMPS{L, T} <: AbstractInfiniteMPS{L}
    const A::AbstractVector{MPSTensor}
    const Λ::AbstractVector{MPSBondTensor}

    function InfiniteMPS{L, T}() where {L, T}
        A = Vector{MPSTensor}(undef, L)
        Λ = Vector{MPSBondTensor}(undef, L)
        return new{L, T}(A, Λ)
    end
    InfiniteMPS(L::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=InfiniteMPS{L, T}()

    function InfiniteMPS{L, T}(A::AbstractVector{<:MPSTensor}, Λ::AbstractVector{<:MPSBondTensor}) where {L, T}
        length(A) == L || throw(ArgumentError("The InfiniteMPS length `L` does not match the number of `MPSTensor` in `A`"))
        length(Λ) == L || throw(ArgumentError("The InfiniteMPS length `L` does not match the number of `MPSBondTensor` in `Λ`"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types $T"))
        return new{L, T}(A, Λ)
    end
    function InfiniteMPS(A::AbstractVector{<:MPSTensor}, Λ::AbstractVector{<:MPSBondTensor})
        L = length(A)
        length(Λ) == L || throw(ArgumentError("The number of `MPSTensor` in `A` does not match the number of `MPSBondTensor` in `Λ`"))
        T = promote_type(mapreduce(eltype, promote_type, A), mapreduce(eltype, promote_type, Λ))
        return new{L, T}(A, Λ)
    end
    # TODO
    # function InfiniteMPS(A::AbstractVector{<:MPSTensor})
    #     L = length(A)
    #     T = mapreduce(eltype, promote_type, A)
    #     Λ = let
    #         # 每个 bond tensor 根据左右 bond dimension 生成 identity matrix
    #     end
    #     return new{L, T}(A, Λ)
    # end
end

const iMPS = InfiniteMPS

function Base.show(io::IO, mps::InfiniteMPS{L, T}) where {L, T}
    # data type
    println(io, "InfiniteMPS{$L, $T}:")
    println(io)

    # schematic diagram: line 1
    print(io, repeat(" ", 3) * "|")
    for l in 1:L-1
        print(io, repeat(" ", 6 + 2 * length(string(l))) * "|")
    end
    println(io)

    # schematic diagram: line 2
    for l in 1:L
        print(io, " — A$l — Λ$l")
    end
    println(io, " —")

    # local tensors
    println(io)
    println(io, "MPS tensors:")
    show(io, mps.A)
    println(io)

    println(io)
    println(io, "MPS bond tensors:")
    show(io, mps.Λ)
    println(io)
end
