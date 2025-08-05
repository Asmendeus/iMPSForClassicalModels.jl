"""
    mutable struct InfiniteMPO{L, T} <: AbstractInfiniteMPO{L, 1}
        const O::AbstractVector{MPOTensor}
     end

Concrete type of iMPO, where `L` is the cell size, `T == Float64` or `ComplexF64` is the number type of local tensors.
"""

mutable struct InfiniteMPO{L, T} <: AbstractInfiniteMPO{L, 1}
    const O::AbstractVector{MPOTensor}

    function InfiniteMPO{L, T}() where {L, T}
        O = Vector{MPOTensor}(undef, L)
        return new{L, T}(O)
    end
    InfiniteMPO(L::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=InfiniteMPO{L, T}()

    function InfiniteMPO{L, T}(O::AbstractVector{<:MPOTensor}) where {L, T}
        length(O) == L || throw(ArgumentError("The InfiniteMPO length `L` does not match the length of the local tensor vector"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types"))
        return new{L, T}(O)
    end
    function InfiniteMPO(O::AbstractVector{<:MPOTensor})
        L = length(O)
        T = mapreduce(eltype, promote_type, O)
        return new{L, T}(O)
    end

end

const iMPO = InfiniteMPO

function Base.show(io::IO, mpo::InfiniteMPO{L, T}) where {L, T}
    # data type
    println(io, "InfiniteMPO{$L, $T}:")
    println(io)

    # schematic diagram: line 1
    print(io, repeat(" ", 3) * "|")
    for l in 1:L-1
        print(io, repeat(" ", 3 + length(string(l))) * "|")
    end
    println(io)

    # schematic diagram: line 2
    for l in 1:L
        print(io, " — O$l")
    end
    println(io, " —")

    # schematic diagram: line 3
    print(io, repeat(" ", 3) * "|")
    for l in 1:L-1
        print(io, repeat(" ", 3 + length(string(l))) * "|")
    end
    println(io)

    # local tensors
    println(io)
    println(io, "MPO tensors:")
    show(io, mpo.O)
end
