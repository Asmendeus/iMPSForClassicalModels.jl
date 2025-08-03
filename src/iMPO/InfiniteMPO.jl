"""
    mutable struct InfiniteMPO{L, T} <: AbstractInfiniteMPO{L, 1}
        const A::AbstractVector{AbstractTensorMap}
     end

Concrete type of iMPO, where `L` is the cell size, `T == Float64` or `ComplexF64` is the number type of local tensors.
"""

mutable struct InfiniteMPO{L, T} <: AbstractInfiniteMPO{L, 1}
    const A::AbstractVector{AbstractTensorMap}

    function InfiniteMPO{L, T}() where {L, T}
        A = Vector{AbstractTensorMap}(undef, L)
        return new{L, T}(A)
    end
    InfiniteMPO(L::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=InfiniteMPO{L, T}()

    function InfiniteMPO{L, T}(A::AbstractVector{<:AbstractTensorMap}) where {L, T}
        length(A) == L || throw(ArgumentError("The InfiniteMPO length `L` does not match the length of the local tensor vector"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types"))
        return new{L, T}(A)
    end
    function InfiniteMPO(A::AbstractVector{<:AbstractTensorMap})
        L = length(A)
        T = mapreduce(eltype, promote_type, A)
        return new{L, T}(A)
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
        print(io, " — A$l")
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
    println(io, "local tensors:")
    show(io, mpo.A)
end
