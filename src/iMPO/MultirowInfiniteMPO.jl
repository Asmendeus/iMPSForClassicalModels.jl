"""
    mutable struct MultirowInfiniteMPO{L, W, T} <: AbstractInfiniteMPO{L, W}
        const O::AbstractMatrix{AbstractTensorMap}
     end

Concrete type of iMPO, where `L` is the cell length and `W` is the cell width, `T == Float64` or `ComplexF64` is the number type of local tensors.
"""
mutable struct MultirowInfiniteMPO{L, W, T} <: AbstractInfiniteMPO{L, W}
    const O::AbstractMatrix{AbstractTensorMap}

    function MultirowInfiniteMPO{L, W, T}() where {L, W, T}
        O = Matrix{AbstractTensorMap}(undef, L, W)
        return new{L, W, T}(O)
    end
    MultirowInfiniteMPO(L::Int, W::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=MultirowInfiniteMPO{L, W, T}()

    function MultirowInfiniteMPO{L, W, T}(O::AbstractMatrix{<:AbstractTensorMap}) where {L, W, T}
        size(A) == (L, W) || throw(ArgumentError("The MultirowInfiniteMPO size `(L, W)` does not match the size of the local tensor matrix"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types"))
        return new{L, W, T}(O)
    end
    function MultirowInfiniteMPO(O::AbstractMatrix{<:AbstractTensorMap})
        L, W = size(O)
        T = mapreduce(eltype, promote_type, O)
        return new{L, W, T}(O)
    end

end

const miMPO = MultirowInfiniteMPO

function MultirowInfiniteMPO(mpo::InfiniteMPO{L, T}) where {L, T}
    return MultirowInfiniteMPO{L, 1, T}([e for _ in 1:1, e in mpo.O])
end
function InfiniteMPO(mpo::MultirowInfiniteMPO{L, 1, T}) where {L, T}
    return InfiniteMPO{L, T}(vec(mpo.O))
end

function Base.show(io::IO, mpo::MultirowInfiniteMPO{L, W, T}) where {L, W, T}
    # data type
    println(io, "MultirowInfiniteMPO{$L, $W, $T}:")
    println(io)

    # schematic diagram: line 1
    print(io, repeat(" ", 3) * "|")
    for l in 1:L-1
        print(io, repeat(" ", 3 + length(string(l))) * "|")
    end
    println(io)

    for w in 1:W
        # schematic diagram: line 2
        for l in 1:L
            print(io, " — $(Char(78+w))$l")
        end
        println(io, " —")

        # schematic diagram: line 3
        print(io, repeat(" ", 3) * "|")
        for l in 1:L-1
            print(io, repeat(" ", 3 + length(string(l))) * "|")
        end
        println(io)
    end

    # local tensors
    println(io)
    println(io, "local tensors:")
    show(io, mpo.O)
end
