"""
    mutable struct MultirowInfiniteMPO{W, L, T} <: AbstractInfiniteMPO{W, L}
        const O::AbstractMatrix{AbstractTensorMap}
     end

Concrete type of iMPO, where `L` is the cell length and `W` is the cell width, `T == Float64` or `ComplexF64` is the number type of local tensors.
"""
mutable struct MultirowInfiniteMPO{W, L, T} <: AbstractInfiniteMPO{W, L}
    const O::AbstractMatrix{AbstractTensorMap}

    function MultirowInfiniteMPO{W, L, T}() where {W, L, T}
        O = Matrix{AbstractTensorMap}(undef, W, L)
        return new{W, L, T}(O)
    end
    MultirowInfiniteMPO(W::Int, L::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=MultirowInfiniteMPO{W, L, T}()

    function MultirowInfiniteMPO{W, L, T}(O::AbstractMatrix{<:AbstractTensorMap}) where {W, L, T}
        size(A) == (W, L) || throw(ArgumentError("The MultirowInfiniteMPO size `(W, L)` does not match the size of the local tensor matrix"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types"))
        return new{W, L, T}(O)
    end
    function MultirowInfiniteMPO(O::AbstractMatrix{<:AbstractTensorMap})
        W, L = size(O)
        T = mapreduce(eltype, promote_type, O)
        return new{W, L, T}(O)
    end

end

const miMPO = MultirowInfiniteMPO

function MultirowInfiniteMPO(mpo::InfiniteMPO{L, T}) where {L, T}
    return MultirowInfiniteMPO{1, L, T}([e for _ in 1:1, e in mpo.O])
end
function InfiniteMPO(mpo::MultirowInfiniteMPO{1, L, T}) where {L, T}
    return InfiniteMPO{L, T}(vec(mpo.O))
end

function Base.show(io::IO, mpo::MultirowInfiniteMPO{W, L, T}) where {W, L, T}
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
