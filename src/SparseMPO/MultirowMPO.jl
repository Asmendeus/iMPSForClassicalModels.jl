"""
    mutable struct MultirowMPO{L, W, T} <: AbstractInfiniteMPS{L}
        const O::AbstractMatrix{AbstractTensorMap}
     end

Concrete type of iMPO, where `L` is the cell length and `W` is the cell width, `T == Float64` or `ComplexF64` is the number type of local tensors.
"""
mutable struct MultirowMPO{L, W, T} <: AbstractInfiniteMPS{L}
    const O::AbstractMatrix{AbstractTensorMap}

    function MultirowMPO{L, W, T}() where {L, W, T}
        O = Matrix{AbstractTensorMap}(undef, L, W)
        return new{L, W, T}(O)
    end
    MultirowMPO(L::Int, W::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=MultirowMPO{L, W, T}()

    function MultirowMPO{L, W, T}(O::AbstractMatrix{<:AbstractTensorMap}) where {L, W, T}
        size(A) == (L, W) || throw(ArgumentError("The MultirowMPO size `(L, W)` does not match the size of the local tensor matrix"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types"))
        return new{L, W, T}(O)
    end
    function MultirowMPO(O::AbstractMatrix{<:AbstractTensorMap})
        L, W = size(O)
        T = mapreduce(eltype, promote_type, O)
        return new{L, W, T}(O)
    end

end

const miMPO = MultirowMPO

function MultirowMPO(mpo::InfiniteMPO{L, T}) where {L, T}
    return MultirowMPO{L, 1, T}([e for _ in 1:1, e in mpo.O])
end
function InfiniteMPO(mpo::MultirowMPO{L, 1, T}) where {L, T}
    return InfiniteMPO{L, T}(vec(mpo.O))
end

function Base.show(io::IO, mpo::MultirowMPO{L, W, T}) where {L, W, T}
    # data type
    println(io, "MultirowMPO{$L, $W, $T}:")
    println(io)

    # graphic presentation: line 1
    print(io, repeat(" ", 3) * "|")
    for l in 1:L-1
        print(io, repeat(" ", 3 + length(string(l))) * "|")
    end
    println(io)

    for w in 1:W
        # graphic presentation: line 2
        for l in 1:L
            print(io, " — $(Char(78+w))$l")
        end
        println(io, " —")

        # graphic presentation: line 3
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
