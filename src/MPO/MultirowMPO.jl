"""
    mutable struct MultirowMPO{W, L, T} <: AbstractInfiniteMPS{L}
        const O::AbstractMatrix{AbstractTensorMap}
     end

Concrete type of iMPO, where `W` is the cell width and `L` is the cell length, `T == Float64` or `ComplexF64` is the number type of local tensors.
"""
mutable struct MultirowMPO{W, L, T} <: AbstractInfiniteMPS{L}
    const O::AbstractMatrix{AbstractTensorMap}

    function MultirowMPO{W, L, T}() where {W, L, T}
        O = Matrix{AbstractTensorMap}(undef, W, L)
        return new{W, L, T}(O)
    end
    MultirowMPO(L::Int, W::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=MultirowMPO{W, L, T}()

    function MultirowMPO{W, L, T}(O::AbstractMatrix{<:AbstractTensorMap}) where {W, L, T}
        size(A) == (W, L) || throw(ArgumentError("The MultirowMPO size `(W, L)` does not match the size of the local tensor matrix"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types"))
        return new{W, L, T}(O)
    end
    function MultirowMPO(O::AbstractMatrix{<:AbstractTensorMap})
        W, L = size(O)
        T = mapreduce(eltype, promote_type, O)
        return new{W, L, T}(O)
    end

end

const mMPO = MultirowMPO

function Base.show(io::IO, mpo::MultirowMPO{W, L, T}) where {W, L, T}
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
