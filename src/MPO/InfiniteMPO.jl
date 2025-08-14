"""
    mutable struct InfiniteMPO{L, T<:Union{Float64, ComplexF64}} <: DenseInfiniteMPS{L, T}
        const A::AbstractVector{<:AbstractMPOTensor}
        const Center::Vector{Int64}
        c::T
    end

Concrete type of iMPO, where `L` is the cell length, `T == Float64` or `ComplexF64` is the number type of local tensors.
All the fields and constructors are exactly the same to those of `InfiniteMPS`, we redefine the type `InfiniteMPO` just for using multiple dispacth when implementing the algebra between `InfiniteMPS` and `InfiniteMPO`.

Details of constructors please see `InfiniteMPS`.
"""
mutable struct InfiniteMPO{L, T<:Union{Float64, ComplexF64}} <: DenseInfiniteMPS{L, T}
    const A::AbstractVector{<:AbstractMPOTensor}
    const Center::Vector{Int64}
    c::T

    function InfiniteMPO{L, T}(A::AbstractVector{<:AbstractMPOTensor}=Vector{MPOTensor}(undef, L), Center::Vector{Int64}=[1, L], c::T=one(T)) where {L, T <: Union{Float64, ComplexF64}}
        L == length(A) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(A))"))
        length(Center) == 2 || throw(ArgumentError("Illegal length of `Center` $(length(Center)), which should be 2"))
        1 ≤ Center[1] && Center[2] ≤ L || throw(ArgumentError("Orthogonal boundaries beyond the range of Imps"))

        if T == ComplexF64  # promote each A
            try
                for i = 1:L
                    eltype(A[i]) != T && (A[i] *= one(T))
                end
            catch
                nothing
            end
        end

        if !(eltype(A) <: AbstractMPOTensor)
            A = convert(Vector{MPOTensor}, A)
        end

        return new{L, T}(A, Center, c)
    end
    InfiniteMPO(L::Int64, T::DataType=Defaults.datatype) = InfiniteMPO{L, T}()

    function InfiniteMPO(A::AbstractVector{<:AbstractMPOTensor}, Center::Vector{Int64}=[1, length(A)], c::T=one(mapreduce(eltype, promote_type, A))) where T
        L = length(A)
        T == mapreduce(eltype, promote_type, A) || throw(ArgumentError("Mismatched datatypes: $(mapreduce(eltype, promote_type, A)) ≠ $T"))
        return InfiniteMPO{L, T}(A, Center, c)
    end
    function InfiniteMPO(A::AbstractVector{<:AbstractTensorMap}, Center::Vector{Int64}=[1, length(A)], c::T=one(mapreduce(eltype, promote_type, A))) where T
        return InfiniteMPO(convert(Vector{MPOTensor}, A), Center, c)
    end
end

const iMPO = InfiniteMPO

"""
    identityInfiniteMPO(::Type{T}, pspace::AbstractVector{<:VectorSpace})

Construct an identity iMPO where the physical spaces are informed by a length `L` vertor `pspace`.

    identityInfiniteMPO(::Type{T}, L::Int64, pspace::VectorSpace)
Assume the all the physical spaces are the same.

    identityInfiniteMPO(obj::DenseInfiniteMPS{L, T})
Deduce the scalar type `T` and physical spaces from an iMPS/iMPO.
"""
function identityInfiniteMPO(::Type{T}, pspace::AbstractVector{<:VectorSpace}) where T <: Union{Float64, ComplexF64}

    L = length(pspace)

    obj = InfiniteMPO(L, T)
    aspace = trivial(pspace[1])

    for si in 1:L
        Ip = normalize(id(pspace[si]))
        Ia = id(aspace)
        @tensor tmp[-1 -2; -3 -4] := Ia[-1 -4] * Ip[-2 -3]
         obj[si] = tmp
    end
    obj.Center[:] = [1, 1]

    return obj
end
function identityInfiniteMPO(::Type{T}, L::Int64, pspace::VectorSpace) where T <: Union{Float64, ComplexF64}
    aspace = trivial(pspace)
    Ip = normalize(id(pspace))
    Ia = id(aspace)
    @tensor tmp[-1 -2; -3 -4] := Ia[-1 -4] * Ip[-2 -3]
    return InfiniteMPO{L, T}(repeat([MPOTensor(tmp),], L), [1, 1], one(T))
end
function identityInfiniteMPO(obj::DenseInfiniteMPS{L, T}) where {L, T}
    pspace = map(i->codomain(obj.A[i], 2), 1:L)
    return identityInfiniteMPO(T, pspace)
end
