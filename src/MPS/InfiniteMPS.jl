"""
     mutable struct InfiniteMPS{L, T<:Union{Float64, ComplexF64}} <: DenseInfiniteMPS{L, T}
          const A::AbstractVector{<:AbstractMPSTensor}
          const Center::Vector{Int64}
          c::T
     end

Concrete type of iMPS, where `L` is the cell length, `T == Float64` or `ComplexF64` is the number type of local tensors.

# Fields
    const A::AbstractVector{<:AbstractMPSTensor}
Length L vector to store the local tensors. Note the vector `A` is immutable while the local tensors in it are mutable.

    const Center::Vector{Int64}
Length 2 vector to label the canonical form. `[a, b]` means left-canonical from `1` to `a-1` and right-canonical from `b+1` to `L`.

    c::T
The global coefficient, i.e. we represented an iMPS with L local tensors and an additional scalar `c`, in order to avoid too large/small local tensors.

# Constructors
    InfiniteMPS{L, T}(A::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L), Center::Vector{Int64}=[1, L], c::T = one(T))
Standard constructor. Note we will promote all local tensors if `T == ComplexF64` while the input local tensors in `A` are of `Float64`.

    InfiniteMPS(L::Int64, T::DataType=Defaults.datatype)
Initialize an iMPS{L, T} with undef local tensors `A` to be filled. Note we initialize `Center = [1, L]` and `c = one(T)`.

    InfiniteMPS(A::AbstractVector{<:AbstractMPSTensor}, Center::Vector{Int64}=[1, L], c::T=one(T))
Return InfiniteMPS{length(A), T}(A, Center, c)

    InfiniteMPS(A::AbstractVector{<:AbstractTensorMap}, Center::Vector{Int64}=[1, L], c::T=one(T))
Automatically convert local tensors to warpper type `MPSTensor`.
"""
mutable struct InfiniteMPS{L, T<:Union{Float64, ComplexF64}} <: DenseInfiniteMPS{L, T}
    const A::AbstractVector{<:AbstractMPSTensor}
    const Center::Vector{Int64}
    c::T

    function InfiniteMPS{L, T}(A::AbstractVector{<:AbstractMPSTensor}=Vector{MPSTensor}(undef, L), Center::Vector{Int64}=[1, L], c::T=one(T)) where {L, T}
        L == length(A) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(A))"))
        length(Center) == 2 || throw(ArgumentError("Illegal length of `Center` $(length(Center)), which should be 2"))
        1 ≤ Center[1] && Center[2] ≤ L || throw(ArgumentError("Orthogonal boundaries beyond the range of Imps"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data type $T, which should be `Float64` or `ComplexF64`"))

        if T == ComplexF64  # promote each A
            for i = 1:L
                 eltype(A[i]) != T && (A[i] *= one(T))
            end
        end

        if !(eltype(A) <: AbstractMPSTensor)
            A = convert(Vector{MPSTensor}, A)
        end

        return new{L, T}(A, Center, c)
    end
    InfiniteMPS(L::Int64, T::DataType=Defaults.datatype) = InfiniteMPS{L, T}()

    function InfiniteMPS(A::AbstractVector{<:AbstractMPSTensor}, Center::Vector{Int64}=[1, length(A)], c::T=one(mapreduce(eltype, promote_type, A))) where T
        L = length(A)
        T == mapreduce(eltype, promote_type, A) || throw(ArgumentError("Mismatched datatypes: $(mapreduce(eltype, promote_type, A)) ≠ $T"))
        return InfiniteMPS{L, T}(A, Center, c)
    end
    function InfiniteMPS(A::AbstractVector{<:AbstractTensorMap}, Center::Vector{Int64}=[1, length(A)], c::T=one(mapreduce(eltype, promote_type, A))) where T
        return InfiniteMPS(convert(Vector{MPSTensor}, A), Center, c)
    end
end

const iMPS = InfiniteMPS

"""
     randInfiniteMPS([::Type{T},] pspace::Vector{VectorSpace}, aspace::Vector{VectorSpace}; kwargs...) -> InfiniteMPS{L, T}

Generate a length `L` random MPS with given length `L` vector `pspace` and `aspace`. `T = Float64`(default) or `ComplexF64` is the number type. Note the canonical center is initialized to the first site.

     randInfiniteMPS([::Type{T},] L::Int64, pspace::VectorSpace, apsace::VectorSpace; kwargs...) -> InfiniteMPS{L, T}

Assume the same `pspace` and `aspace`, except for the boundary bond, which is assumed to be trivial.

    randInfiniteMPS([::Type{T},] pdim::Vector{Int64}, adim::Vector{Int64}; kwargs...) -> InfiniteMPS{L, T}

Assume all spaces are trivial space, i.e., ℝ for `T == Float64` and ℂ for `T == ComplexF64`.

    randInfiniteMPS([::Type{T},] L::Int64, pdim::Int64, adim::Int64; kwargs...)

Assume the same trivial `pspace` and `aspace`.
"""
