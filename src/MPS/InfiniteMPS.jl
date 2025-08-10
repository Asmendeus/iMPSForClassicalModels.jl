"""
    mutable struct InfiniteMPS{L, T} <: DenseInfiniteMPS{L}
        const A::AbstractVector{MPSTensor}
    end

Concrete type of iMPS, where `L` is the cell size, `T == Float64` or `ComplexF64` is the number type of local tensors.

Graphic presentation:

    ... -- A1 -- ... -- AL -- ...
           |            |

# Constructors
    InfiniteMPS{L, T}(::AbstractVector{<:MPSTensor})
    InfiniteMPS(::AbstractVector{<:MPSTensor})
    InfiniteMPS{L, T}()
    InfiniteMPS(L, T=Float64)
"""

mutable struct InfiniteMPS{L, T} <: DenseInfiniteMPS{L}
    const A::AbstractVector{MPSTensor}

    function InfiniteMPS{L, T}(A::AbstractVector{<:MPSTensor}) where {L, T}
        length(A) == L || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(A))"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types $T"))
        return new{L, T}(A)
    end
    function InfiniteMPS(A::AbstractVector{<:MPSTensor})
        L = length(A)
        T = mapreduce(eltype, promote_type, A)
        return new{L, T}(A)
    end

    function InfiniteMPS{L, T}() where {L, T}
        A = Vector{MPSTensor}(undef, L)
        return new{L, T}(A)
    end
    InfiniteMPS(L::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=InfiniteMPS{L, T}()
end

const iMPS = InfiniteMPS

"""
    randInfiniteMPS([::Type{T},], pspace::AbstractVector{<:VectorSpace}, aspace::AbstractVector{<:VectorSpace})::InfiniteMPSMPS{L, T}

Generate a cell length `L` random infinite MPS with given length `L` vector `pspace` and vector `aspace`. `T = Float64`(default) or `ComplexF64` is the number type.
Notice that 1-st and (L+1)-th auxiliary space are the same, so we only take the 1 ~ L to make up the `aspace`.

    randInfiniteMPS([::Type{T},], d::AbstractVector{Int64}, D::AbstractVector{Int64})::InfiniteMPSMPS{L, T}

Assume no symmetry is added. `T = Float64` corresponds to `CartesianSpace`, while `T = ComplexF64` corresponds to `ComplexSpace`.

    randInfiniteMPS([::Type{T},], L::Int64, pspace::VectorSpace, apsace::VectorSpace)::InfiniteMPSMPS{L, T}

Assume the same `pspace` and `aspace`.

    randInfiniteMPS([::Type{T},], L::Int64, d::Int64, D::Int64)::InfiniteMPSMPS{L, T}

Assume the same `pspace` with physical dimension `d` and `aspace` with bond dimension `D`.
"""
function randInfiniteMPS(::Type{T}, pspace::AbstractVector{<:VectorSpace}, aspace::AbstractVector{<:VectorSpace}) where T <: Union{Float64, ComplexF64}
    (L = length(pspace) == length(aspace)) || throw(ArgumentError("Mismatched lengths of pspace and aspace: $(length(pspace)) ≠ $(length(aspace))"))
end
function randInfiniteMPS(::Type{T}, d::AbstractVector{Int64}, D::AbstractVector{Int64}) where T <: Union{Float64, ComplexF64}
    (L = length(d) == length(D)) || throw(ArgumentError("Mismatched lengths of pspace and aspace: $(length(d)) ≠ $(length(D))"))
    space = T == Float64 ? ℝ : ℂ
    return randInfiniteMPS(T, map(l->space^d[l], 1:L), map(l->space^D[l], 1:L))
end
function randInfiniteMPS(::Type{T}, L::Int64, pspace::VectorSpace, apsace::VectorSpace) where T <: Union{Float64, ComplexF64}
    return randInfiniteMPS(T, repeat([pspace,], L), repeat([aspace,], L))
end
function randInfiniteMPS(::Type{T}, L::Int64, d::Int64, D::Int64) where T <: Union{Float64, ComplexF64}
    return randInfiniteMPS(T, repeat([d,], L), repeat([D,], L))
end
