"""
    mutable struct MultirowRightEnvironment{L, N} <: AbstractEnvironment{L, N}
        const AR::AbstractVector{RightIsometricTensor}
        const O::AbstractMatrix{MPOTensor}
        const BR::AbstractVector{AdjointRightIsometricTensor}
    end

Wrapper type for Right environment tensor's generating environment.

Graphic presentation:

    —— BR1 —— ... —— BRL ——
       |             |
    —— P1  —— ... —— PL  ——    (N-2)-th
       ⋮              ⋮
    —— O1  —— ... —— OL  ——    1-st
       |             |
    —— AR1 —— ... —— ARL ——

# Constructors
    MultirowRightEnvironment{L, N}(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
    MultirowRightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
    MultirowRightEnvironment(AR::RightIsometricTensor, O::AbstractMatrix{<:MPOTensor}, BR::AdjointRightIsometricTensor)
"""
mutable struct MultirowRightEnvironment{L, N} <: AbstractEnvironment{L, N}
    const AR::AbstractVector{RightIsometricTensor}
    const O::AbstractMatrix{MPOTensor}
    const BR::AbstractVector{AdjointRightIsometricTensor}

    MultirowRightEnvironment{L, N}(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}) where {L, N} = new{L, N}(AR, O, BR)
    function MultirowRightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
        L = length(AR)
        N = size(O, 2) + 2
        (L == size(O, 1) == length(BR)) || throw(ArgumentError("The lengths of `AR`, `O` and `BR` do not match"))
        return new{L, N}(AR, O, BR)
    end
    function MultirowRightEnvironment(AR::RightIsometricTensor, O::AbstractMatrix{<:MPOTensor}, BR::AdjointRightIsometricTensor)
        return new{1, N}([AR,], O, [BR,])
    end
end
