"""
    mutable struct RightEnvironment{L} <: AbstractEnvironment{L, 3}
        const AR::AbstractVector{RightIsometricTensor}
        const O::AbstractVector{MPOTensor}
        const BR::AbstractVector{AdjointRightIsometricTensor}
    end

Wrapper type for Right environment tensor's generating environment.

Graphic presentation:

    —— BR1 —— ... —— BRL ——
       |             |
    —— O1  —— ... —— OL  ——
       |             |
    —— AR1 —— ... —— ARL ——

# Constructors
    RightEnvironment{L}(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
    RightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
    RightEnvironment(AR::RightIsometricTensor, O::MPOTensor, BR::AdjointRightIsometricTensor)
"""
mutable struct RightEnvironment{L} <: AbstractEnvironment{L, 3}
    const AR::AbstractVector{RightIsometricTensor}
    const O::AbstractVector{MPOTensor}
    const BR::AbstractVector{AdjointRightIsometricTensor}

    RightEnvironment{L}(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}) where L = new{L}(AR, O, BR)
    function RightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
        L = length(AR)
        (L == length(O) == length(BR)) || throw(ArgumentError("The lengths of `AR`, `O` and `BR` do not match"))
        return new{L}(AR, O, BR)
    end
    function RightEnvironment(AR::RightIsometricTensor, O::MPOTensor, BR::AdjointRightIsometricTensor)
        return new{1}([AR,], [O,], [BR,])
    end
end
