"""
    mutable struct MultirowLeftEnvironment{L, N} <: AbstractEnvironment{L, N}
        const AL::AbstractVector{LeftIsometricTensor}
        const O::AbstractMatrix{MPOTensor}
        const BL::AbstractVector{AdjointLeftIsometricTensor}
    end

Wrapper type for left environment tensor's generating environment.

Graphic presentation:

    —— BL1 —— ... —— BLL ——
       |             |
    —— P1  —— ... —— PL  ——    (N-2)-th
       ⋮              ⋮
    —— O1  —— ... —— OL  ——    1-st
       |             |
    —— AL1 —— ... —— ALL ——

# Constructors
    MultirowLeftEnvironment{L, N}(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
    MultirowLeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
    MultirowLeftEnvironment(AL::LeftIsometricTensor, O::AbstractMatrix{<:MPOTensor}, BL::AdjointLeftIsometricTensor)
"""
mutable struct MultirowLeftEnvironment{L, N} <: AbstractEnvironment{L, N}
    const AL::AbstractVector{LeftIsometricTensor}
    const O::AbstractMatrix{MPOTensor}
    const BL::AbstractVector{AdjointLeftIsometricTensor}

    MultirowLeftEnvironment{L, N}(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}) where {L, N} = new{L, N}(AL, O, BL)
    function MultirowLeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
        L = length(AL)
        N = size(O, 2) + 2
        (L == size(O, 1) == length(BL)) || throw(ArgumentError("The lengths of `AL`, `O` and `BL` do not match"))
        return new{L, N}(AL, O, BL)
    end
    function MultirowLeftEnvironment(AL::LeftIsometricTensor, O::AbstractMatrix{<:MPOTensor}, BL::AdjointLeftIsometricTensor)
        return new{1, N}([AL,], O, [BL,])
    end
end
