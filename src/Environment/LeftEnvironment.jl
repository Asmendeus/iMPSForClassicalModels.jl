"""
    mutable struct LeftEnvironment{L} <: AbstractEnvironment{L, 3}
        const AL::AbstractVector{LeftIsometricTensor}
        const O::AbstractVector{MPOTensor}
        const BL::AbstractVector{AdjointLeftIsometricTensor}
    end

Wrapper type for left environment tensor's generating environment.

Graphic presentation:

    —— BL1 —— ... —— BLL ——
       |             |
    —— O1  —— ... —— OL  ——
       |             |
    —— AL1 —— ... —— ALL ——

# Constructors
    LeftEnvironment{L}(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
    LeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
    LeftEnvironment(AL::LeftIsometricTensor, O::MPOTensor, BL::AdjointLeftIsometricTensor)
"""
mutable struct LeftEnvironment{L} <: AbstractEnvironment{L, 3}
    const AL::AbstractVector{LeftIsometricTensor}
    const O::AbstractVector{MPOTensor}
    const BL::AbstractVector{AdjointLeftIsometricTensor}

    LeftEnvironment{L}(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}) where L = new{L}(AL, O, BL)
    function LeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
        L = length(AL)
        (L == length(O) == length(BL)) || throw(ArgumentError("The lengths of `AL`, `O` and `BL` do not match"))
        return new{L}(AL, O, BL)
    end
    function LeftEnvironment(AL::LeftIsometricTensor, O::MPOTensor, BL::AdjointLeftIsometricTensor)
        return new{1}([AL,], [O,], [BL,])
    end
end
