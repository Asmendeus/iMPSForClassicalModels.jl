"""
    mutable struct LeftEnvironment{L} <: AbstractEnvironment{L, 3}
        const AL::AbstractVector{LeftIsometricMPSTensor}
        const O::AbstractVector{MPOTensor}
        const BL::AbstractVector{AdjointLeftIsometricMPSTensor}
    end

Wrapper type for left environment tensor's generating environment.

Graphic presentation:

    —— BL1 —— ... —— BLL ——
       |
    —— O1  —— ... —— OL  ——
       |
    —— AL1 —— ... —— ALL ——

# Constructors
    LeftEnvironment{L}(AL::AbstractVector{<:LeftIsometricMPSTensor}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricMPSTensor})
    LeftEnvironment(AL::AbstractVector{<:LeftIsometricMPSTensor}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricMPSTensor})
    LeftEnvironment(AL::LeftIsometricMPSTensor, O::MPOTensor, BL::AdjointLeftIsometricMPSTensor)
"""
mutable struct LeftEnvironment{L} <: AbstractEnvironment{L, 3}
    const AL::AbstractVector{LeftIsometricMPSTensor}
    const O::AbstractVector{MPOTensor}
    const BL::AbstractVector{AdjointLeftIsometricMPSTensor}

    LeftEnvironment{L}(AL::AbstractVector{<:LeftIsometricMPSTensor}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricMPSTensor}) where L = new{L}(AL, O, BL)
    function LeftEnvironment(AL::AbstractVector{<:LeftIsometricMPSTensor}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricMPSTensor})
        L = length(AL)
        (L == length(O) == length(BL)) || throw(ArgumentError("The lengths of `AL`, `O` and `BL` do not match"))
        return new{L}(AL, O, BL)
    end
    function LeftEnvironment(AL::LeftIsometricMPSTensor, O::MPOTensor, BL::AdjointLeftIsometricMPSTensor)
        return new{1}([AL,], [O,], [BL,])
    end
end
