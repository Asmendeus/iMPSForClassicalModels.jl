"""
    mutable struct LeftEnvironment{L, R₁, R₂} <: AbstractEnvironment{3}
        const AL::AbstractVector{LeftIsometricTensor{R₁}}
        const O::AbstractVector{MPSTensor{R₂}}
        const BL::AbstractVector{AdjointLeftIsometricTensor{R₁}}
    end

Wrapper type for left environment tensor's generating environment.

Graphic presentation:

    —— AL1 —— ... —— ALL ——     —— AL1 —— ... —— ALL ——
       |             |             | |           | |
    —— O1  —— ... —— OL  ——     —— O1| —— ... —— OL| ——
       |             |             | |           | |
    —— BL1 —— ... —— BLL ——     —— BL1 —— ... —— BLL ——

# Constructors
    LeftEnvironment{R₁, R₂, L}(AL::AbstractVector{<:LeftIsometricTensor{R₁}}, O::AbstractVector{<:MPSTensor{R₂}}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R₁}})
    LeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractVector{<:MPSTensor{R₂}}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
    LeftEnvironment(AL::LeftIsometricTensor, O::MPOTensor, BL::AdjointLeftIsometricTensor)
"""
mutable struct LeftEnvironment{L, R₁, R₂} <: AbstractEnvironment{3}
    const AL::AbstractVector{LeftIsometricTensor{R₁}}
    const O::AbstractVector{MPSTensor{R₂}}
    const BL::AbstractVector{AdjointLeftIsometricTensor{R₁}}

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
