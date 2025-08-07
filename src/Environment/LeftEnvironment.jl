"""
    mutable struct LeftEnvironment <: AbstractEnvironment{3}
        const BL::AbstractVector{<:AdjointLeftIsometricTensor}
        const O::AbstractVector{<:MPOTensor}
        const AL::AbstractVector{<:LeftIsometricTensor}
    end

Wrapper type for left environment tensor's generating environment.

Graphic presentation:

    —— BL ——
       |
    —— O  ——
       |
    —— AL ——

# Constructors
    LeftEnvironment(BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractVector{<:MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor})
    LeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
    LeftEnvironment(BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractVector{<:MPOTensor})
    LeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractVector{<:MPOTensor})
    LeftEnvironment(O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor})
    LeftEnvironment(O::AbstractVector{<:MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
"""
mutable struct LeftEnvironment <: AbstractEnvironment{3}
    const BL::AbstractVector{<:AdjointLeftIsometricTensor}
    const O::AbstractVector{<:MPOTensor}
    const AL::AbstractVector{<:LeftIsometricTensor}

    LeftEnvironment(BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractVector{<:MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor}) = new(BL, O, AL)
    LeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}) = new(BL, O, AL)
    LeftEnvironment(BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractVector{<:MPOTensor}) = new(BL, O, AL)
    LeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractVector{<:MPOTensor}) = new(BL, O, AL)
    LeftEnvironment(O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor}) = new(BL, O, AL)
    LeftEnvironment(O::AbstractVector{<:MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}) = new(BL, O, AL)
end
