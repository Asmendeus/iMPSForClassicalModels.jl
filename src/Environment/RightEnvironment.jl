"""
    mutable struct RightEnvironment <: AbstractEnvironment{3}
        const BR::AbstractVector{<:AdjointRightIsometricTensor}
        const O::AbstractVector{<:MPOTensor}
        const AR::AbstractVector{<:RightIsometricTensor}
    end

Wrapper type for Right environment tensor's generating environment.

Graphic presentation:

    —— BR ——
       |
    —— O  ——
       |
    —— AR ——

# Constructors
    RightEnvironment(BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractVector{<:MPOTensor}, AR::AbstractVector{<:RightIsometricTensor})
    RightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
    RightEnvironment(BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor}, O::AbstractVector{<:MPOTensor})
    RightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractVector{<:MPOTensor})
    RightEnvironment(O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor})
    RightEnvironment(O::AbstractVector{<:MPOTensor}, AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
"""
mutable struct RightEnvironment <: AbstractEnvironment{3}
    const BR::AbstractVector{<:AdjointRightIsometricTensor}
    const O::AbstractVector{<:MPOTensor}
    const AR::AbstractVector{<:RightIsometricTensor}

    RightEnvironment(BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractVector{<:MPOTensor}, AR::AbstractVector{<:RightIsometricTensor}) = new(BR, O, AR)
    RightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}) = new(BR, O, AR)
    RightEnvironment(BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor}, O::AbstractVector{<:MPOTensor}) = new(BR, O, AR)
    RightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractVector{<:MPOTensor}) = new(BR, O, AR)
    RightEnvironment(O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor}) = new(BR, O, AR)
    RightEnvironment(O::AbstractVector{<:MPOTensor}, AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}) = new(BR, O, AR)
end
