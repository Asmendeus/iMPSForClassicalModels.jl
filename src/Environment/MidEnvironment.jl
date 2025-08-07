"""
    mutable struct MidEnvironment <: AbstractEnvironment{3}
        const BC::AbstractVector{<:AdjointMPSTensor}
        const O::AbstractVector{<:MPOTensor}
        const AC::AbstractVector{<:MPSTensor}
    end

Wrapper type for mid environment tensor's generating environment.

Graphic presentation:

    —— BC ——
       |
    —— O  ——
       |
    —— AC ——

# Constructors
    MidEnvironment(BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractVector{<:MPOTensor}, AC::AbstractVector{<:MPSTensor})
    MidEnvironment(AC::AbstractVector{<:MPSTensor}, O::AbstractVector{<:MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor})
    MidEnvironment(BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor}, O::AbstractVector{<:MPOTensor})
    MidEnvironment(AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractVector{<:MPOTensor})
    MidEnvironment(O::AbstractVector{<:MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor})
    MidEnvironment(O::AbstractVector{<:MPOTensor}, AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor})
"""
mutable struct MidEnvironment <: AbstractEnvironment{3}
    const BC::AbstractVector{<:AdjointMPSTensor}
    const O::AbstractVector{<:MPOTensor}
    const AC::AbstractVector{<:MPSTensor}

    MidEnvironment(BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractVector{<:MPOTensor}, AC::AbstractVector{<:MPSTensor}) = new(BC, O, AC)
    MidEnvironment(AC::AbstractVector{<:MPSTensor}, O::AbstractVector{<:MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor}) = new(BC, O, AC)
    MidEnvironment(BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor}, O::AbstractVector{<:MPOTensor}) = new(BC, O, AC)
    MidEnvironment(AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractVector{<:MPOTensor}) = new(BC, O, AC)
    MidEnvironment(O::AbstractVector{<:MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor}) = new(BC, O, AC)
    MidEnvironment(O::AbstractVector{<:MPOTensor}, AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor}) = new(BC, O, AC)
end
