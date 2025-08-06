"""
    —— BC ——
       |
    —— O  ——
       |
    —— AC ——
"""
mutable struct MidEnvironment <: AbstractEnvironment
    const BC::AbstractVector{AdjointMPSTensor}
    const O::AbstractVector{MPOTensor}
    const AC::AbstractVector{MPSTensor}
end
