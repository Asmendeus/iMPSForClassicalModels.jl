"""
    —— BL ——
       |
    —— O  ——
       |
    —— AL ——
"""
mutable struct LeftEnvironment <: AbstractEnvironment
    const BL::AbstractVector{AdjointLeftIsometricTensor}
    const O::AbstractVector{MPOTensor}
    const AL::AbstractVector{LeftIsometricTensor}
end
