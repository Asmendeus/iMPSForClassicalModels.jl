"""
    —— BR ——
       |
    —— O  ——
       |
    —— AR ——
"""
mutable struct RightEnvironment <: AbstractEnvironment
    const BR::AbstractVector{AdjointRightIsometricTensor}
    const O::AbstractVector{MPOTensor}
    const AR::AbstractVector{RightIsometricTensor}
end
