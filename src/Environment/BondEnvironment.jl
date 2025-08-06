"""
    —— CB ——

    —— CA ——
"""
mutable struct BondEnvironment <: AbstractEnvironment
    const CB::AbstractVector{AdjointMPSBondTensor}
    const CA::AbstractVector{MPSBondTensor}
end
