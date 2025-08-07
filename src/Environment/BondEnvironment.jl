"""
    mutable struct BondEnvironment <: AbstractEnvironment{2}
        const CB::AbstractVector{<:AdjointMPSBondTensor}
        const CA::AbstractVector{<:MPSBondTensor}
    end

Wrapper type for mid environment tensor's generating environment.

Graphic presentation:

    —— CB ——

    —— CA ——

# Constructors
    BondEnvironment(CB::AbstractVector{<:AdjointMPSBondTensor}, CA::AbstractVector{<:MPSBondTensor})
    BondEnvironment(CA::AbstractVector{<:MPSBondTensor}, CB::AbstractVector{<:AdjointMPSBondTensor})
"""
mutable struct BondEnvironment <: AbstractEnvironment{2}
    const CB::AbstractVector{<:AdjointMPSBondTensor}
    const CA::AbstractVector{<:MPSBondTensor}

    BondEnvironment(CB::AbstractVector{<:AdjointMPSBondTensor}, CA::AbstractVector{<:MPSBondTensor}) = new(CB, CA)
    BondEnvironment(CA::AbstractVector{<:MPSBondTensor}, CB::AbstractVector{<:AdjointMPSBondTensor}) = new(CB, CA)
end
