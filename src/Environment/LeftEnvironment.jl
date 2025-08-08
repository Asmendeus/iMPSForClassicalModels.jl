"""
    mutable struct LeftEnvironment{L, R} <: AbstractEnvironment{3}
        const AL::AbstractVector{LeftIsometricTensor{R}}
        const O::AbstractVector{MPOTensor}
        const BL::AbstractVector{AdjointLeftIsometricTensor{R}}
    end

Wrapper type for left environment tensor's generating environment.

Graphic presentation:

    -- AL1 -- ... -- ALL --     -- AL1 -- ... -- ALL --
       |             |             | |           | |
    —— O1  —— ... —— OL  ——     —— O1| —— ... —— OL| ——
       |             |             | |           | |
    -- BL1 -- ... -- BLL --     -- BL1 -- ... -- BLL --

# Constructors
    LeftEnvironment{L, R}(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}})
    LeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}})
    LeftEnvironment(AL::LeftIsometricTensor{R}, O::MPOTensor, BL::AdjointLeftIsometricTensor{R})
"""
mutable struct LeftEnvironment{L, R} <: AbstractEnvironment{3}
    const AL::AbstractVector{LeftIsometricTensor{R}}
    const O::AbstractVector{MPOTensor}
    const BL::AbstractVector{AdjointLeftIsometricTensor{R}}

    function LeftEnvironment{L, R}(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}}) where {L, R}
        (L == length(AL) == length(O) == length(BL)) || throw(ArgumentError("Mismatched lengths: ($L, $(length(AL)), $(length(O)), $(length(BL)))"))
        R == 2 && throw(ArgumentError("Not support bond tensors in forming left environment"))

        pspace_ALdn = [space(AL[l], 2) for l in 1:L]
        pspace_Oup = [space(O[l], 3) for l in 1:L]
        pspace_ALdn == pspace_Oup || throw(SpaceMismatch("Mismatched physical spaces of MPS/MPO tensors and MPO Tensors: $(pspace_ALdn) ≠ $(pspace_Oup))"))

        pspace_BLup = [space(BL[l], R) for l in 1:L]
        pspace_Odn = [space(O[l], 2) for l in 1:L]
        pspace_BLup == pspace_Odn || throw(SpaceMismatch("Mismatched physical spaces of MPO tensors and adjoint MPS/MPO tensors: $(pspace_BLup) ≠ $(pspace_Odn))"))

        if R > 3
            pspace_ALup = [otimes([space(AL[l], i) for i in 3:R-1]...) for l in 1:L]
            pspace_BLdn = [otimes([space(BL[l], i) for i in 1:R-3]...) for l in 1:L]
            pspace_ALup == pspace_BLdn || throw(SpaceMismatch("Mismatched physical spaces of MPS/MPO tensors and adjoint MPS/MPO tensors: $(pspace_BLup) ≠ $(pspace_Odn))"))
        end

        aspace_ALlt = [space(AL[l], 1) for l in 2:L]
        aspace_ALrt = [space(AL[l], R) for l in 1:L-1]
        aspace_ALlt == aspace_ALrt || throw(SpaceMismatch("Mismatched intra-MPS/MPO virtual spaces of MPS/MPO tensors: $(aspace_ALlt) ≠ $(aspace_ALrt))"))

        aspace_Olt = [space(O[l], 1) for l in 2:L]
        aspace_Ort = [space(O[l], 4) for l in 1:L-1]
        aspace_Olt == aspace_Ort || throw(SpaceMismatch("Mismatched intra-MPO virtual spaces of MPO tensors: $(aspace_Olt) ≠ $(aspace_Ort))"))

        aspace_BLlt = [space(BL[l], R-1) for l in 2:L]
        aspace_BLrt = [space(BL[l], R-2) for l in 1:L-1]
        aspace_BLlt == aspace_BLrt || throw(SpaceMismatch("Mismatched intra-MPS/MPO virtual spaces of adjoint MPS/MPO tensors: $(aspace_BLlt) ≠ $(aspace_BLrt))"))

        return new{L, R}(AL, O, BL)
    end
    function LeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}}) where R
        L = length(AL)
        return LeftEnvironment{L, R}(AL, O, BL)
    end
    function LeftEnvironment(AL::LeftIsometricTensor{R}, O::MPOTensor, BL::AdjointLeftIsometricTensor{R}) where R
        return LeftEnvironment{1, R}([AL,], [O,], [BL,])
    end
end

const LeftMPSEnvironment{L} = LeftEnvironment{L, 3}
const LeftMPOEnvironment{L} = LeftEnvironment{L, 4}
