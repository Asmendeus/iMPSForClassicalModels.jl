"""
    mutable struct RightEnvironment{L, R} <: AbstractEnvironment{3}
        const AR::AbstractVector{RightIsometricTensor{R}}
        const O::AbstractVector{MPOTensor}
        const BR::AbstractVector{AdjointRightIsometricTensor{R}}
    end

Wrapper type for right environment tensor's generating environment.

Graphic presentation:

    -- AR1 -- ... -- ARL --     -- AR1 -- ... -- ARL --
       |             |             | |           | |
    —— O1  —— ... —— OL  ——     —— O1| —— ... —— OL| ——
       |             |             | |           | |
    -- BR1 -- ... -- BRL --     -- BR1 -- ... -- BRL --

# Constructors
    RightEnvironment{L, R}(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}})
    RightEnvironment(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}})
    RightEnvironment(AR::RightIsometricTensor{R}, O::MPOTensor, BR::AdjointRightIsometricTensor{R})
"""
mutable struct RightEnvironment{L, R} <: AbstractEnvironment{3}
    const AR::AbstractVector{RightIsometricTensor{R}}
    const O::AbstractVector{MPOTensor}
    const BR::AbstractVector{AdjointRightIsometricTensor{R}}

    function RightEnvironment{L, R}(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}}) where {L, R}
        (L == length(AR) == length(O) == length(BR)) || throw(ArgumentError("Mismatched lengths: ($L, $(length(AR)), $(length(O)), $(length(BR)))"))
        R == 2 && throw(ArgumentError("Not support bond tensors in forming right environment"))

        pspace_ARdn = [space(AR[l], 2) for l in 1:L]
        pspace_Oup = [space(O[l], 3) for l in 1:L]
        pspace_ARdn == pspace_Oup || throw(SpaceMismatch("Mismatched physical spaces of MPS/MPO tensors and MPO Tensors: $(pspace_ARdn) ≠ $(pspace_Oup))"))

        pspace_BRup = [space(BR[l], R) for l in 1:L]
        pspace_Odn = [space(O[l], 2) for l in 1:L]
        pspace_BRup == pspace_Odn || throw(SpaceMismatch("Mismatched physical spaces of MPO tensors and adjoint MPS/MPO tensors: $(pspace_BRup) ≠ $(pspace_Odn))"))

        if R > 3
            pspace_ARup = [otimes([space(AR[l], i) for i in 3:R-1]...) for l in 1:L]
            pspace_BRdn = [otimes([space(BR[l], i) for i in 1:R-3]...) for l in 1:L]
            pspace_ARup == pspace_BRdn || throw(SpaceMismatch("Mismatched physical spaces of MPS/MPO tensors and adjoint MPS/MPO tensors: $(pspace_BRup) ≠ $(pspace_Odn))"))
        end

        aspace_ARlt = [space(AR[l], 1) for l in 2:L]
        aspace_ARrt = [space(AR[l], R) for l in 1:L-1]
        aspace_ARlt == aspace_ARrt || throw(SpaceMismatch("Mismatched intra-MPS/MPO virtual spaces of MPS/MPO tensors: $(aspace_ARlt) ≠ $(aspace_ARrt))"))

        aspace_Olt = [space(O[l], 1) for l in 2:L]
        aspace_Ort = [space(O[l], 4) for l in 1:L-1]
        aspace_Olt == aspace_Ort || throw(SpaceMismatch("Mismatched intra-MPO virtual spaces of MPO tensors: $(aspace_Olt) ≠ $(aspace_Ort))"))

        aspace_BRlt = [space(BR[l], R-1) for l in 2:L]
        aspace_BRrt = [space(BR[l], R-2) for l in 1:L-1]
        aspace_BRlt == aspace_BRrt || throw(SpaceMismatch("Mismatched intra-MPS/MPO virtual spaces of adjoint MPS/MPO tensors: $(aspace_BRlt) ≠ $(aspace_BRrt))"))

        return new{L, R}(AR, O, BR)
    end
    function RightEnvironment(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}}) where R
        L = length(AR)
        return RightEnvironment{L, R}(AR, O, BR)
    end
    function RightEnvironment(AR::RightIsometricTensor{R}, O::MPOTensor, BR::AdjointRightIsometricTensor{R}) where R
        return RightEnvironment{1, R}([AR,], [O,], [BR,])
    end
end

const RightMPSEnvironment{L} = RightEnvironment{L, 3}
const RightMPOEnvironment{L} = RightEnvironment{L, 4}
