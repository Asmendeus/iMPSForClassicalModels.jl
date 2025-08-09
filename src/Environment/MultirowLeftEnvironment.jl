"""
    mutable struct MultirowLeftEnvironment{N, L, R} <: AbstractEnvironment{N}
        const AL::AbstractVector{LeftIsometricTensor{R}}
        const O::AbstractMatrix{MPOTensor}
        const BL::AbstractVector{AdjointLeftIsometricTensor{R}}
    end

Wrapper type for general left environment tensor's generating environment.

Graphic presentation:
                                  (a)           (b)
                                   |             |
    -- AL1 -- ... -- ALL --     -- AL1 -- ... -- ALL --
       |             |             |             |
    —— O11 —— ... —— O1L ——     —— O11 —— ... —— O1L ——
       |             |             |             |
       ⋮              ⋮             ⋮              ⋮
       |             |             |             |
    —— OW1 —— ... —— OWL ——     —— OW1 —— ... —— OWL ——         W = N - 2
       |             |             |             |
    -- BL1 -- ... -- BLL --     -- BL1 -- ... -- BLL --
                                   |             |
                                  (a)           (b)

# Constructors
    MultirowLeftEnvironment{N, L, R}(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}})
    MultirowLeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}})
    MultirowLeftEnvironment(AL::LeftIsometricTensor{R}, O::AbstractVecOrMat{<:MPOTensor}, BL::AdjointLeftIsometricTensor{R})
"""
mutable struct MultirowLeftEnvironment{N, L, R} <: AbstractEnvironment{N}
    const AL::AbstractVector{LeftIsometricTensor{R}}
    const O::AbstractMatrix{MPOTensor}
    const BL::AbstractVector{AdjointLeftIsometricTensor{R}}

    function MultirowLeftEnvironment{N, L, R}(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}}) where {N, L, R}
        (L == length(AL) == size(O, 2) == length(BL)) || throw(ArgumentError("Mismatched lengths: ($L, $(length(AL)), $(length(O)), $(length(BL)))"))
        (N == size(O, 1) + 2) || throw(ArgumentError("Mismatched widths: ($N, $(size(O, 1) + 2))"))
        R == 2 && throw(ArgumentError("Not support bond tensors in forming multirow left environment"))

        pspace_ALdn = [space(AL[l], 2) for l in 1:L]
        pspace_Oup = [space(O[1, l], 3) for l in 1:L]
        pspace_ALdn == pspace_Oup || throw(SpaceMismatch("Mismatched physical spaces of MPS/MPO tensors and MPO Tensors: $(pspace_ALdn) ≠ $(pspace_Oup))"))

        for w in 1:N-3
            pspace_Odn_w = [space(O[w, l], 2) for l in 1:L]
            pspace_Oup_w = [space(O[w+1, l], 3) for l in 1:L]
            pspace_Odn_w == pspace_Oup_w || throw(SpaceMismatch("Mismatched intra-MPO physical spaces of MPO tensors: $(pspace_Odn_w) ≠ $(pspace_Oup_w))"))
        end

        pspace_BLup = [space(BL[l], R) for l in 1:L]
        pspace_Odn = [space(O[end, l], 2) for l in 1:L]
        pspace_BLup == pspace_Odn || throw(SpaceMismatch("Mismatched physical spaces of MPO tensors and adjoint MPS/MPO tensors: $(pspace_BLup) ≠ $(pspace_Odn))"))

        if R > 3
            pspace_ALup = [otimes([space(AL[l], i) for i in 3:R-1]...) for l in 1:L]
            pspace_BLdn = [otimes([space(BL[l], i) for i in 1:R-3]...) for l in 1:L]
            pspace_ALup == pspace_BLdn || throw(SpaceMismatch("Mismatched physical spaces of MPS/MPO tensors and adjoint MPS/MPO tensors: $(pspace_BLup) ≠ $(pspace_Odn))"))
        end

        aspace_ALlt = [space(AL[l], 1) for l in 2:L]
        aspace_ALrt = [space(AL[l], R) for l in 1:L-1]
        aspace_ALlt == aspace_ALrt || throw(SpaceMismatch("Mismatched intra-MPS/MPO virtual spaces of MPS/MPO tensors: $(aspace_ALlt) ≠ $(aspace_ALrt))"))

        for w in 1:N-2
            aspace_Olt_w = [space(O[w, l], 1) for l in 2:L]
            aspace_Ort_w = [space(O[w, l], 4) for l in 1:L-1]
            aspace_Olt_w == aspace_Ort_w || throw(SpaceMismatch("Mismatched intra-MPO virtual spaces of MPO tensors: $(aspace_Olt_w) ≠ $(aspace_Ort_w))"))
        end

        aspace_BLlt = [space(BL[l], R-1) for l in 2:L]
        aspace_BLrt = [space(BL[l], R-2) for l in 1:L-1]
        aspace_BLlt == aspace_BLrt || throw(SpaceMismatch("Mismatched intra-MPS/MPO virtual spaces of adjoint MPS/MPO tensors: $(aspace_BLlt) ≠ $(aspace_BLrt))"))

        return new{N, L, R}(AL, O, BL)
    end
    function MultirowLeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}}) where R
        L = length(AL)
        N = size(O, 1) + 2
        return MultirowLeftEnvironment{N, L, R}(AL, O, BL)
    end
    function MultirowLeftEnvironment(AL::LeftIsometricTensor{R}, O::AbstractVecOrMat{<:MPOTensor}, BL::AdjointLeftIsometricTensor{R}) where R
        N = size(O, 1) + 2
        return MultirowLeftEnvironment{N, 1, R}([AL,], reshape(O, (N-2, 1)), [BL,])
    end
end

const MultirowLeftMPSEnvironment{N, L} = MultirowLeftEnvironment{N, L, 3}
const MultirowLeftMPOEnvironment{N, L} = MultirowLeftEnvironment{N, L, 4}
