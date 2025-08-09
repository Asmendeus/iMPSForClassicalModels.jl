"""
    mutable struct MultirowRightEnvironment{N, L, R} <: AbstractEnvironment{N}
        const AR::AbstractVector{RightIsometricTensor{R}}
        const O::AbstractMatrix{MPOTensor}
        const BR::AbstractVector{AdjointRightIsometricTensor{R}}
    end

Wrapper type for general right environment tensor's generating environment.

Graphic presentation:
                                  (a)           (b)
                                   |             |
    -- AR1 -- ... -- ARL --     -- AR1 -- ... -- ARL --
       |             |             |             |
    —— O11 —— ... —— O1L ——     —— O11 —— ... —— O1L ——
       |             |             |             |
       ⋮              ⋮             ⋮              ⋮
       |             |             |             |
    —— OW1 —— ... —— OWL ——     —— OW1 —— ... —— OWL ——         W = N - 2
       |             |             |             |
    -- BR1 -- ... -- BRL --     -- BR1 -- ... -- BRL --
                                   |             |
                                  (a)           (b)

# Constructors
    MultirowRightEnvironment{N, L, R}(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}})
    MultirowRightEnvironment(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}})
    MultirowRightEnvironment(AR::RightIsometricTensor{R}, O::AbstractVecOrMat{<:MPOTensor}, BR::AdjointRightIsometricTensor{R})
"""
mutable struct MultirowRightEnvironment{N, L, R} <: AbstractEnvironment{N}
    const AR::AbstractVector{RightIsometricTensor{R}}
    const O::AbstractMatrix{MPOTensor}
    const BR::AbstractVector{AdjointRightIsometricTensor{R}}

    function MultirowRightEnvironment{N, L, R}(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}}) where {N, L, R}
        (L == length(AR) == size(O, 2) == length(BR)) || throw(ArgumentError("Mismatched lengths: ($L, $(length(AR)), $(length(O)), $(length(BR)))"))
        (N == size(O, 1) + 2) || throw(ArgumentError("Mismatched widths: ($N, $(size(O, 1) + 2))"))
        R == 2 && throw(ArgumentError("Not support bond tensors in forming multirow right environment"))

        pspace_ALdn = [space(AR[l], 2) for l in 1:L]
        pspace_Oup = [space(O[1, l], 3) for l in 1:L]
        pspace_ALdn == pspace_Oup || throw(SpaceMismatch("Mismatched physical spaces of MPS/MPO tensors and MPO Tensors: $(pspace_ALdn) ≠ $(pspace_Oup))"))

        for w in 1:N-3
            pspace_Odn_w = [space(O[w, l], 2) for l in 1:L]
            pspace_Oup_w = [space(O[w+1, l], 3) for l in 1:L]
            pspace_Odn_w == pspace_Oup_w || throw(SpaceMismatch("Mismatched intra-MPO physical spaces of MPO tensors: $(pspace_Odn_w) ≠ $(pspace_Oup_w))"))
        end

        pspace_BLup = [space(BR[l], R) for l in 1:L]
        pspace_Odn = [space(O[end, l], 2) for l in 1:L]
        pspace_BLup == pspace_Odn || throw(SpaceMismatch("Mismatched physical spaces of MPO tensors and adjoint MPS/MPO tensors: $(pspace_BLup) ≠ $(pspace_Odn))"))

        if R > 3
            pspace_ALup = [otimes([space(AR[l], i) for i in 3:R-1]...) for l in 1:L]
            pspace_BLdn = [otimes([space(BR[l], i) for i in 1:R-3]...) for l in 1:L]
            pspace_ALup == pspace_BLdn || throw(SpaceMismatch("Mismatched physical spaces of MPS/MPO tensors and adjoint MPS/MPO tensors: $(pspace_BLup) ≠ $(pspace_Odn))"))
        end

        aspace_ALlt = [space(AR[l], 1) for l in 2:L]
        aspace_ALrt = [space(AR[l], R) for l in 1:L-1]
        aspace_ALlt == aspace_ALrt || throw(SpaceMismatch("Mismatched intra-MPS/MPO virtual spaces of MPS/MPO tensors: $(aspace_ALlt) ≠ $(aspace_ALrt))"))

        for w in 1:N-2
            aspace_Olt_w = [space(O[w, l], 1) for l in 2:L]
            aspace_Ort_w = [space(O[w, l], 4) for l in 1:L-1]
            aspace_Olt_w == aspace_Ort_w || throw(SpaceMismatch("Mismatched intra-MPO virtual spaces of MPO tensors: $(aspace_Olt_w) ≠ $(aspace_Ort_w))"))
        end

        aspace_BLlt = [space(BR[l], R-1) for l in 2:L]
        aspace_BLrt = [space(BR[l], R-2) for l in 1:L-1]
        aspace_BLlt == aspace_BLrt || throw(SpaceMismatch("Mismatched intra-MPS/MPO virtual spaces of adjoint MPS/MPO tensors: $(aspace_BLlt) ≠ $(aspace_BLrt))"))

        return new{N, L, R}(AR, O, BR)
    end
    function MultirowRightEnvironment(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}}) where R
        L = length(AR)
        N = size(O, 1) + 2
        return MultirowRightEnvironment{N, L, R}(AR, O, BR)
    end
    function MultirowRightEnvironment(AR::RightIsometricTensor{R}, O::AbstractVecOrMat{<:MPOTensor}, BR::AdjointRightIsometricTensor{R}) where R
        N = size(O, 1) + 2
        return MultirowRightEnvironment{N, 1, R}([AR,], reshape(O, (N-2, 1)), [BR,])
    end
end

const MultirowRightMPSEnvironment{N, L} = MultirowRightEnvironment{N, L, 3}
const MultirowRightMPOEnvironment{N, L} = MultirowRightEnvironment{N, L, 4}
