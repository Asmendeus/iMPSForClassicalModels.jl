"""
    mutable struct TransferMatrix{L, R} <: AbstractTransferMatrix
        const A::AbstractVector{<:Union{LocalTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}
        const B::AbstractVector{<:Union{AdjointLocalTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}}
    end

Wrapper type for transfer matrix.

Graphic presentation:

    -- A1 -- ... -- AL --     -- A1 -- ... -- AL --
       |            |            ||           ||
    -- B1 -- ... -- BL --     -- B1 -- ... -- BL --

# Constructors
    TransferMatrix{L, R}(A::AbstractVector{<:Union{LocalTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}, B::AbstractVector{<:Union{AdjointLocalTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}})
    TransferMatrix(A::AbstractVector{<:Union{LocalTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}, B::AbstractVector{<:Union{AdjointLocalTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}})
    TransferMatrix(A::Union{LocalTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}, B::Union{AdjointLocalTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}})
"""
mutable struct TransferMatrix{L, R} <: AbstractTransferMatrix{L}
    A::AbstractVector{<:Union{LocalTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}
    B::AbstractVector{<:Union{AdjointLocalTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}}

    function TransferMatrix{L, R}(A::AbstractVector{<:Union{LocalTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}, B::AbstractVector{<:Union{AdjointLocalTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}}) where {L, R}
        (L == length(A) == length(B)) || throw(ArgumentError("Mismatched lengths: ($L, $(length(A)), $(length(B)))"))
        R == 2 && throw(ArgumentError("Not support bond tensors in forming transfer matrix"))

        pspace_A = [otimes([space(A[l], i) for i in vcat([2,], 3:R-1)]...) for l in 1:L]
        pspace_B = [otimes([space(B[l], i) for i in vcat([R,], 1:R-3)]...) for l in 1:L]
        pspace_A == pspace_B || throw(SpaceMismatch("Mismatched physical spaces: $(pspace_A) ≠ $(pspace_B))"))

        aspace_Alt = [space(A[l], 1) for l in 2:L]
        aspace_Art = [space(A[l], R) for l in 1:L-1]
        aspace_Alt == aspace_Art || throw(SpaceMismatch("Mismatched virtual spaces of MPS/MPO tensors: $(aspace_Alt) ≠ $(aspace_Art))"))

        aspace_Blt = [space(B[l], R-1) for l in 2:L]
        aspace_Brt = [space(B[l], R-2) for l in 1:L-1]
        aspace_Blt == aspace_Brt || throw(SpaceMismatch("Mismatched virtual spaces of adjoint MPS/MPO tensors: $(aspace_Blt) ≠ $(aspace_Brt))"))

        return new{L, R}(A, B)
    end
    function TransferMatrix(A::AbstractVector{<:Union{LocalTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}, B::AbstractVector{<:Union{AdjointLocalTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}}) where R
        L = length(A)
        return TransferMatrix{L, R}(A, B)
    end
    function TransferMatrix(A::Union{LocalTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}, B::Union{AdjointLocalTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}) where R
        return TransferMatrix{1, R}([A,], [B,])
    end
end

const MPSTransferMatrix{L} = TransferMatrix{L, 3}
const MPOTransferMatrix{L} = TransferMatrix{L, 4}
const MPSOrMPOTransferMatrix{L} = Union{MPSTransferMatrix{L}, MPOTransferMatrix{L}}
