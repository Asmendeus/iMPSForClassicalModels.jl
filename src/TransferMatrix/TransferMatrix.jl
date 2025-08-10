"""
    mutable struct TransferMatrix{L, R} <: AbstractTransferMatrix{L, R}
        const A::AbstractVector{LocalTensor{R}}
        const B::AbstractVector{AdjointLocalTensor{R}}
    end

Wrapper type for transfer matrix.

Graphic presentation:
                                    (a1)           (aL)      bonds with same tag are connected
                                     |              |
    -- A[1] -- ... -- A[L] --     -- A[1] -- ... -- A[L] --
       |              |              |              |
    -- B[1] -- ... -- B[L] --     -- B[1] -- ... -- B[L] --
                                     |              |
                                    (a1)           (aL)

# Constructors
    TransferMatrix{R}(A::AbstractVector{<:LocalTensor{R}}, B::AbstractVector{<:AdjointLocalTensor{R}})
    TransferMatrix(A::AbstractVector{<:LocalTensor{R}}, B::AbstractVector{<:AdjointLocalTensor{R}})
    TransferMatrix(A::LocalTensor{R}, B::AdjointLocalTensor{R})
"""
mutable struct TransferMatrix{L, R} <: AbstractTransferMatrix{L, R}
    const A::AbstractVector{LocalTensor{R}}
    const B::AbstractVector{AdjointLocalTensor{R}}

    function TransferMatrix{L, R}(A::AbstractVector{<:LocalTensor{R}}, B::AbstractVector{<:AdjointLocalTensor{R}}) where {L, R}
        (L == length(A) == length(B)) || throw(ArgumentError("Mismatched lengths: ($L, $(length(A)), $(length(B)))"))
        R == 2 && throw(ArgumentError("Not support bond tensors in forming transfer matrix"))

        for l in 1:L
            pspace_A_l = [space(A[l], i) for i in vcat([2,], 3:R-1)]
            pspace_B_l = [space(B[l], i) for i in vcat([R,], 1:R-3)]
            pspace_A_l == pspace_B_l || throw(SpaceMismatch("Mismatched $(l)-th physical space: $(pspace_A_l) ≠ $(pspace_B_l))"))
        end

        for l in 1:L
            aspace_Art_l = space(A[l], R)
            aspace_Alt_l = space(A[mod(l, L)+1], 1)
            aspace_Art_l == aspace_Alt_l || throw(SpaceMismatch("Mismatched virtual space between $(l)-th and $(mod(l, L)+1)-th local tensors: $(aspace_Art_l) ≠ $(aspace_Alt_l))"))

            aspace_Brt_l = space(B[l], R-2)
            aspace_Blt_l = space(B[mod(l, L)+1], R-1)
            aspace_Brt_l == aspace_Blt_l || throw(SpaceMismatch("Mismatched virtual space between $(l)-th and $(mod(l, L)+1)-th adjoint local tensors: $(aspace_Brt_l) ≠ $(aspace_Blt_l))"))
        end

        return new{L, R}(A, B)
    end
    TransferMatrix(A::AbstractVector{<:LocalTensor{R}}, B::AbstractVector{<:AdjointLocalTensor{R}}) where R = TransferMatrix{length(A), R}(A, B)
    TransferMatrix(A::LocalTensor{R}, B::AdjointLocalTensor{R}) where R = TransferMatrix{1, R}([A,], [B,])
end

const MPSTransferMatrix{L} = TransferMatrix{L, 3}
const MPOTransferMatrix{L} = TransferMatrix{L, 4}
