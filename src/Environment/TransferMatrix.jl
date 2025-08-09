"""
     abstract type AbstractTransferMatrix{R} <: AbstractEnvironment{2}

Wrapper type for transfer matrix.
"""
abstract type AbstractTransferMatrix{R} <: AbstractEnvironment{2} end

const AbstractMPSTransferMatrix = AbstractTransferMatrix{3}
const AbstractMPOTransferMatrix = AbstractTransferMatrix{4}

"""
    mutable struct TransferMatrix{R} <: AbstractEnvironment{2}
        const A::LocalTensor{R}
        const B::AdjointLocalTensor{R}
    end

Wrapper type for transfer matrix.

Graphic presentation:

    -- A --     -- A --
       |           ||
    -- B --     -- B --

# Constructors
    TransferMatrix{R}(A::LocalTensor{R}, B::AdjointLocalTensor{R})
    TransferMatrix(A::LocalTensor{R}, B::AdjointLocalTensor{R})
"""
mutable struct TransferMatrix{R} <: AbstractTransferMatrix{R}
    A::LocalTensor{R}
    B::AdjointLocalTensor{R}

    function TransferMatrix{R}(A::LocalTensor{R}, B::AdjointLocalTensor{R}) where R
        R == 2 && throw(ArgumentError("Not support bond tensors in forming transfer matrix"))

        pspace_A = [space(A, i) for i in vcat([2,], 3:R-1)]
        pspace_B = [space(B, i) for i in vcat([R,], 1:R-3)]
        pspace_A == pspace_B || throw(SpaceMismatch("Mismatched physical spaces: $(pspace_A) ≠ $(pspace_B))"))

        return new{R}(A, B)
    end
    TransferMatrix(A::LocalTensor{R}, B::AdjointLocalTensor{R}) where R = TransferMatrix{R}(A, B)
end

const MPSTransferMatrix = TransferMatrix{3}
const MPOTransferMatrix = TransferMatrix{4}
