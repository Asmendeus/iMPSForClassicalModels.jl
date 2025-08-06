"""
    mutable struct TransferMatrix <: AbstractTransferMatrix
        const B::AbstractVector{Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}}
        const A::AbstractVector{Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}}
    end

Wrapper type for transfer matrix.

Graphic presentation:

    —— B ——
       |
    —— A ——
"""
mutable struct TransferMatrix <: AbstractTransferMatrix
    const B::AbstractVector{Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}}
    const A::AbstractVector{Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}}
end
