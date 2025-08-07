"""
    mutable struct TransferMatrix <: AbstractTransferMatrix
        const B::Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}
        const A::Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}
    end

Wrapper type for transfer matrix.

Graphic presentation:

    —— B ——
       |
    —— A ——

# Constructors
    TransferMatrix(B::Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}, A::Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor})
    TransferMatrix(A::Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}, B::Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor})
"""
mutable struct TransferMatrix <: AbstractTransferMatrix
    const B::Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}
    const A::Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}

    function TransferMatrix(B::Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}, A::Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor})
        return new(B, A)
    end
    function TransferMatrix(A::Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}, B::Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor})
        return new(B, A)
    end
end
