"""
    mutable struct TransferMatrix <: AbstractTransferMatrix
        const B::Union{AdjointMPSTensor, AdjointLeftIsometricMPSTensor, AdjointRightIsometricMPSTensor}
        const A::Union{MPSTensor, LeftIsometricMPSTensor, RightIsometricMPSTensor}
    end

Wrapper type for transfer matrix.

Graphic presentation:

    —— B ——
       |
    —— A ——

# Constructors
    TransferMatrix(B::Union{AdjointMPSTensor, AdjointLeftIsometricMPSTensor, AdjointRightIsometricMPSTensor}, A::Union{MPSTensor, LeftIsometricMPSTensor, RightIsometricMPSTensor})
    TransferMatrix(A::Union{MPSTensor, LeftIsometricMPSTensor, RightIsometricMPSTensor}, B::Union{AdjointMPSTensor, AdjointLeftIsometricMPSTensor, AdjointRightIsometricMPSTensor})
"""
mutable struct TransferMatrix <: AbstractTransferMatrix
    const B::Union{AdjointMPSTensor, AdjointLeftIsometricMPSTensor, AdjointRightIsometricMPSTensor}
    const A::Union{MPSTensor, LeftIsometricMPSTensor, RightIsometricMPSTensor}

    function TransferMatrix(B::Union{AdjointMPSTensor, AdjointLeftIsometricMPSTensor, AdjointRightIsometricMPSTensor}, A::Union{MPSTensor, LeftIsometricMPSTensor, RightIsometricMPSTensor})
        return new(B, A)
    end
    function TransferMatrix(A::Union{MPSTensor, LeftIsometricMPSTensor, RightIsometricMPSTensor}, B::Union{AdjointMPSTensor, AdjointLeftIsometricMPSTensor, AdjointRightIsometricMPSTensor})
        return new(B, A)
    end
end
