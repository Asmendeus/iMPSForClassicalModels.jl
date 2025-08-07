"""
    mutable struct TransferMatrix{TA, TB} <: AbstractTransferMatrix where {TA<:Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}, TB<:Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}}
        A::TA
        B::TB
    end

Wrapper type for transfer matrix.

Graphic presentation:

    —— B ——
       |
    —— A ——

# Constructors
    TransferMatrix{TA, TB}(A::TA, B::TB) where {Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}, TB<:Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}}
    TransferMatrix(A::Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}, B::Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor})
"""
mutable struct TransferMatrix{TA, TB} <: AbstractTransferMatrix where {TA<:Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}, TB<:Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}}
    A::TA
    B::TB

    function TransferMatrix{TA, TB}(A::TA, B::TB) where {TA<:Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}, TB<:Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}}
        space(A, 2) == space(B, 3) || throw(SpaceMismatch("$(space(A, 2)) ≠ $(space(B, 3))"))
        return new{TA, TB}(A, B)
    end
    TransferMatrix(A::Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}, B::Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}) = TransferMatrix{typeof(A), typeof(B)}(A, B)
end

function Base.show(io::IO, T::TransferMatrix{TA, TB}) where {TA, TB}
    # data type
    println(io, "TransferMatrix{$TA, $TB}:")
    println(io)

    # graphic presentation
    println(io, " —— B ——")
    println(io, repeat(" ", 4) * "|")
    println(io, " —— A ——")
    println(io)

    # local tensors
    println(io, "A:")
    show(io, T.A)
    println(io)
    println(io)
    println(io, "B:")
    show(io, T.B)
    println(io)
end