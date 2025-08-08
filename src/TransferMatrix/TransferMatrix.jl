"""
    mutable struct TransferMatrix{L, R} <: AbstractTransferMatrix
        const A::AbstractVector{<:Union{MPSTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}
        const B::AbstractVector{<:Union{AdjointMPSTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}}
    end

Wrapper type for transfer matrix.

Graphic presentation:

    -- A1 -- ... -- AL --     -- A1 -- ... -- AL --
       |            |            ||           ||
    -- B1 -- ... -- BL --     -- B1 -- ... -- BL --

# Constructors
    TransferMatrix{L, R}(A::AbstractVector{<:Union{MPSTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}, B::AbstractVector{<:Union{AdjointMPSTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}})
    TransferMatrix(A::AbstractVector{<:Union{MPSTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}, B::AbstractVector{<:Union{AdjointMPSTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}})
    TransferMatrix(A::Union{MPSTensor, LeftIsometricTensor, RightIsometricTensor}, B::Union{AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor})
"""
mutable struct TransferMatrix{L, R} <: AbstractTransferMatrix{L}
    A::AbstractVector{<:Union{MPSTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}
    B::AbstractVector{<:Union{AdjointMPSTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}}

    # function TransferMatrix(A::Union{MPSTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}, B::Union{AdjointMPSTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}) where R
    #     R == 2 && throw(ArgumentError("Not support bond tensors in forming transfer matrix"))
    #     pspace_A = [space(A, i) for i in vcat([2,], 3:R-1)]
    #     pspace_B = [space(B, i) for i in vcat([R,], 1:R-3)]
    #     pspace_A == pspace_B || throw(SpaceMismatch("$(otimes(pspace_A...)) ≠ $(otimes(pspace_B...))"))
    #     return new{R, typeof(A), typeof(B)}(A, B)
    # end
    function TransferMatrix{L, R}(A::AbstractVector{<:Union{MPSTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}, B::AbstractVector{<:Union{AdjointMPSTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}}) where {L, R}
    end
    function TransferMatrix(A::AbstractVector{<:Union{MPSTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}}, B::AbstractVector{<:Union{AdjointMPSTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}}) where R
    end
    function TransferMatrix(A::Union{MPSTensor{R}, LeftIsometricTensor{R}, RightIsometricTensor{R}}, B::Union{AdjointMPSTensor{R}, AdjointLeftIsometricTensor{R}, AdjointRightIsometricTensor{R}}) where R
    end
end
