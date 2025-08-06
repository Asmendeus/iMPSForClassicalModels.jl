"""
    abstract type AbstractMPSBondTensor <: AbstractTensorWrapper

Bond tensor (matrix) between two MPS tensors.
"""
abstract type AbstractMPSBondTensor <: AbstractTensorWrapper end

# =================================
"""
    struct MPSBondTensor <: AbstractMPSBondTensor
        A::AbstractTensorMap
    end

Wrapper type for rank-2 bond tensors.

Convention (' marks codomain):

    1' —— C —— 2

# Constructor(s)
    MPSBondTensor(::AbstractTensorMap)
"""
struct MPSBondTensor <: AbstractMPSBondTensor
    A::AbstractTensorMap

    function MPSBondTensor(A::AbstractTensorMap)
        (length(codomain(A)) == 1 && length(domain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `MPSBondTensor` space"))
        return new(A)
    end
end

# ============ Adjoint ============
"""
    struct AdjointMPSBondTensor <: AbstractMPSBondTensor
        A::AbstractTensorMap
    end

Lazy wrapper type for rank-2 adjoint bond tensors.

Convention (' marks codomain):

    2 —— C —— 1'

# Constructor(s)
    AdjointMPSBondTensor(::AbstractTensorMap)
"""
struct AdjointMPSBondTensor <: AbstractMPSBondTensor
    A::AbstractTensorMap

    function AdjointMPSBondTensor(A::AbstractTensorMap)
        (length(codomain(A)) == 1 && length(domain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointMPSBondTensor` space"))
        return new(A)
    end
end
adjoint(A::MPSBondTensor) = AdjointMPSBondTensor(A.A')
adjoint(A::AdjointMPSBondTensor)::MPSBondTensor = A.A'
