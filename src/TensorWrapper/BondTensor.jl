"""
    abstract type AbstractBondTensor <: AbstractTensorWrapper

Bond tensor (matrix) between two MPS tensors.
"""
abstract type AbstractBondTensor <: AbstractTensorWrapper end

# =================================
"""
    struct BondTensor <: AbstractBondTensor
        A::AbstractTensorMap
    end

Wrapper type for rank-2 bond tensors.

Convention (' marks codomain):

    1' —— C —— 2

# Constructor(s)
    BondTensor(::AbstractTensorMap)
"""
struct BondTensor <: AbstractBondTensor
    A::AbstractTensorMap

    function BondTensor(A::AbstractTensorMap)
        (length(codomain(A)) == 1 && length(domain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `BondTensor` space"))
        return new(A)
    end
end

# ============ Adjoint ============
"""
    struct AdjointBondTensor <: AbstractBondTensor
        A::AbstractTensorMap
    end

Lazy wrapper type for rank-2 adjoint bond tensors.

Convention (' marks codomain):

    2 —— C —— 1'

# Constructor(s)
    AdjointBondTensor(::AbstractTensorMap)
"""
struct AdjointBondTensor <: AbstractBondTensor
    A::AbstractTensorMap

    function AdjointBondTensor(A::AbstractTensorMap)
        (length(codomain(A)) == 1 && length(domain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointBondTensor` space"))
        return new(A)
    end
end
adjoint(A::BondTensor) = AdjointBondTensor(A.A')
adjoint(A::AdjointBondTensor)::BondTensor = A.A'
