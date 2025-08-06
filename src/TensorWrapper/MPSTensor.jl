"""
    abstract type AbstractMPSTensor <: AbstractTensorWrapper

Elements of rank-3 MPS.
"""
abstract type AbstractMPSTensor <: AbstractTensorWrapper end

# =================================
"""
    struct MPSTensor <: AbstractMPSTensor
        A::AbstractTensorMap
    end

Wrapper type for rank-3 MPS local tensors.

Convention (' marks codomain):

          2'
          |
    1' —— A —— 3

# Constructor(s)
    MPSTensor(::AbstractTensorMap)
"""
struct MPSTensor <: AbstractMPSTensor
    A::AbstractTensorMap

    function MPSTensor(A::AbstractTensorMap)
        (length(codomain(A)) == 2 && length(domain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `MPSTensor` space"))
        return new(A)
    end
end

"""
    struct LeftIsometricTensor <: AbstractMPSTensor
        A::AbstractTensorMap
    end

Wrapper type for rank-3 MPS local tensor with property of left isometry.

Convention (' marks codomain):

          2'
          |
    1' —— AL —— 3

Left isometry means

     —— AL* ——    |
    |   |       = I    (identity matrix)
     —— AL  ——    |

# Constructor(s)
    LeftIsometricTensor(::AbstractTensorMap)
"""
struct LeftIsometricTensor <: AbstractMPSTensor
    A::AbstractTensorMap

    function LeftIsometricTensor(A::AbstractMPSTensor)
        (length(codomain(A)) == 2 && length(domain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `LeftIsometricTensor` space"))
        return new(A)
    end
end

"""
    struct RightIsometricTensor <: AbstractMPSTensor
        A::AbstractTensorMap
    end

Wrapper type for rank-3 MPS local tensor with property of right isometry.

Convention (' marks codomain):

          2'
          |
    1' —— AR —— 3

Right isometry means

     —— AR* ——    |
        |     | = I    (identity matrix)
     —— AR  ——    |

# Constructor(s)
    RightIsometricTensor(::AbstractTensorMap)
"""
struct RightIsometricTensor <: AbstractMPSTensor
    A::AbstractTensorMap

    function RightIsometricTensor(A::AbstractMPSTensor)
        (length(codomain(A)) == 2 && length(domain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `RightIsometricTensor` space"))
        return new(A)
    end
end

# ============ Adjoint ============
"""
    struct AdjointMPSTensor <: AbstractMPSTensor
        A::AbstractTensorMap
    end

Lazy wrapper type for rank-3 adjoint MPS local tensors.

Convention (' marks codomain):

    2 —— B —— 1'
         |
         3

# Constructor(s)
    AdjointMPSTensor(::AbstractTensorMap)
"""
struct AdjointMPSTensor <: AbstractMPSTensor
    A::AbstractTensorMap

    function AdjointMPSTensor(A::AbstractTensorMap)
        (length(codomain(A)) == 1 && length(domain(A)) == 2) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointMPSTensor` space"))
        return new(A)
    end
end
adjoint(A::MPSTensor) = AdjointMPSTensor(A.A')
adjoint(A::AdjointMPSTensor)::MPSTensor = A.A'

"""
    struct AdjointLeftIsometricTensor <: AbstractMPSTensor
        A::AbstractTensorMap
    end

Lazy wrapper type for rank-3 adjoint MPS local tensors with property of left isometry.

Convention (' marks codomain):

    2 —— BL —— 1'
         |
         3

Left isometry means

     —— BL  ——    |
    |   |      = I    (identity matrix)
     —— BL* ——    |

# Constructor(s)
    AdjointLeftIsometricTensor(::AbstractTensorMap)
"""
struct AdjointLeftIsometricTensor <: AbstractMPSTensor
    A::AbstractTensorMap

    function AdjointLeftIsometricTensor(A::AbstractMPSTensor)
        (length(codomain(A)) == 1 && length(domain(A)) == 2) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointLeftIsometricTensor` space"))
        return new(A)
    end
end
adjoint(A::LeftIsometricTensor) = AdjointLeftIsometricTensor(A.A')
adjoint(A::AdjointLeftIsometricTensor)::LeftIsometricTensor = A.A'

"""
    struct AdjointRightIsometricTensor <: AbstractMPSTensor
        A::AbstractTensorMap
    end

Lazy wrapper type for rank-3 adjoint MPS local tensors with property of right isometry.

Convention (' marks codomain):

    2 —— BR —— 1'
         |
         3

Right isometry means

     —— BR  ——    |
        |     | = I    (identity matrix)
     —— BR* ——    |

# Constructor(s)
    AdjointRightIsometricTensor(::AbstractTensorMap)
"""
struct AdjointRightIsometricTensor <: AbstractMPSTensor
    A::AbstractTensorMap

    function AdjointRightIsometricTensor(A::AbstractMPSTensor)
        (length(codomain(A)) == 1 && length(domain(A)) == 2) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointRightIsometricTensor` space"))
        return new(A)
    end
end
adjoint(A::RightIsometricTensor) = AdjointRightIsometricTensor(A.A')
adjoint(A::AdjointRightIsometricTensor)::RightIsometricTensor = A.A'
