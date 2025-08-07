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
    struct LeftIsometricMPSTensor <: AbstractMPSTensor
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
    LeftIsometricMPSTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct LeftIsometricMPSTensor <: AbstractMPSTensor
    A::AbstractTensorMap

    function LeftIsometricMPSTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (length(codomain(A)) == 2 && length(domain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `LeftIsometricMPSTensor` space"))
        (check && !isLeftIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not left isometric"))
        return new(A)
    end
end

"""
    struct RightIsometricMPSTensor <: AbstractMPSTensor
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
    RightIsometricMPSTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct RightIsometricMPSTensor <: AbstractMPSTensor
    A::AbstractTensorMap

    function RightIsometricMPSTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (length(codomain(A)) == 2 && length(domain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `RightIsometricMPSTensor` space"))
        (check && !isRightIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not right isometric"))
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
    struct AdjointLeftIsometricMPSTensor <: AbstractMPSTensor
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
    AdjointLeftIsometricMPSTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct AdjointLeftIsometricMPSTensor <: AbstractMPSTensor
    A::AbstractTensorMap

    function AdjointLeftIsometricMPSTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (length(codomain(A)) == 1 && length(domain(A)) == 2) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointLeftIsometricMPSTensor` space"))
        (check && !isLeftIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not left isometric"))
        return new(A)
    end
end
adjoint(A::LeftIsometricMPSTensor) = AdjointLeftIsometricMPSTensor(A.A')
adjoint(A::AdjointLeftIsometricMPSTensor)::LeftIsometricMPSTensor = A.A'

"""
    struct AdjointRightIsometricMPSTensor <: AbstractMPSTensor
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
    AdjointRightIsometricMPSTensor(::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
"""
struct AdjointRightIsometricMPSTensor <: AbstractMPSTensor
    A::AbstractTensorMap

    function AdjointRightIsometricMPSTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (length(codomain(A)) == 1 && length(domain(A)) == 2) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointRightIsometricMPSTensor` space"))
        (check && !isRightIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not right isometric"))
        return new(A)
    end
end
adjoint(A::RightIsometricMPSTensor) = AdjointRightIsometricMPSTensor(A.A')
adjoint(A::AdjointRightIsometricMPSTensor)::RightIsometricMPSTensor = A.A'
