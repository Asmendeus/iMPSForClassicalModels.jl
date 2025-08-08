"""
    struct RightIsometricTensor{R} <: AbstractMPSTensor{R}
        A::AbstractTensorMap
    end

Wrapper type for R-leg right isometric tensor.

Convention (' marks codomain):

                          3            3 ... (R-1)
                          |            \\ | /
    1' -- AR -- 3   1' -- AR -- 4   1' -- AR -- R
          |               |               |
          2'              2'              2'

Right isometry means

     -- AR  --     |
        ||    |  = I    (identity matrix)
     -- AR* --     |

# Constructors
    RightIsometricTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
    RightIsometricTensor{R}(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct RightIsometricTensor{R} <: AbstractMPSTensor{R}
    A::AbstractTensorMap

    function RightIsometricTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (numout(A) == 2 && numin(A)> 0) || throw(ArgumentError("The space $(space(A)) does not conform to `RightIsometricTensor` space"))
        (check && !isRightIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not right isometric"))
        R = numind(A)
        return new{R}(A)
    end
    function RightIsometricTensor{R}(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true) where R
        (numout(A) == 2 && numin(A)== R-2 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `RightIsometricTensor{$R}` space"))
        (check && !isRightIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not right isometric"))
        return new{R}(A)
    end
end

# ============ Adjoint ============
"""
    struct AdjointRightIsometricTensor{R} <: AbstractMPSTensor{R}
        A::AbstractTensorMap
    end

Lazy wrapper type for R-leg right isometric tensor.

Convention (' marks codomain):

         3               4                   R
         |               |                   |
    2 -- BR -- 1'   3 -- BR -- 2'   (R-1) -- BR -- (R-2)'
                         |                 / | \\
                         1'              1' ... (R-3)'

Right isometry means

     -- BR* --    |
        ||    | = I    (identity matrix)
     -- BR  --    |

# Constructors
    AdjointRightIsometricTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
    AdjointRightIsometricTensor{R}(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct AdjointRightIsometricTensor{R} <: AbstractMPSTensor{R}
    A::AbstractTensorMap

    function AdjointRightIsometricTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (numin(A)== 2 && numout(A) > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointRightIsometricTensor` space"))
        (check && !isRightIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not right isometric"))
        R = numind(A)
        return new{R}(A)
    end
    function AdjointRightIsometricTensor{R}(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true) where R
        (numin(A)== 2 && numout(A) == R-2 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointRightIsometricTensor{$R}` space"))
        (check && !isRightIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not right isometric"))
        return new{R}(A)
    end
end
adjoint(A::RightIsometricTensor) = AdjointRightIsometricTensor(A.A')
adjoint(A::AdjointRightIsometricTensor)::RightIsometricTensor = A.A'
