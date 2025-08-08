"""
    struct RightIsometricTensor{R} <: AbstractLocalTensor{R}
        A::AbstractTensorMap
    end

Wrapper type for R-leg right isometric tensor.

Convention (' marks codomain):

                                          3            3 ... (R-1)
                                          |            \\ | /
    1' -- AR -- 2   1' -- AR -- 3   1' -- AR -- 4   1' -- AR -- R
                          |               |               |
                          2'              2'              2'

Right isometry means

     -- AR  --     |    -- AR  --     |     -- AR  --     |
              |  = I       |     |  = I        ||    |  = I    (identity matrix)
     -- AR* --     |    -- AR* --     |     -- AR* --     |

# Constructors
    RightIsometricTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
    RightIsometricTensor{R}(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct RightIsometricTensor{R} <: AbstractLocalTensor{R}
    A::AbstractTensorMap

    function RightIsometricTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (check && !isRightIsometric(A, false; tol=tol)) && throw(ArgumentError("The tensor is not right isometric"))
        if numout(A) == 1 && numin(A) == 1
            return new{2}(A)
        elseif numout(A) == 2 && numin(A) > 0
            R = numind(A)
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `RightIsometricTensor` space"))
        end
    end
    function RightIsometricTensor{R}(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true) where R
        (check && !isRightIsometric(A, false; tol=tol)) && throw(ArgumentError("The tensor is not right isometric"))
        if numout(A) == 1 && numin(A) == 1
            return new{2}(A)
        elseif numout(A) == 2 && numin(A) == R-2 > 0
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `RightIsometricTensor{$R}` space"))
        end
    end
end

const RightIsometricBondTensor = RightIsometricTensor{2}
const RightIsometricMPSTensor = RightIsometricTensor{3}
const RightIsometricMPOTensor = RightIsometricTensor{4}
const RightIsometricMPSOrMPOTensor = Union{RightIsometricMPSTensor, RightIsometricMPOTensor}

# ============ Adjoint ============
"""
    struct AdjointRightIsometricTensor{R} <: AbstractLocalTensor{R}
        A::AbstractTensorMap
    end

Lazy wrapper type for R-leg right isometric tensor.

Convention (' marks codomain):

                         3               4                   R
                         |               |                   |
    2 -- BR -- 1'   2 -- BR -- 1'   3 -- BR -- 2'   (R-1) -- BR -- (R-2)'
                                         |                 / | \\
                                         1'              1' ... (R-3)'

Right isometry means

     -- BR* --     |    -- BR* --     |     -- BR* --     |
              |  = I       |     |  = I        ||    |  = I    (identity matrix)
     -- BR  --     |    -- BR  --     |     -- BR  --     |

# Constructors
    AdjointRightIsometricTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
    AdjointRightIsometricTensor{R}(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct AdjointRightIsometricTensor{R} <: AbstractLocalTensor{R}
    A::AbstractTensorMap

    function AdjointRightIsometricTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (check && !isRightIsometric(A, true; tol=tol)) && throw(ArgumentError("The tensor is not right isometric"))
        if numin(A) == 1 && numout(A) == 1
            return new{2}(A)
        elseif numin(A)== 2 && numout(A) > 0
            R = numind(A)
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `AdjointRightIsometricTensor` space"))
        end
    end
    function AdjointRightIsometricTensor{R}(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true) where R
        (check && !isRightIsometric(A, true; tol=tol)) && throw(ArgumentError("The tensor is not right isometric"))
        if numin(A) == 1 && numout(A) == 1
            return new{2}(A)
        elseif numin(A)== 2 && numout(A) == R-2 > 0
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `AdjointRightIsometricTensor{$R}` space"))
        end
    end
end
adjoint(A::RightIsometricTensor) = AdjointRightIsometricTensor(A.A')
adjoint(A::AdjointRightIsometricTensor) = RightIsometricTensor(A.A')

const AdjointRightIsometricBondTensor = AdjointRightIsometricTensor{2}
const AdjointRightIsometricMPSTensor = AdjointRightIsometricTensor{3}
const AdjointRightIsometricMPOTensor = AdjointRightIsometricTensor{4}
const AdjointRightIsometricMPSOrMPOTensor = Union{AdjointRightIsometricMPSTensor, AdjointRightIsometricMPOTensor}
