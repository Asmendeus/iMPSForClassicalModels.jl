"""
    struct LeftIsometricTensor{R} <: AbstractLocalTensor{R}
        A::AbstractTensorMap
    end

Wrapper type for R-leg left isometric tensor.

Convention (' marks codomain):

                                          3            3 ... (R-1)
                                          |            \\ | /
    1' -- AL -- 2   1' -- AL -- 3   1' -- AL -- 4   1' -- AL -- R
                          |               |               |
                          2'              2'              2'

Left isometry means

     -- AL  --   |       -- AL  --   |       -- AL  --   |
    |          = I      |   |      = I      |   ||     = I    (identity matrix)
     -- AL* --   |       -- AL* --   |       -- AL* --   |

# Constructors
    LeftIsometricTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
    LeftIsometricTensor{R}(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct LeftIsometricTensor{R} <: AbstractLocalTensor{R}
    A::AbstractTensorMap

    function LeftIsometricTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (check && !isLeftIsometric(A, false; tol=tol)) && throw(ArgumentError("The tensor is not left isometric"))
        if numout(A) == 1 && numin(A) == 1
            return new{2}(A)
        elseif numout(A) == 2 && numin(A) > 0
            R = numind(A)
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `LeftIsometricTensor` space"))
        end
    end
    function LeftIsometricTensor{R}(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true) where R
        (check && !isLeftIsometric(A, false; tol=tol)) && throw(ArgumentError("The tensor is not left isometric"))
        if numout(A) == 1 && numin(A) == 1
            return new{2}(A)
        elseif (numout(A) == 2 && numin(A) == R-2 > 0)
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `LeftIsometricTensor{$R}` space"))
        end
    end
end

const LeftIsometricBondTensor = LeftIsometricTensor{2}
const LeftIsometricMPSTensor = LeftIsometricTensor{3}
const LeftIsometricMPOTensor = LeftIsometricTensor{4}
const LeftIsometricMPSOrMPOTensor = Union{LeftIsometricMPSTensor, LeftIsometricMPOTensor}

# ============ Adjoint ============
"""
    struct AdjointLeftIsometricTensor{R} <: AbstractLocalTensor{R}
        A::AbstractTensorMap
    end

Lazy wrapper type for adjoint R-leg left isometric tensor.

Convention (' marks codomain):

                         3               4                   R
                         |               |                   |
    2 -- BL -- 1'   2 -- BL -- 1'   3 -- BL -- 2'   (R-1) -- BL -- (R-2)'
                                         |                 / | \\
                                         1'              1' ... (R-3)'

Left isometry means

     -- BL* --   |       -- BL* --   |       -- BL* --   |
    |          = I      |   |      = I      |   ||     = I    (identity matrix)
     -- BL  --   |       -- BL  --   |       -- BL  --   |

# Constructors
    AdjointLeftIsometricTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
    AdjointLeftIsometricTensor{R}(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct AdjointLeftIsometricTensor{R} <: AbstractLocalTensor{R}
    A::AbstractTensorMap

    function AdjointLeftIsometricTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (check && !isLeftIsometric(A, true; tol=tol)) && throw(ArgumentError("The tensor is not left isometric"))
        if numin(A) == 1 && numout(A) == 1
            return new{2}(A)
        elseif numin(A) == 2 && numout(A) > 0
            R = numind(A)
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `AdjointLeftIsometricTensor` space"))
        end
    end
    function AdjointLeftIsometricTensor{R}(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true) where R
        (check && !isLeftIsometric(A, true; tol=tol)) && throw(ArgumentError("The tensor is not left isometric"))
        if numin(A) == 1 && numout(A) == 1
            return new{2}(A)
        elseif numin(A) == 2 && numout(A) == R-2 > 0
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `AdjointLeftIsometricTensor{$R}` space"))
        end
    end
end
adjoint(A::LeftIsometricTensor) = AdjointLeftIsometricTensor(A.A')
adjoint(A::AdjointLeftIsometricTensor) = LeftIsometricTensor(A.A')

const AdjointLeftIsometricBondTensor = AdjointLeftIsometricTensor{2}
const AdjointLeftIsometricMPSTensor = AdjointLeftIsometricTensor{3}
const AdjointLeftIsometricMPOTensor = AdjointLeftIsometricTensor{4}
const AdjointLeftIsometricMPSOrMPOTensor = Union{AdjointLeftIsometricMPSTensor, AdjointLeftIsometricMPOTensor}
