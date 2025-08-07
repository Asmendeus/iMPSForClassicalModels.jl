"""
    struct LeftIsometricTensor{R} <: AbstractMPSTensor
        A::AbstractTensorMap
    end

Wrapper type for R-leg left isometric tensor.

Convention (' marks codomain):

                          3            3 ... (R-1)
                          |            \\ | /
    1' -- AL -- 3   1' -- AL -- 4   1' -- AL -- R
          |               |               |
          2'              2'              2'

Left isometry means

     -- AL  --   |
    |   ||     = I    (identity matrix)
     -- AL* --   |

# Constructors
    LeftIsometricTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
    LeftIsometricTensor{R}(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct LeftIsometricTensor{R} <: AbstractMPSTensor
    A::AbstractTensorMap

    function LeftIsometricTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (length(codomain(A)) == 2 && length(domain(A)) > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `LeftIsometricTensor` space"))
        (check && !isLeftIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not left isometric"))
        R = 2 + length(domain(A))
        return new{R}(A)
    end
    function LeftIsometricTensor{R}(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (length(codomain(A)) == 2 && length(domain(A)) == R-2 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `LeftIsometricTensor{$R}` space"))
        (check && !isLeftIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not left isometric"))
        return new{R}(A)
    end
end

# ============ Adjoint ============
"""
    struct AdjointLeftIsometricTensor{R} <: AbstractMPSTensor
        A::AbstractTensorMap
    end

Lazy wrapper type for adjoint R-leg left isometric tensor.

Convention (' marks codomain):

         3               4                   R
         |               |                   |
    2 -- BL -- 1'   3 -- BL -- 2'   (R-1) -- BL -- (R-2)'
                         |                 / | \\
                         1'              1' ... (R-3)'

Left isometry means

     —— BL* ——   |
    |   ||     = I    (identity matrix)
     —— BL  ——   |

# Constructor(s)
    AdjointLeftIsometricTensor(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
    AdjointLeftIsometricTensor{R}(::AbstractTensorMap; tol::Float64=1e-8, check::Bool=true)
"""
struct AdjointLeftIsometricTensor{R} <: AbstractMPSTensor
    A::AbstractTensorMap

    function AdjointLeftIsometricTensor(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (length(domain(A)) == 2 && length(codomain(A)) > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointLeftIsometricTensor` space"))
        (check && !isLeftIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not left isometric"))
        R = 2 + length(codomain(A))
        return new{R}(A)
    end
    function AdjointLeftIsometricTensor{R}(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm, check::Bool=true)
        (length(domain(A)) == 2 && length(codomain(A)) == R-2 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointLeftIsometricTensor{$R}` space"))
        (check && !isLeftIsometric(A; tol=tol)) && throw(ArgumentError("The tensor is not left isometric"))
        return new{R}(A)
    end
end
adjoint(A::LeftIsometricTensor) = AdjointLeftIsometricTensor(A.A')
adjoint(A::AdjointLeftIsometricTensor)::LeftIsometricTensor = A.A'
