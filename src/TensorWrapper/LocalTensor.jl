"""
    abstract type AbstractLocalTensor{R} <: AbstractTensorWrapper{R}

R-leg local tensor.
"""
abstract type AbstractLocalTensor{R} <: AbstractTensorWrapper{R} end

const AbstractBondTensor = AbstractLocalTensor{2}
const AbstractMPSTensor = AbstractLocalTensor{3}
const AbstractMPOTensor = AbstractLocalTensor{4}

# =================================
"""
    struct LocalTensor{R} <: AbstractLocalTensor{R}
        A::AbstractTensorMap
    end

Wrapper type for R-leg local tensor.

Convention (' marks codomain):
                                          3            3 ... (R-1)
                                          |            \\ | /
    1' -- A -- 2    1' -- A -- 3    1' -- A -- 4    1' -- A -- R
                          |               |               |
                          2'              2'              2'

# Constructors
    LocalTensor(::AbstractTensorMap)
    LocalTensor{R}(::AbstractTensorMap)
"""
struct LocalTensor{R} <: AbstractLocalTensor{R}
    A::AbstractTensorMap

    function LocalTensor(A::AbstractTensorMap)
        if numout(A) == 1 && numin(A) == 1
            return new{2}(A)
        elseif numout(A) == 2 && numin(A) > 0
            R = numind(A)
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `LocalTensor` space"))
        end
    end
    function LocalTensor{R}(A::AbstractTensorMap) where R
        if R == 2
            (numout(A) == 1 && numin(A) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `LocalTensor{$R}` space"))
        elseif R > 2
            (numout(A) == 2 && numin(A) == R-2 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `LocalTensor{$R}` space"))
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `LocalTensor{$R}` space"))
        end
        return new{R}(A)
    end
end

const BondTensor = LocalTensor{2}
const MPSTensor = LocalTensor{3}
const MPOTensor = LocalTensor{4}

# ============ Adjoint ============
"""
    struct AdjointLocalTensor{R} <: AbstractLocalTensor{R}
        A::AbstractTensorMap
    end

Lazy wrapper type for adjoint R-leg local tensor.

Convention (' marks codomain):

                         3               4                   R
                         |               |                   |
    2 -- B -- 1'    2 -- B -- 1'    3 -- B -- 2'    (R-1) -- B  -- (R-2)'
                                         |                 / | \\
                                         1'              1' ... (R-3)'

# Constructors
    AdjointLocalTensor(::AbstractTensorMap)
    AdjointLocalTensor{R}(::AbstractTensorMap)
"""
struct AdjointLocalTensor{R} <: AbstractLocalTensor{R}
    A::AbstractTensorMap

    function AdjointLocalTensor(A::AbstractTensorMap)
        if numin(A) == 1 && numout(A) == 1
            return new{2}(A)
        elseif numin(A) == 2 && numout(A) > 0
            R = numind(A)
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `AdjointLocalTensor` space"))
        end
    end
    function AdjointLocalTensor{R}(A::AbstractTensorMap) where R
        if R == 2
            (numin(A) == 1 && numout(A) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointLocalTensor{$R}` space"))
        elseif R > 2
            (numin(A) == 2 && numout(A) == R-2 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointLocalTensor{$R}` space"))
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `AdjointLocalTensor{$R}` space"))
        end
        return new{R}(A)
    end
end
adjoint(A::LocalTensor) = AdjointLocalTensor(A.A')
adjoint(A::AdjointLocalTensor)::LocalTensor = A.A'

const AdjointBondTensor = AdjointLocalTensor{2}
const AdjointMPSTensor = AdjointLocalTensor{3}
const AdjointMPOTensor = AdjointLocalTensor{4}
