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
