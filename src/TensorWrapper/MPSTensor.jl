"""
    abstract type AbstractMPSTensor{R} <: AbstractTensorWrapper

R-leg local tensor both of MPS and MPO.
"""
abstract type AbstractMPSTensor{R} <: AbstractTensorWrapper end
numind(::AbstractMPSTensor{R}) where R = R

# =================================
"""
    struct MPSTensor{R} <: AbstractMPSTensor{R}
        A::AbstractTensorMap
    end

Wrapper type for R-leg local tensor.
In particular, R == 2 for bond tensor.

Convention (' marks codomain):
                                          3            3 ... (R-1)
                                          |            \\ | /
    1' -- A -- 2    1' -- A -- 3    1' -- A -- 4    1' -- A -- R
                          |               |               |
                          2'              2'              2'

# Constructors
    MPSTensor(::AbstractTensorMap)
    MPSTensor{R}(::AbstractTensorMap)
"""
struct MPSTensor{R} <: AbstractMPSTensor{R}
    A::AbstractTensorMap

    function MPSTensor(A::AbstractTensorMap)
        if numout(A) == 1 && numin(A) == 1
            return new{2}(A)
        elseif numout(A) == 2 && numin(A) > 0
            R = numind(A)
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `MPSTensor` space"))
        end
    end
    function MPSTensor{R}(A::AbstractTensorMap) where R
        if R == 2
            (numout(A) == 1 && numin(A) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `MPSTensor{$R}` space"))
        elseif R > 2
            (numout(A) == 2 && numin(A) == R-2 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `MPSTensor{$R}` space"))
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `MPSTensor{$R}` space"))
        end
        return new{R}(A)
    end
end

# ============ Adjoint ============
"""
    struct AdjointMPSTensor{R} <: AbstractMPSTensor{R}
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
    AdjointMPSTensor(::AbstractTensorMap)
    AdjointMPSTensor{R}(::AbstractTensorMap)
"""
struct AdjointMPSTensor{R} <: AbstractMPSTensor{R}
    A::AbstractTensorMap

    function AdjointMPSTensor(A::AbstractTensorMap)
        if numin(A) == 1 && numout(A) == 1
            return new{2}(A)
        elseif numin(A) == 2 && numout(A) > 0
            R = numind(A)
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `AdjointMPSTensor` space"))
        end
    end
    function AdjointMPSTensor{R}(A::AbstractTensorMap) where R
        if R == 2
            (numin(A) == 1 && numout(A) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointMPSTensor{$R}` space"))
        elseif R > 2
            (numin(A) == 2 && numout(A) == R-2 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointMPSTensor{$R}` space"))
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `AdjointMPSTensor{$R}` space"))
        end
        return new{R}(A)
    end
end
adjoint(A::MPSTensor) = AdjointMPSTensor(A.A')
adjoint(A::AdjointMPSTensor)::MPSTensor = A.A'
