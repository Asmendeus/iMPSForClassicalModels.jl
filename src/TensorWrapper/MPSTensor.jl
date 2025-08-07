"""
    abstract type AbstractMPSTensor <: AbstractTensorWrapper

Elements both of MPS and MPO.
"""
abstract type AbstractMPSTensor <: AbstractTensorWrapper end

# =================================
"""
    struct MPSTensor{R} <: AbstractMPSTensor
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
struct MPSTensor{R} <: AbstractMPSTensor
    A::AbstractTensorMap

    function MPSTensor(A::AbstractTensorMap)
        if length(codomain(A)) == 1 && length(domain(A)) == 1
            return new{2}(A)
        elseif length(codomain(A)) == 2 && length(domain(A)) > 0
            R = 2 + length(domain(A))
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `MPSTensor` space"))
        end
    end
    function MPSTensor{R}(A::AbstractTensorMap)
        if R == 2
            (length(codomain(A)) == 1 && length(domain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `MPSTensor{$R}` space"))
        elseif R > 2
            (length(codomain(A)) == 2 && length(domain(A)) == R-2 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `MPSTensor{$R}` space"))
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `MPSTensor{$R}` space"))
        end
        return new{R}(A)
    end
end

# ============ Adjoint ============
"""
    struct AdjointMPSTensor{R} <: AbstractMPSTensor
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
struct AdjointMPSTensor{R} <: AbstractMPSTensor
    A::AbstractTensorMap

    function AdjointMPSTensor(A::AbstractTensorMap)
        if length(domain(A)) == 1 && length(codomain(A)) == 1
            return new{2}(A)
        elseif length(domain(A)) == 2 && length(codomain(A)) > 0
            R = 2 + length(codomain(A))
            return new{R}(A)
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `AdjointMPSTensor` space"))
        end
    end
    function AdjointMPSTensor{R}(A::AbstractTensorMap)
        if R == 2
            (length(domain(A)) == 1 && length(codomain(A)) == 1) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointMPSTensor{$R}` space"))
        elseif R > 2
            (length(domain(A)) == 2 && length(codomain(A)) == R-2 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `AdjointMPSTensor{$R}` space"))
        else
            throw(ArgumentError("The space $(space(A)) does not conform to `AdjointMPSTensor{$R}` space"))
        end
        return new{R}(A)
    end
end
adjoint(A::MPSTensor) = AdjointMPSTensor(A.A')
adjoint(A::AdjointMPSTensor)::MPSTensor = A.A'
