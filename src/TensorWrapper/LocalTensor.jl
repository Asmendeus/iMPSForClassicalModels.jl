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


# ========== BondTensor multiplication ==========
# Notice A * B = (A' * B')'   if {typeof(A), typeof(B)} = {BondTensor, MPSTensor}
function Base.:(*)(A::BondTensor, B::LocalTensor{R}) where R
    if R == 2
        @tensor E[-1; -2] := A.A[-1 1] * B.A[1 -2]
    # ==== Retain to improve performance ===={
    elseif R == 3
        @tensor E[-1 -2; -3] := A.A[-1 1] * B.A[1 -2 -3]
    elseif R == 4
        @tensor E[-1 -2; -3 -4] := A.A[-1 1] * B.A[1 -2 -3 -4]
    # ==== Retain to improve performance ====}
    else
        # R ≥ 3
        str_domain_E = prod([" -$i" for i in 3:R])
        str_legs = "1" * prod([" -$i" for i in 2:R])
        str_expr = "@tensor E[-1 -2;$str_domain_E] := A.A[-1 1] * B.A[$str_legs]"
        eval(Meta.parse(str_expr))
    end
    return typeof(B)(E)
end
function Base.:(*)(A::LocalTensor{R}, B::BondTensor) where R
    if R == 2
        @tensor E[-1; -2] := A.A[-1 1] * B.A[1 -2]
    # ==== Retain to improve performance ===={
    elseif R == 3
        @tensor E[-1 -2; -3] := A.A[-1 -2 1] * B.A[1 -3]
    elseif R == 4
        @tensor E[-1 -2; -3 -4] := A.A[-1 -2 -3 1] * B.A[1 -4]
    # ==== Retain to improve performance ====}
    else
        # R ≥ 3
        str_domain_E = prod([" -$i" for i in 3:R])
        str_legs = prod(["-$i " for i in 1:R-1]) * " 1"
        str_expr = "@tensor E[-1 -2;$str_domain_E] := A.A[$str_legs] * B.A[1 -$R]"
        eval(Meta.parse(str_expr))
    end
    return typeof(A)(E)
end
