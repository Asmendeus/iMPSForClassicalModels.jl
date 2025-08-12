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


# ========== AdjointBondTensor multiplication ==========
# Notice A * B = (A' * B')'   if {typeof(A), typeof(B)} = {AdjointBondTensor, AdjointMPSTensor}
function Base.:(*)(A::AdjointBondTensor, B::AdjointLocalTensor{R}) where R
    if R == 2
        @tensor E[-1; -2] := A.A[1 -2] * B.A[-1 1]
    # ==== Retain to improve performance ===={
    elseif R == 3
        @tensor E[-1; -2 -3] := A.A[1 -2] * B.A[-1 1 -3]
    elseif R == 4
        @tensor E[-1 -2; -3 -4] := A.A[1 -3] * B.A[-1 -2 1 -4]
    # ==== Retain to improve performance ====}
    else
        # R ≥ 3
        str_legs_E = "-1" * prod([" -$i" for i in 2:R-2]) * "; -$(R-1) -$R"
        str_legs_B = prod(["-$i " for i in 1:R-2]) * "1 -$R"
        str_expr = "@tensor E[$str_legs_E] := A.A[1 -$(R-1)] * B.A[$str_legs_B]"
        eval(Meta.parse(str_expr))
    end
    return typeof(B)(E)
end
function Base.:(*)(A::AdjointLocalTensor{R}, B::AdjointBondTensor) where R
    if R == 2
        @tensor E[-1; -2] := A.A[1 -2] * B.A[-1 1]
    # ==== Retain to improve performance ===={
    elseif R == 3
        @tensor E[-1; -2 -3] := A.A[1 -2 -3] * B.A[-1 1]
    elseif R == 4
        @tensor E[-1 -2; -3 -4] := A.A[-1 1 -3 -4] * B.A[-2 1]
    # ==== Retain to improve performance ====}
    else
        # R ≥ 3
        str_legs_E = "-1" * prod([" -$i" for i in 2:R-2]) * "; -$(R-1) -$R"
        str_legs_A = prod(["-$i " for i in 1:R-3]) * "1 -$(R-1) -$R"
        str_expr = "@tensor E[$str_legs_E] := A.A[$str_legs_A] * B.A[-$(R-2) 1]"
        eval(Meta.parse(str_expr))
    end
    return typeof(A)(E)
end
