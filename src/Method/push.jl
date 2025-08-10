# For TransferMatrix
"""
    function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::LocalTensor{R}, B::AdjointLocalTensor{R}) -> X::LeftEnvironmentTensor{2}

Contract left boundary tensor and transfer matrix left environment
     -- A --      --
    |   |        |
    X₀  |     =  X
    |   |        |
     -- B --      --
"""
function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::LocalTensor{R}, B::AdjointLocalTensor{R}) where R
    if R == 2
        throw(ArgumentError("`TransferMatrix{L, 2}` is illegal."))
    # ==== Retain to improve performance ===={
    elseif R == 3
        @tensor X[-1; -2] := X₀.A[1 2] * A.A[2 3 -2] * B.A[-1 1 3]
        return typeof(X₀)(X)
    elseif R == 4
        @tensor X[-1; -2] := X₀.A[1 2] * A.A[2 3 4 -2] * B.A[4 -1 1 3]
        return typeof(X₀)(X)
    # ==== Retain to improve performance ====}
    else
        # R ≥ 3
        str_legs_A = "2 3" * prod([" $i" for i in 4:R]) * " -2"
        str_legs_B = prod(["$i " for i in 4:R]) * "-1 1 3"
        str_expr = "@tensor X[-1; -2] := X₀.A[1 2] * A.A[$str_legs_A] * B.A[$str_legs_B]"
        eval(Meta.parse(str_expr))
        return typeof(X₀)(X)
    end
end

"""
    function pushright(X₀::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::LocalTensor{R}, B::AdjointLocalTensor{R}) -> X::RightEnvironmentTensor{2}

Contract transfer matrix right environment and right boundary tensor
    -- A --       --
       |   |        |
       |   X₀  =    X
       |   |        |
    -- B --       --
"""
function pushright(X₀::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::LocalTensor{R}, B::AdjointLocalTensor{R}) where R
    if R == 2
        throw(ArgumentError("`TransferMatrix{L, 2}` is illegal."))
    # ==== Retain to improve performance ===={
    elseif R == 3
        @tensor X[-1; -2] := X₀.A[2 3] * A.A[-1 1 2] * B.A[3 -2 1]
        return typeof(X₀)(X)
    elseif R == 4
        @tensor X[-1; -2] := X₀.A[2 3] * A.A[-1 1 4 2] * B.A[4 3 -2 1]
        return typeof(X₀)(X)
    # ==== Retain to improve performance ====}
    else
        # R ≥ 3
        str_legs_A = "-1 1" * prod([" $i" for i in 4:R]) * " 2"
        str_legs_B = prod(["$i " for i in 4:R]) * "3 -2 1"
        str_expr = "@tensor X[-1; -2] := X₀.A[2 3] * A.A[$str_legs_A] * B.A[$str_legs_B]"
        eval(Meta.parse(str_expr))
        return typeof(X₀)(X)
    end
end

# For Environments
