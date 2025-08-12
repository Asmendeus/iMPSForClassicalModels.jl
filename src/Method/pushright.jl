"""
    function pushright(X₀::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::LocalTensor{R}, B::AdjointLocalTensor{R})
        -> X::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}

Contract transfer matrix environment and right environment tensor
    -- A --       --
       |   |        |
       |   X₀  =    X
       |   |        |
    -- B --       --
"""
function pushright(X₀::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::LocalTensor{R}, B::AdjointLocalTensor{R}) where R
    if R == 2
        @tensor X[-1; -2] := X₀.A[1 2] * A.A[-1 1] * B.A[2 -2]
        return typeof(X₀)(X)
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

"""
    function pushright(X₀::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::AbstractVector{<:LocalTensor{R}}, B::AbstractVector{<:AdjointLocalTensor{R}})
        -> X::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}
    function pushright(X₀::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, t::TransferMatrix{L, R})
        -> X::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}

Iteratively contract transfer matrix and right environment tensor
    -- A[1] -- ... -- A[L] --       --
       |              |      |        |
       |              |      X₀  =    X
       |              |      |        |
    -- B[1] -- ... -- B[L] --       --
"""
function pushright(X₀::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::AbstractVector{<:LocalTensor{R}}, B::AbstractVector{<:AdjointLocalTensor{R}}) where R
    (L = length(A)) == length(B) || throw(ArgumentError("Mismatched lengths: $(length(A)) ≠ $(length(B))"))
    X = deepcopy(X₀)
    for l in L:-1:1
        X = pushright(X, A[l], B[l])
    end
    return X
end
function pushright(X₀::Union{RightEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, t::TransferMatrix{L, R}) where {L, R}
    return pushright(X₀, t.A, t.B)
end

"""
    function pushright(X₀::RightEnvironmentTensor{N}, A::LocalTensor{R}, O::AbstractVector{MPOTensor}, B::AdjointLocalTensor{R})
        -> X::RightEnvironmentTensor{N}
    function pushright(X₀::RightEnvironmentTensor{3}, A::LocalTensor{R}, O::MPOTensor, B::AdjointLocalTensor{R})
        -> X::RightEnvironmentTensor{3}

Contract channel environment and right environment tensor
                __          __
    -- A ----- |  |     -- |  |
       |       |  |        |  |
    —— O[1] —— |  |     —— |  |
       |       |  |        |  |
       ⋮       |X₀|  =    ⋮ | X|
       |       |  |        |  |
    —— O[W] —— |  |     —— |  |
       |       |  |        |  |
    -- B ----- |  |     -- |  |
                ‾‾          ‾‾
"""
function pushright(X₀::RightEnvironmentTensor{N}, A::LocalTensor{R}, O::AbstractVector{MPOTensor}, B::AdjointLocalTensor{R}) where {N, R}
    N == length(O) + 2 || throw(ArgumentError("Mismatched number of RightEnvironmentTensors' legs and MPOTensors"))
    if R == 2
        throw(ArgumentError("Illegal behavior: boundary tensors are `BondTensor` and `AdjointBondTensor`"))
    # ==== Retain to improve performance ===={
    elseif R == 3 && N == 3
        @tensor X[-1 -2; -3] := X₀.A[1 2 3] * A.A[-1 4 1] * O[1][-2 5 4 2] * B.A[3 -3 5]
        return typeof(X₀)(X)
    elseif R == 4 && N == 3
        @tensor X[-1 -2; -3] := X₀.A[1 2 3] * A.A[-1 4 6 1] * O[1][-2 5 4 2] * B.A[6 3 -3 5]
        return typeof(X₀)(X)
    # ==== Retain to improve performance ====}
    else
        # R ≥ 3
        str_legs_X = prod(["-$i " for i in 1:N-1]) * "; -$N"
        str_legs_X₀ = "1" * prod([" $i" for i in 2:N])
        str_legs_A = "-1 $(N+1)" * prod([" $i" for i in 2*N:2*N+R-4]) * " 1"
        str_O = prod(["* O[$w][-$(w+1) $(N+w+1) $(N+w) $(w+1)]" for w in 1:N-2])
        str_legs_B = prod(["$i " for i in 2*N:2*N+R-4]) * "$N -$N $(2*N-1)"
        str_expr = "@tensor X[$str_legs_X] := X₀.A[$str_legs_X₀] * A.A[$str_legs_A] $str_O * B.A[$str_legs_B]"
        eval(Meta.parse(str_expr))
        return typeof(X₀)(X)
    end
end
function pushright(X₀::RightEnvironmentTensor{3}, A::LocalTensor{R}, O::MPOTensor, B::AdjointLocalTensor{R}) where R
    return pushright(X₀, A, [O,], B)
end

"""
    function pushright(X₀::RightEnvironmentTensor{N}, A::AbstractVector{<:LocalTensor{R}}, O::AbstractMatrix{MPOTensor}, B::AbstractVector{<:AdjointLocalTensor{R}})
        -> X::RightEnvironmentTensor{N}
    function pushright(X₀::RightEnvironmentTensor{N}, c::ChannelEnvironment{N, L, R})
        -> X::RightEnvironmentTensor{N}

Contract channel environment tensor and right environment
                                     __          __
    -- A[1] ----- ... -- A[L] ----- |  |     -- |  |
       |                 |          |  |        |  |
    —— O[1, 1] —— ... —— O[1, L] —— |  |     —— |  |
       |                 |          |  |        |  |
       ⋮                  ⋮          |X₀|  =   ⋮ | X|
       |                 |          |  |        |  |
    —— O[W, 1] —— ... —— O[W, L] —— |  |     —— |  |
       |                 |          |  |        |  |
    -- B[1] ----- ... -- B[L] ----- |  |     -- |  |
                                     ‾‾          ‾‾
"""
function pushright(X₀::RightEnvironmentTensor{N}, A::AbstractVector{<:LocalTensor{R}}, O::AbstractMatrix{MPOTensor}, B::AbstractVector{<:AdjointLocalTensor{R}}) where {N, R}
    N == size(O, 1) + 2 || throw(ArgumentError("Mismatched widths: $N ≠ $(size(O, 1)+2)"))
    (L = length(A)) == size(O, 2) == length(B) || throw(ArgumentError("Mismatched lengths: ($L, $(length(A)), $(size(O, 2)), $(length(B)))"))
    X = deepcopy(X₀)
    for l in L:-1:1
        X = pushright(X, A[l], O[:, l], B[l])
    end
    return X
end
function pushright(X₀::RightEnvironmentTensor{N}, c::ChannelEnvironment{N, L, R})
    return pushright(X₀, c.A, c.O, c.B)
end
