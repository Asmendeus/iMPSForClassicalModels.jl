"""
    function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::LocalTensor{R}, B::AdjointLocalTensor{R})
        -> X::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}

Contract left environment tensor and transfer matrix environment
     -- A --      --
    |   |        |
    X₀  |     =  X
    |   |        |
     -- B --      --
"""
function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::LocalTensor{R}, B::AdjointLocalTensor{R}) where R
    if R == 2
        @tensor X[-1; -2] := X₀.A[1 2] * A.A[2 -2] * B.A[-1 1]
        return typeof(X₀)(X)
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
    function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::AbstractVector{<:LocalTensor{R}}, B::AbstractVector{<:AdjointLocalTensor{R}})
        -> X::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}
    function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, t::TransferMatrix{L, R})
        -> X::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}

Iteratively contract left environment tensor and transfer matrix
     -- A[1] -- ... -- A[L] --      --
    |   |              |           |
    X₀  |              |        =  X
    |   |              |           |
     -- B[1] -- ... -- B[L] --      --
"""
function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::AbstractVector{<:LocalTensor{R}}, B::AbstractVector{<:AdjointLocalTensor{R}}) where R
    (L = length(A)) == length(B) || throw(ArgumentError("Mismatched lengths: $(length(A)) ≠ $(length(B))"))
    X = deepcopy(X₀)
    for l in 1:L
        X = pushleft(X, A[l], B[l])
    end
    return X
end
function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, t::TransferMatrix{L, R}) where {L, R}
    return pushleft(X₀, t.A, t.B)
end

"""
    function pushleft(X₀::LeftEnvironmentTensor{N}, A::LocalTensor{R}, O::AbstractVector{MPOTensor}, B::AdjointLocalTensor{R})
        -> X::LeftEnvironmentTensor{N}
    function pushleft(X₀::LeftEnvironmentTensor{3}, A::LocalTensor{R}, O::MPOTensor, B::AdjointLocalTensor{R})
        -> X::LeftEnvironmentTensor{3}

Contract left environment tensor and channel environment
     __                  __
    |  | -- A -----     |  | --
    |  |    |           |  |
    |  | —— O[1] ——     |  | ——
    |  |    |           |  |
    |X₀|    ⋮        =  | X| ⋮
    |  |    |           |  |
    |  | —— O[W] ——     |  | ——
    |  |    |           |  |
    |  | -- B -----     |  | --
     ‾‾                  ‾‾
"""
function pushleft(X₀::LeftEnvironmentTensor{N}, A::LocalTensor{R}, O::AbstractVector{MPOTensor}, B::AdjointLocalTensor{R}) where {N, R}
    N == length(O) + 2 || throw(ArgumentError("Mismatched number of LeftEnvironmentTensors' legs and MPOTensors"))
    if R == 2
        throw(ArgumentError("Illegal behavior: boundary tensors are `BondTensor` and `AdjointBondTensor`"))
    # ==== Retain to improve performance ===={
    elseif R == 3 && N == 3
        @tensor X[-1; -2 -3] := X₀.A[1 2 3] * A.A[3 5 -3] * O[1][2 4 5 -2] * B.A[-1 1 4]
        return typeof(X₀)(X)
    elseif R == 4 && N == 3
        @tensor X[-1; -2 -3] := X₀.A[1 2 3] * A.A[3 5 6 -3] * O[1][2 4 5 -2] * B.A[6 -1 1 4]
        return typeof(X₀)(X)
    # ==== Retain to improve performance ====}
    else
        # R ≥ 3
        str_legs_X = "-1;" * prod([" -$i" for i in 2:N])
        str_legs_X₀ = "1" * prod([" $i" for i in 2:N])
        str_legs_A = "$N $(2*N-1)" * prod([" $i" for i in 2*N:2*N+R-4]) * " -$N"
        str_O = prod(["* O[$w][$(N-w) $(2*N-w-1) $(2*N-w) -$(N-w)]" for w in 1:N-2])
        str_legs_B = prod(["$i " for i in 2*N:2*N+R-4]) * "-1 1 $(N+1)"
        str_expr = "@tensor X[$str_legs_X] := X₀.A[$str_legs_X₀] * A.A[$str_legs_A] $str_O * B.A[$str_legs_B]"
        eval(Meta.parse(str_expr))
        return typeof(X₀)(X)
    end
end
function pushleft(X₀::LeftEnvironmentTensor{3}, A::LocalTensor{R}, O::MPOTensor, B::AdjointLocalTensor{R}) where R
    return pushleft(X₀, A, [O,], B)
end

"""
    function pushleft(X₀::LeftEnvironmentTensor{N}, A::AbstractVector{<:LocalTensor{R}}, O::AbstractMatrix{MPOTensor}, B::AbstractVector{<:AdjointLocalTensor{R}})
        -> X::LeftEnvironmentTensor{N}
    function pushleft(X₀::LeftEnvironmentTensor{N}, c::ChannelEnvironment{N, L, R})
        -> X::LeftEnvironmentTensor{N}

Contract left environment tensor and channel environment
     __                                       __
    |  | -- A[1] ----- ... -- A[L] -----     |  | --
    |  |    |                 |              |  |
    |  | —— O[1, 1] —— ... —— O[1, L] ——     |  | ——
    |  |    |                 |              |  |
    |X₀|    ⋮                  ⋮           =  | X| ⋮
    |  |    |                 |              |  |
    |  | —— O[W, 1] —— ... —— O[W, L] ——     |  | ——
    |  |    |                 |              |  |
    |  | -- B[1] ----- ... -- B[L] -----     |  | --
     ‾‾                                       ‾‾
"""
function pushleft(X₀::LeftEnvironmentTensor{N}, A::AbstractVector{<:LocalTensor{R}}, O::AbstractMatrix{MPOTensor}, B::AbstractVector{<:AdjointLocalTensor{R}}) where {N, R}
    N == size(O, 1) + 2 || throw(ArgumentError("Mismatched widths: $N ≠ $(size(O, 1)+2)"))
    (L = length(A)) == size(O, 2) == length(B) || throw(ArgumentError("Mismatched lengths: ($L, $(length(A)), $(size(O, 2)), $(length(B)))"))
    X = deepcopy(X₀)
    for l in 1:L
        X = pushleft(X, A[l], O[:, l], B[l])
    end
    return X
end
function pushleft(X₀::LeftEnvironmentTensor{N}, c::ChannelEnvironment{N, L, R}) where {N, L, R}
    return pushleft(X₀, c.A, c.O, c.B)
end
