"""
    function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::LocalTensor{R}, B::AdjointLocalTensor{R})
        -> X::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}

Contract left boundary tensor and transfer matrix
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
    function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::AbstractVector{<:LocalTensor{R}}, B::AbstractVector{<:AdjointLocalTensor{R}})
        -> X::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}

Iteratively contract left boundary tensor and transfer matrix
     -- A[1] -- ... -- A[L] --      --
    |   |              |           |
    X₀  |              |        =  X
    |   |              |           |
     -- B[1] -- ... -- B[L] --      --
"""
function pushleft(X₀::Union{LeftEnvironmentTensor{2}, BondTensor, AdjointBondTensor}, A::AbstractVector{<:LocalTensor{R}}, B::AbstractVector{<:AdjointLocalTensor{R}}) where R
    (L = length(A)) == length(B) || throw(ArgumentError("Mismatched lengths: $(length(A)) ≠ $(length(B))"))
    X = X₀
    for l in 1:L
        X = pushleft(X, A[l], B[l])
    end
    return X
end
