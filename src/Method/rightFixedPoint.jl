"""
    function rightFixedPoint(env::S, X₀::T, alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()) where {S, T}

A series of functions for solving right/lower fixed point equations or maximum eigenequations.

# Arguments
`env::S`: environment for solving fixed point equations or maximum eigenequations
`X₀::T`: initial state
`alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()`: `SimpleIteration` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# ===============================================
    function rightFixedPoint(A::Vector{<:LocalTensor{R}},
                X₀::Vector{<:BondTensor}=_default_X₀_rigthFixedPoint(A),
                alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi();
                kwargs...) where R

# Arguments
`A::Vector{<:LocalTensor{R}}`: vector of LocalTensor{R}
`X₀::Vector{<:BondTensor}=_default_X₀_rightFixedPoint(A)`: initial tensors, guessed solution of fixed point equations
`alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()`: `SimpleIteration` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# Keyword Arguments
`kwargs`: keyword arguments for `rightorth`

# Return
`λ::Vector{<:Number}`: length `L` vector, coefficients of solution tensors of `L` fixed point equations
`X::Vector{RightEnvironmentTensor{2}}`: length `L` vector, solution tensors of `L` fixed point equations
`AR::Vector{LocalTensor{R}}`: rigtht orthogonal local tensors at fixed point

# Fixed point equations

    for l = 1:L, the l-th fixed point equation is

        -- A[l] -- X[l] --  =  λ[l] * -- X[l-1] -- AR[l] --
           |                                       |

    where `X[0] = X[L]`, and `AR` are right-orthogonal local tensors and `X` are bond tensors normalized with normalization coefficients `λ`.
    Here we replace the original iterative equations with new equation iterative equations to speed up:

        ---- A[l] ----                --
             |        |                 |
             |        X[l]  =  λ[l] *   X[l-1]
             |        |                 |
        ---- AR[l]* --                --

# Maximum eigenequations:

    for l in 1:L, the l-th maximum eigenequation is

        -- A[l+1] --- A[l+2] --- ... -- A[L] --- A[1] --- ... -- A[l] --- X[l]                     -- X[l] ---
           |          |                 |        |               |        |                                   |
           |          |                 |        |               |        |                                   |
           |          |                 |        |               |        |      =  (∏_l λ[l])^2 *            |
           |          |                 |        |               |        |                                   |
           |          |                 |        |               |        |                                   |
        -- A[l+1]* -- A[l+2]* -- ... -- A[L]* -- A[1]* -- ... -- A[l]* -- X[l]*                    -- X[l]* --

    where `(∏_l λ[l])^2` is the dominant eigenvalues and `X[l]X[l]†` is the dominant eigenvectors of the l-th eigenequation.
    We solve only for `l = L`, and use SVD decomposition `X[l]X[l]† = U S V†` to solve the `X[l] = U √S`, where U = V and all singular values are positive.
    Then the complete `λ` and `X` are generated from the fixed point equations above.
"""
function rightFixedPoint(A::Vector{<:LocalTensor{R}},
            X₀::Vector{<:BondTensor}=_default_X₀_rightFixedPoint(A),
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi();
            kwargs...) where R
    (L = length(A)) == length(X₀) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(X₀))"))
    if alg isa SimpleIteration
        func = [x -> pushright(x, A[l], rightorth((A[l] * x)')[2]; kwargs...) for l in 1:L]
        λ, X, info = iterate(func, X₀, false, alg)
        AR = [rightorth(A[l] * X[l])[2] for l in 1:L]
        return λ, X, AR, info
    elseif alg isa KrylovKit.KrylovAlgorithm
        func = x -> pushright(x, A, adjoint.(A))
        vals, vecs, info = eigsolve(func, X₀[end], 1, :LM, alg)

        λ = Vector{typeof(vals[1])}(undef, L)
        X = Vector{BondTensor}(undef, L)
        AR = Vector{LocalTensor{R}}(undef, L)

        u, s, _, _ = tsvd(vecs[1].A)

        X[end] = BondTensor(u * sqrt(s))
        X[end] /= sign_first_element(X[end])
        # Notice AR[end] and λ[end] need to be re-updated
        for l in L:-1:0
            l_now = mod(l-1, L) + 1
            l₋ = mod(l-2, L) + 1
            x, AR[l_now] = rightorth(A[l_now] * X[l_now]; kwargs...)
            λ[l_now] = norm(x) * sign_first_element(x)
            X[l₋] = x / λ[l_now]
        end

        if alg isa Lanczos
            info = convert(LanczosInfo, info)
        elseif alg isa Arnoldi
            info = convert(ArnoldiInfo, info)
        else
            @warn "$(typeof(alg)) is available but its information type is not defined, please contact the author."
        end

        return λ, X, AR, info
    end
end
function rightFixedPoint(A::Vector{<:LocalTensor{R}}, alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}) where R
    return rightFixedPoint(A, _default_X₀_rightFixedPoint(A), alg)
end

"""
    function rightFixedPoint(B::Vector{<:AdjointLocalTensor{R}},
                X₀::Vector{<:AdjointBondTensor}=_default_X₀_rightFixedPoint(B),
                alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()) where R

# Arguments
`B::Vector{<:AdjointLocalTensor{R}}`: vector of AdjointLocalTensor{R}
`X₀::Vector{AdjointBondTensor}=_default_X₀_rightFixedPoint(B)`: initial tensors, guessed solution of fixed point equations
`alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()`: `SimpleIteration` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# Keyword Arguments
`kwargs`: keyword arguments for `rightorth`

# Return
`λ::Vector{<:Number}`: length `L` vector, coefficients of solution tensors of `L` fixed point equations
`X::Vector{RightEnvironmentTensor{2}}`: length `L` vector, solution tensors of `L` fixed point equations
`BR::Vector{AdjointLocalTensor{R}}`: right orthogonal adjoint local tensors at fixed point

# Fixed point equations

    for l = 1:L, the l-th fixed point equation is

           |                                       |
        -- B[l] -- X[l] --  =  λ[l] * -- X[l-1] -- BR[l] --

    where `X[0] = X[L]`, and `BR` are right-orthogonal adjoint local tensors and `X` are adjoint bond tensors normalized with normalization coefficients `λ`.
    Here we replace the original iterative equations with new equation iterative equations to speed up:

        ---- BR[l]* --                --
             |        |                 |
             |        X[l]  =  λ[l] *   X[l-1]
             |        |                 |
        ---- B[l] ----                --

# Maximum eigenequations

    for l = 1:L, the l-th eigenequation is

        -- B[l+1]* -- B[l+2]* -- ... -- B[L]* -- B[1]* -- ... -- B[l]* -- X[l]*                    -- X[l]*
           |          |                 |        |               |         |                          |
           |          |                 |        |               |         |                          |
           |          |                 |        |               |         |     =  (∏_l λ[l])^2 *    |
           |          |                 |        |               |         |                          |
           |          |                 |        |               |         |                          |
        -- B[l+1] --- B[l+2] --- ... -- B[L] --- B[1] --- ... -- B[l] --- X[l]                     -- X[l]

    where `(∏_l λ[l])^2` is the dominant eigenvalues and `X[l]†X[l]` is the dominant eigenvectors of the l-th eigenequation.
    We solve only for `l = L`, and use SVD decomposition `X[l]†X[l] = U S V†` to solve the `X[l] = √S V†`, where U = V and all singular values are positive.
    Then the complete `λ` and `X` are generated from the fixed point equations above.
"""
function rightFixedPoint(B::Vector{<:AdjointLocalTensor{R}},
            X₀::Vector{<:AdjointBondTensor}=_default_X₀_rightFixedPoint(B),
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi();
            kwargs...) where R
    (L = length(B)) == length(X₀) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(X₀))"))
    if alg isa SimpleIteration
        func = [x -> pushright(x, rightorth((B[l] * x)')[2], B[l]; kwargs...) for l in 1:L]
        λ, X, info = iterate(func, X₀, false, alg)
        BR = [rightorth(B[l] * X[l])[2] for l in 1:L]
        return λ, X, BR, info
    elseif alg isa KrylovKit.KrylovAlgorithm
        func = x -> pushright(x, adjoint.(B), B)
        vals, vecs, info = eigsolve(func, X₀[end], 1, :LM, alg)

        λ = Vector{typeof(vals[1])}(undef, L)
        X = Vector{AdjointBondTensor}(undef, L)
        BR = Vector{AdjointLocalTensor{R}}(undef, L)

        _, s, vd, _ = tsvd(vecs[1].A)

        X[end] = AdjointBondTensor(sqrt(s) * vd)
        X[end] /= sign_first_element(X[end])
        # Notice BR[end] and λ[end] need to be re-updated
        for l in L:-1:0
            l_now = mod(l-1, L) + 1
            l₋ = mod(l-2, L) + 1
            x, BR[l_now] = rightorth(B[l_now] * X[l_now]; kwargs...)
            λ[l_now] = norm(x) * sign_first_element(x)
            X[l₋] = x / λ[l_now]
        end

        if alg isa Lanczos
            info = convert(LanczosInfo, info)
        elseif alg isa Arnoldi
            info = convert(ArnoldiInfo, info)
        else
            @warn "$(typeof(alg)) is available but its information type is not defined, please contact the author."
        end

        return λ, X, BR, info
    end
end
function rightFixedPoint(B::Vector{<:AdjointLocalTensor{R}}, alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}) where R
    return rightFixedPoint(B, _default_X₀_rightFixedPoint(B), alg)
end

"""
    function rightFixedPoint(t::TransferMatrix{L, R},
                X₀::Vector{<:RightEnvironmentTensor{2}}=_default_X₀_rightFixedPoint(t),
                alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=SimpleIteration(; tol=repeat([Defaults.tol,], L))) where {L, R}

# Arguments
`t::TransferMatrix{L, R}`: a transfer matrix wrapper
`X₀::Vector{RightEnvironmentTensor{2}}=_default_X₀_rightFixedPoint(t)`: initial tensors, guessed solution of fixed point equations
`alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=SimpleIteration()`: `SimpleIteration` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# Return
`λ::Vector{<:Number}`: length `L` vector, coefficients of solution tensors of `L` fixed point equations
`X::Vector{RightEnvironmentTensor{2}}`: length `L` vector, solution tensors of `L` fixed point equations

# Fixed point equations

    for l = 1:L, the l-th fixed point equation is

        ---- A[l] --                --
             |      |                 |
             |      X[l]  =  λ[l] *   X[l-1]
             |      |                 |
        ---- B[l] --                --

    where `X[0] = X[L]`, and `X` are right environment tensors normalized with normalization coefficients `λ`.

# Maximum eigenequations

    for l in 1:L, the l-th maximum eigenequation is

        ---- A[l+1] -- A[l+2] -- ... -- A[L] -- A[1] -- ... -- A[l] --                      --
             |         |                |       |              |      |                       |
             |         |                |       |              |      X[l]  =  (∏_l λ[l]) *   X[l]
             |         |                |       |              |      |                       |
        ---- B[l+1] -- B[l+2] -- ... -- B[L] -- B[1] -- ... -- B[l] --                      --

    where `(∏_l λ[l])` is the dominant eigenvalues and `X` are the dominant eigenvectors.
    We solve only for `l = L`, and then the complete `λ` and `X` are generated from the fixed point equations above.
"""
function rightFixedPoint(t::TransferMatrix{L, R},
            X₀::Vector{<:RightEnvironmentTensor{2}}=_default_X₀_rightFixedPoint(t),
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=SimpleIteration(; tol=repeat([Defaults.tol,], L))) where {L, R}
    L == length(X₀) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(X₀))"))
    if alg isa SimpleIteration
        func = [x -> pushright(x, t.A[l], t.B[l]) for l in 1:L]
        λ, X, info = iterate(func, X₀, false, alg)
        return λ, X, info
    elseif alg isa KrylovKit.KrylovAlgorithm
        func = x -> pushright(x, t.A, t.B)
        vals, vecs, info = eigsolve(func, X₀[end], 1, :LM, alg)

        λ = Vector{typeof(vals[1])}(undef, L)
        X = Vector{RightEnvironmentTensor{2}}(undef, L)

        X[end] = vecs[1]
        X[end] /= sign_first_element(X[end])
        for l in L:-1:1
            l₋ = mod(l-2, L) + 1
            vec = pushright(X[l], t.A[l], t.B[l])
            λ[l] = norm(vec) * sign_first_element(vec)
            X[l₋] = vec / λ[l]
        end

        if alg isa Lanczos
            info = convert(LanczosInfo, info)
        elseif alg isa Arnoldi
            info = convert(ArnoldiInfo, info)
        else
            @warn "$(typeof(alg)) is available but its information type is not defined, please contact the author."
        end

        return λ, X, info
    end
end
function rightFixedPoint(t::TransferMatrix{L, R}, alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}) where {L, R}
    return rightFixedPoint(t, _default_X₀_rightFixedPoint(t), alg)
end


# ========== Auxiliary functions generating default X₀ ==========
function _default_X₀_rightFixedPoint(A::Vector{<:LocalTensor{R}}) where R
    L = length(A)
    tensortype = BondTensor
    datatype = domain(A[1], R-2) isa CartesianSpace ? Float64 : ComplexF64
    X₀ = tensortype[TensorMap(rand, datatype, domain(A[l], R-2), domain(A[l], R-2)) for l in 1:L]
    return X₀
end
function _default_X₀_rightFixedPoint(B::Vector{<:AdjointLocalTensor{R}}) where R
    L = length(B)
    tensortype = AdjointBondTensor
    datatype = codomain(B[end], R-2) isa CartesianSpace ? Float64 : ComplexF64
    X₀ = tensortype[TensorMap(rand, datatype, codomain(B[l], R-2), codomain(B[l], R-2)) for l in 1:L]
    return X₀
end
function _default_X₀_rightFixedPoint(t::TransferMatrix{L, R}) where {L, R}
    tensortype = RightEnvironmentTensor{2}
    datatype = domain(t.A[end], R-2) isa CartesianSpace ? Float64 : ComplexF64
    X₀ = tensortype[TensorMap(rand, datatype, domain(t.A[l], R-2), codomain(t.B[l], R-2)) for l in 1:L]
    return X₀
end
