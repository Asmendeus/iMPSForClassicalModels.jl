"""
    function leftFixedPoint(env::S, X₀::T, alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()) where {S, T}

A series of functions for solving left/upper fixed point equations or maximum eigenequations.

# Arguments
`env::S`: environment for solving fixed point equations or maximum eigenequations
`X₀::T`: initial state
`alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()`: `SimpleIteration` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# ===============================================
    function leftFixedPoint(A::Vector{<:LocalTensor{R}},
                X₀::Vector{<:BondTensor}=_default_X₀_leftFixedPoint(A),
                alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi();
                kwargs...) where R

# Arguments
`A::Vector{<:LocalTensor{R}}`: vector of LocalTensor{R}
`X₀::Vector{<:BondTensor}=_default_X₀_leftFixedPoint(A)`: initial tensors, guessed solution of fixed point equations
`alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()`: `SimpleIteration` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# Keyword Arguments
`kwargs`: keyword arguments for `leftorth`

# Return
`λ::Vector{<:Number}`: length `L` vector, coefficients of solution tensors of `L` fixed point equations
`X::Vector{LeftEnvironmentTensor{2}}`: length `L` vector, solution tensors of `L` fixed point equations
`BL::Vector{LocalTensor{R}}`: left orthogonal local tensors at fixed point

# Fixed point equations

    for l = 1:L, the l-th fixed point equation is

        -- X[l] -- A[l] --  =  λ[l] * -- BL[l] -- X[l+1] --
                   |                     |

    where `X[L+1] = X[1]`, and `BL` are left-orthogonal local tensors and `X` are bond tensors normalized with normalization coefficients `λ`.
    Here we replace the original iterative equations with new equation iterative equations to speed up:

         ---- A[l] ----             --
        |     |                    |
        X[l]  |          =  λ[l] * X[l+1]
        |     |                    |
         ---- BL[l]* --             --

# Maximum eigenequations:

    for l in 1:L, the l-th maximum eigenequation is

        X[l] --- A[l] --- A[l+1] --- ... -- A[L] --- A[1] --- ... -- A[l-1] ---                    X[l] ---
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |           =  (∏_l λ[l])^2 * |
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |                             |
        X[l]* -- A[l]* -- A[l+1]* -- ... -- A[L]* -- A[1]* -- ... -- A[l-1]* --                    X[l]* --

    where `(∏_l λ[l])^2` is the dominant eigenvalues and `X[l]†X[l]` is the dominant eigenvectors of the l-th eigenequation.
    We solve only for `l = 1`, and use SVD decomposition `X[l]†X[l] = U S V†` to solve the `X[l] = √S V†`, where U = V and all singular values are positive.
    Then the complete `λ` and `X` are generated from the fixed point equations above.
"""
function leftFixedPoint(A::Vector{<:LocalTensor{R}},
            X₀::Vector{<:BondTensor}=_default_X₀_leftFixedPoint(A),
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi();
            kwargs...) where R
    (L = length(A)) == length(X₀) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(X₀))"))
    if alg isa SimpleIteration
        func = [x -> pushleft(x, A[l], leftorth((x * A[l])')[1]; kwargs...) for l in 1:L]
        λ, X, info = iterate(func, X₀, true, alg)
        BL = [leftorth(X[l] * A[l])[1] for l in 1:L]
        return λ, X, BL, info
    elseif alg isa KrylovKit.KrylovAlgorithm
        func = x -> pushleft(x, A, adjoint.(A))
        vals, vecs, info = eigsolve(func, X₀[1], 1, :LM, alg)

        λ = Vector{typeof(vals[1])}(undef, L)
        X = Vector{BondTensor}(undef, L)
        BL = Vector{LocalTensor{R}}(undef, L)

        _, s, vd, _ = tsvd(vecs[1].A)

        X[1] = BondTensor(sqrt(s) * vd)
        X[1] /= sign_first_element(X[1])
        # Notice BL[1] and λ[1] need to be re-updated
        for l in 1:L+1
            l_now = mod(l-1, L) + 1
            l₊ = mod(l, L) + 1
            BL[l_now], x = leftorth(X[l_now] * A[l_now]; kwargs...)
            λ[l_now] = norm(x) * sign_first_element(x)
            X[l₊] = x / λ[l_now]
        end

        if alg isa Lanczos
            info = convert(LanczosInfo, info)
        elseif alg isa Arnoldi
            info = convert(ArnoldiInfo, info)
        else
            @warn "$(typeof(alg)) is available but its information type is not defined, please contact the author."
        end

        return λ, X, BL, info
    end
end
function leftFixedPoint(A::Vector{<:LocalTensor{R}}, alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}) where R
    return leftFixedPoint(A, _default_X₀_leftFixedPoint(A), alg)
end

"""
    function leftFixedPoint(B::Vector{<:AdjointLocalTensor{R}},
                X₀::Vector{<:AdjointBondTensor}=_default_X₀_leftFixedPoint(B),
                alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()) where R

# Arguments
`B::Vector{<:AdjointLocalTensor{R}}`: vector of AdjointLocalTensor{R}
`X₀::Vector{AdjointBondTensor}=_default_X₀_leftFixedPoint(B)`: initial tensors, guessed solution of fixed point equations
`alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()`: `SimpleIteration` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# Keyword Arguments
`kwargs`: keyword arguments for `leftorth`

# Return
`λ::Vector{<:Number}`: length `L` vector, coefficients of solution tensors of `L` fixed point equations
`X::Vector{LeftEnvironmentTensor{2}}`: length `L` vector, solution tensors of `L` fixed point equations
`BL::Vector{AdjointLocalTensor{R}}`: left orthogonal adjoint local tensors at fixed point

# Fixed point equations

    for l = 1:L, the l-th fixed point equation is

                   |                     |
        -- X[l] -- B[l] --  =  λ[l] * -- BL[l] -- X[l+1] --

    where `X[L+1] = X[1]`, and `BL` are left-orthogonal adjoint local tensors and `X` are adjoint bond tensors normalized with normalization coefficients `λ`.
    Here we replace the original iterative equations with new equation iterative equations to speed up:

         ---- BL[l]* --             --
        |     |                    |
        X[l]  |          =  λ[l] * X[l+1]
        |     |                    |
         ---- B[l] ----             --

# Maximum eigenequations

    for l = 1:L, the l-th eigenequation is

        X[l]* -- B[l]* -- B[l+1]* -- ... -- B[L]* -- B[1]* -- ... -- B[l-1]* --                    X[l]* --
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |           =  (∏_l λ[l])^2 * |
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |                             |
        X[l] --- B[l] --- B[l+1] --- ... -- B[L] --- B[1] --- ... -- B[l-1] ---                    X[l] ---

    where `(∏_l λ[l])^2` is the dominant eigenvalues and `X[l]X[l]†` is the dominant eigenvectors of the l-th eigenequation.
    We solve only for `l = 1`, and use SVD decomposition `X[l]X[l]† = U S V†` to solve the `X[l] = U √S`, where U = V and all singular values are positive.
    Then the complete `λ` and `X` are generated from the fixed point equations above.
"""
function leftFixedPoint(B::Vector{<:AdjointLocalTensor{R}},
            X₀::Vector{<:AdjointBondTensor}=_default_X₀_leftFixedPoint(B),
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi();
            kwargs...) where R
    (L = length(B)) == length(X₀) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(X₀))"))
    if alg isa SimpleIteration
        func = [x -> pushleft(x, leftorth((x * B[l])')[1], B[l]; kwargs...) for l in 1:L]
        λ, X, info = iterate(func, X₀, true, alg)
        BL = [leftorth(X[l] * B[l])[1] for l in 1:L]
        return λ, X, BL, info
    elseif alg isa KrylovKit.KrylovAlgorithm
        func = x -> pushleft(x, adjoint.(B), B)
        vals, vecs, info = eigsolve(func, X₀[1], 1, :LM, alg)

        λ = Vector{typeof(vals[1])}(undef, L)
        X = Vector{AdjointBondTensor}(undef, L)
        BL = Vector{AdjointLocalTensor{R}}(undef, L)

        u, s, _, _ = tsvd(vecs[1].A)

        X[1] = AdjointBondTensor(u * sqrt(s))
        X[1] /= sign_first_element(X[1])
        # Notice BL[1] and λ[1] need to be re-updated
        for l in 1:L+1
            l_now = mod(l-1, L) + 1
            l₊ = mod(l, L) + 1
            BL[l_now], x = leftorth(X[l_now] * B[l_now]; kwargs...)
            λ[l_now] = norm(x) * sign_first_element(x)
            X[l₊] = x / λ[l_now]
        end

        if alg isa Lanczos
            info = convert(LanczosInfo, info)
        elseif alg isa Arnoldi
            info = convert(ArnoldiInfo, info)
        else
            @warn "$(typeof(alg)) is available but its information type is not defined, please contact the author."
        end

        return λ, X, BL, info
    end
end
function leftFixedPoint(B::Vector{<:AdjointLocalTensor{R}}, alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}) where R
    return leftFixedPoint(B, _default_X₀_leftFixedPoint(B), alg)
end

"""

    function leftFixedPoint(t::TransferMatrix{L, R},
                X₀::Vector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(t),
                alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()) where {L, R}

# Arguments
`t::TransferMatrix{L, R}`: a transfer matrix wrapper
`X₀::Vector{LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(t)`: initial tensors, guessed solution of fixed point equations
`alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()`: `SimpleIteration` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# Return
`λ::Vector{<:Number}`: length `L` vector, coefficients of solution tensors of `L` fixed point equations
`X::Vector{LeftEnvironmentTensor{2}}`: length `L` vector, solution tensors of `L` fixed point equations

# Fixed point equations

    for l = 1:L, the l-th fixed point equation is

         ---- A[l] --             --
        |     |                  |
        X[l]  |        =  λ[l] * X[l+1]
        |     |                  |
         ---- B[l] --             --

    where `X[L+1] = X[1]`, and `X` are left environment tensors normalized with normalization coefficients `λ`.

# Maximum eigenequations

    for l in 1:L, the l-th maximum eigenequation is

         ---- A[l] -- A[l+1] -- ... -- A[L] -- A[1] -- ... -- A[l-1] --                   --
        |     |       |                |       |              |                          |
        X[l]  |       |                |       |              |          =  (∏_l λ[l]) * X[l]
        |     |       |                |       |              |                          |
         ---- B[l] -- B[l+1] -- ... -- B[L] -- B[1] -- ... -- B[l-1] --                   --

    where `(∏_l λ[l])` is the dominant eigenvalues and `X` are the dominant eigenvectors.
    We solve only for `l = 1`, and then the complete `λ` and `X` are generated from the fixed point equations above.
"""
function leftFixedPoint(t::TransferMatrix{L, R},
            X₀::Vector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(t),
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()) where {L, R}
    L == length(X₀) || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(X₀))"))
    if alg isa SimpleIteration
        func = [x -> pushleft(x, t.A[l], t.B[l]) for l in 1:L]
        λ, X, info = iterate(func, X₀, true, alg)
        return λ, X, info
    elseif alg isa KrylovKit.KrylovAlgorithm
        func = x -> pushleft(x, t.A, t.B)
        vals, vecs, info = eigsolve(func, X₀[1], 1, :LM, alg)

        λ = Vector{typeof(vals[1])}(undef, L)
        X = Vector{LeftEnvironmentTensor{2}}(undef, L)

        X[1] = vecs[1]
        X[1] /= sign_first_element(X[1])
        for l in 1:L
            l₊ = mod(l, L) + 1
            vec = pushleft(X[l], t.A[l], t.B[l])
            λ[l] = norm(vec) * sign_first_element(vec)
            X[l₊] = vec / λ[l]
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
function leftFixedPoint(t::TransferMatrix{L, R}, alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}) where {L, R}
    return leftFixedPoint(t, _default_X₀_leftFixedPoint(t), alg)
end


# ========== Auxiliary functions generating default X₀ ==========
function _default_X₀_leftFixedPoint(A::Vector{<:LocalTensor{R}}) where R
    L = length(A)
    tensortype = BondTensor
    datatype = codomain(A[1], 1) isa CartesianSpace ? Float64 : ComplexF64
    X₀ = tensortype[TensorMap(rand, datatype, codomain(A[l], 1), codomain(A[l], 1)) for l in 1:L]
    return X₀
end
function _default_X₀_leftFixedPoint(B::Vector{<:AdjointLocalTensor{R}}) where R
    L = length(B)
    tensortype = AdjointBondTensor
    datatype = domain(B[1], 1) isa CartesianSpace ? Float64 : ComplexF64
    X₀ = tensortype[TensorMap(rand, datatype, domain(B[l], 1), domain(B[l], 1)) for l in 1:L]
    return X₀
end
function _default_X₀_leftFixedPoint(t::TransferMatrix{L, R}) where {L, R}
    tensortype = LeftEnvironmentTensor{2}
    datatype = codomain(t.A[1], 1) isa CartesianSpace ? Float64 : ComplexF64
    X₀ = tensortype[TensorMap(rand, datatype, domain(t.B[l], 1), codomain(t.A[l], 1)) for l in 1:L]
    return X₀
end