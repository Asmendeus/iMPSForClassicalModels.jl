"""
    function rightFixedPoint(env::S, X₀::T, alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()) where {S, T}

A series of functions for solving right/lower fixed point equations or maximum eigenequations.

# Arguments
`env::S`: environment for solving fixed point equations or maximum eigenequations
`X₀::T`: initial state
`alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()`: `SimpleIteration` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# ===============================================
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
              |         |                |       |              |      X[l]  =  (∏_l λ[l]) *   X[l-1]
              |         |                |       |              |      |                       |
         ---- B[l+1] -- B[l+2] -- ... -- B[L] -- B[1] -- ... -- B[l] --                      --

    where `(∏_l λ[l])` is the dominant eigenvalues and `X` are the dominant eigenvectors.
    We solve only for `l = L`, and then the complete `λ` and `X` are generated from the fixed point equations above.
"""
function rightFixedPoint(t::TransferMatrix{L, R},
            X₀::Vector{<:RightEnvironmentTensor{2}}=_default_X₀_rightFixedPoint(t),
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=SimpleIteration(; tol=repeat([Defaults.tol,], L))) where {L, R}

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
function _default_X₀_rightFixedPoint(t::TransferMatrix{L, R}) where {L, R}
    tensortype = RightEnvironmentTensor{2}
    datatype = domain(t.A[end], R-2) isa CartesianSpace ? Float64 : ComplexF64
    X₀ = tensortype[TensorMap(rand, datatype, domain(t.A[l], R-2), codomain(t.B[l], R-2)) for l in 1:L]
    return X₀
end
