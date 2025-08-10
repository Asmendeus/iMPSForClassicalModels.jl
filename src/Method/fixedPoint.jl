"""
    function leftFixedPoint(t::TransferMatrix{L, R}, mode::Int64;
                alg::Union{SimpleIteration, Arnoldi, Lanczos}=SimpleIteration(),
                X₀::Vector{<:AbstractTensorWrapper{2}}=_default_X₀_leftFixedPoint(L, mode))

# Arguments
`t::TransferMatrix{L, R}`: a transfer matrix wrapper
`mode::Int64`: a flag finding which fixed point, ranging from {0, 1, 2}

# Keyword Arguments
`alg::Union{SimpleIteration, Arnoldi, Lanczos}=SimpleIteration()`: `SimpleIteration` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`
`X₀::Vector{<:AbstractTensorWrapper{2}}`: initial tensors, guessed solution of fixed point equations
    mode == 0 -> X₀::Vector{LeftEnvironmentTensor{2}}=[LeftEnvironmentTensor{2}(TensorMap(rand, ))]
    mode == 1 -> X₀::Vector{BondTensor}
    mode == 2 -> X₀::Vector{AdjointBondTensor}

# Return
`λ::Vector{<:Number}`: length `L` vector, coefficients of solution tensors of `L` fixed point equations
`X::Vector{<:AbstractTensorWrapper{2}}`: length `L` vector, solution tensors of `L` fixed point equations
    mode == 0 -> X₀::Vector{LeftEnvironmentTensor{2}}
    mode == 1 -> X₀::Vector{BondTensor}
    mode == 2 -> X₀::Vector{AdjointBondTensor}
(mode == 1) `AL::Vector{LocalTensor{R}}`: left orthogonal local tensors at fixed point
(mode == 2) `BL::Vector{AdjointLocalTensor{R}}`: left orthogonal adjoint local tensors at fixed point

# Fixed point equations
1. mode == 0

    for l = 1:L, the l-th fixed point equation is

         ---- A[l] --             --
        |     |                  |
        X[l]  |        =  λ[l] * X[l+1]
        |     |                  |
         ---- B[l] --             --

    where `X[L+1] = X[1]`, and `X` are left environment tensors normalized with normalization coefficients `λ`.
    The `L` fixed point equations can also be solved by the following largest eigenvalue equation

         ---- A[l] -- A[l+1] -- ... -- A[L] -- A[1] -- ... -- A[l-1] --                   --
        |     |       |                |       |              |                          |
        X[l]  |       |                |       |              |          =  (∏_l λ[l]) * X[l]
        |     |       |                |       |              |                          |
         ---- B[l] -- B[l+1] -- ... -- B[L] -- B[1] -- ... -- B[l-1] --                   --

    where `(∏_l λ[l])` is the dominant eigenvalues and `X` are the dominant eigenvectors.
    We can solve only for `l = 1`, and then the complete `λ` and `X` are generated from the fixed point equations above.

2. mode == 1

    for l = 1:L, the l-th fixed point equation is

        -- X[l] -- A[l] --  =  λ[l] * -- AL[l] -- X[l+1] --
                   |                     |        |

    where `X[L+1] = X[1]`, and `AL` are left-orthogonal local tensors and `X` are bond tensors normalized with normalization coefficients `λ`.
    Here we replace the original iterative equations with new equation iterative equations to speed up:

         ---- A[l] ----             --
        |     |                    |
        X[l]  |          =  λ[l] * X[l+1]
        |     |                    |
         ---- AL[l]* --             --

    For this case, adjoint local tensors `B` in the transfer matrix are used as the initial `AL*`.
    Such `L` fixed point equations correspond to `L` largest eigenvalue equations,

        X[l] --- A[l] --- A[l+1] --- ... -- A[L] --- A[1] --- ... -- A[l-1] ---                    X[l] ---
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |           =  (∏_l λ[l])^2 * |
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |                             |
        X[l]* -- A[l]* -- A[l+1]* -- ... -- A[L]* -- A[1]* -- ... -- A[l-1]* --                    X[l]* --

    where `(∏_l λ[l])^2` is the dominant eigenvalues and `X[l]†X[l]` is the dominant eigenvectors of the l-th eigenequation.
    We can solve only for `l = 1`, and use SVD decomposition `X[l]†X[l] = U S V†` to solve the `X[l] = √S V†`, where U = V and all singular values are positive.
    Then the complete `λ` and `X` are generated from the fixed point equations above.

3. mode == 2

    for l = 1:L, the l-th fixed point equation is

                   |                     |        |
        -- X[l] -- B[l] --  =  λ[l] * -- BL[l] -- X[l+1] --

    where `X[L+1] = X[1]`, and `BL` are left-orthogonal adjoint local tensors and `X` are adjoint bond tensors normalized with normalization coefficients `λ`.
    Here we replace the original iterative equations with new equation iterative equations to speed up:

         ---- BL[l]* --             --
        |     |                    |
        X[l]  |          =  λ[l] * X[l+1]
        |     |                    |
         ---- B[l] ----             --

    For this case, local tensors `A` in the transfer matrix are used as the initial `BL*`.
    Such `L` fixed point equations correspond to `L` largest eigenvalue equations,

        X[l]* -- B[l]* -- B[l+1]* -- ... -- B[L]* -- B[1]* -- ... -- B[l-1]* --                    X[l]* --
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |           =  (∏_l λ[l])^2 * |
        |        |        |                 |        |               |                             |
        |        |        |                 |        |               |                             |
        X[l] --- B[l] --- B[l+1] --- ... -- B[L] --- B[1] --- ... -- B[l-1] ---                    X[l] ---

    where `(∏_l λ[l])^2` is the dominant eigenvalues and `X[l]X[l]†` is the dominant eigenvectors of the l-th eigenequation.
    We can solve only for `l = 1`, and use SVD decomposition `X[l]X[l]† = U S V†` to solve the `X[l] = U √S`, where U = V and all singular values are positive.
    Then the complete `λ` and `X` are generated from the fixed point equations above.

"""
function leftFixedPoint(t::TransferMatrix{L, R}, mode::Int64;
            alg::Union{SimpleIteration, Arnoldi, Lanczos}=SimpleIteration(),
            X₀::Vector{<:AbstractTensorWrapper{2}}=_default_X₀_leftFixedPoint(L, mode))

    if mode == 0
    elseif mode == 1
    elseif mode == 2
    else
        throw(ArgumentError("Undefined behavior: mode = $mode ∉ {0, 1, 2}"))
    end
end

"""

"""
function rightFixedPoint(t::TransferMatrix{L, R}, mode::Int64;
            alg::Union{SimpleIteration, Arnoldi, Lanczos}=SimpleIteration(),
            X₀::Vector{<:AbstractTensorWrapper{2}}=_default_X₀_rightFixedPoint(L, mode))
    if mode == 0
    elseif mode == 1
    elseif mode == 2
    else
        throw(ArgumentError("Undefined behavior: mode = $mode ∉ {0, 1, 2}"))
    end
end

# ========== Auxiliary functions generating default X₀ ==========
function _default_X₀_leftFixedPoint(L::Int, mode::Int64)
end
function _default_X₀_rightFixedPoint(L::Int, mode::Int64)
end
