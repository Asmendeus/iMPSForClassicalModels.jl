"""
    midFixedPoint(env::S, X₀::T, alg::EigenAlgorithm=Defaults.alg_eig) where {S, T}

A series of functions for solving mid fixed point equations or maximum eigenequations.

# Arguments
`env::S`: environment for solving fixed point equations or maximum eigenequations
`X₀::T`: initial state
`alg::EigenAlgorithm=Defaults.alg_eig`: `SimpleIterator` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# Return in common
`λ::Vector{<:Number}`: length `L` vector, coefficients of solution tensors of `L` fixed point equations
`X::Vector{<:AbstractTensorWrapper}`: length `L` vector, solution tensors of `L` fixed point equations
`info`: information of the algorithm

# Notes
By convention, we appoint that:
the eigenvector `X` must meet `norm(X) = 1` and have a positive real first element

# ===============================================
    midFixedPoint(env::BondEnvironment{N},
                X₀::Union{BondTensor, AdjointBondTensor}=_default_X₀_midFixedPoint(env, false),
                alg::EigenAlgorithm=Defaults.alg_eig;
                which::Symbol=:LM) where N

# Arguments
`env::BondEnvironment{N}`: bond environment
`X₀::Union{BondTensor, AdjointBondTensor}=_default_X₀_midFixedPoint(env, false)`: initial tensor, guessed solution of fixed point equation
`alg::EigenAlgorithm=Defaults.alg_eig`: `SimpleIterator` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# Keyword Arguments
`which::Symbol=:LM`: Which eigenvalue to solve by Krylov algorithm

# Return
`λ::Number`: coefficient of solution tensor of fixed point equation
`X::Union{BondTensor, AdjointBondTensor}`: solution tensor of fixed point equation
`info`: information of the algorithm

# Fixed point equation & Maximum eigenequation
    if X₀ isa BondTensor
         __           __
        |  | -- X -- |  |
        |  |         |  |
        |  | ——————— |  |
        |  |         |  |
        |FL|    ⋮    |FR|  =  λ * -- X --
        |  |         |  |
        |  | ——————— |  |
        |  |         |  |
        |  | --   -- |  |
         ‾‾           ‾‾
    else X₀ isa AdjointBondTensor
         __           __
        |  | --   -- |  |
        |  |         |  |
        |  | ——————— |  |
        |  |         |  |
        |FL|    ⋮    |FR|  =  λ * -- X --
        |  |         |  |
        |  | ——————— |  |
        |  |         |  |
        |  | -- X -- |  |
         ‾‾           ‾‾
"""
function midFixedPoint(env::BondEnvironment{N},
            X₀::Union{BondTensor, AdjointBondTensor}=_default_X₀_midFixedPoint(env, false),
            alg::EigenAlgorithm=Defaults.alg_eig;
            which::Symbol=:LM) where N
    if alg isa SimpleIterator
        func = x -> pushmid(x, env.FL, env.FR)
        λ, X, info = iterate(func, X₀, alg)
        return λ, X, info
    elseif alg isa KrylovKit.KrylovAlgorithm
        func = x -> pushmid(x, env.FL, env.FR)
        vals, vecs, info = eigsolve(func, X₀, 1, which, alg)

        λ = vals[1]
        X = vecs[1] / lambda(vecs[1])

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
function midFixedPoint(env::BondEnvironment{N}, alg::EigenAlgorithm; which::Symbol=:LM) where N
    return midFixedPoint(env, _default_X₀_midFixedPoint(env, false), alg; which=which)
end

"""
    midFixedPoint(env::CenterEnvironment{N},
                X₀::Union{LocalTensor{R}, AdjointLocalTensor{R}}=_default_X₀_midFixedPoint(env, false),
                alg::EigenAlgorithm=Defaults.alg_eig;
                which::Symbol=:LM) where {N, R}

# Arguments
`env::CenterEnvironment{N}`: center environment
`X₀::Union{LocalTensor{R}, AdjointLocalTensor{R}}=_default_X₀_midFixedPoint(env, false)`: initial tensor, guessed solution of fixed point equation
`alg::EigenAlgorithm=Defaults.alg_eig`: `SimpleIterator` for `iterate`, while `Arnoldi` or `Lanczos` for `eigsolve`

# Keyword Arguments
`which::Symbol=:LM`: Which eigenvalue to solve by Krylov algorithm

# Return
`λ::Number`: coefficient of solution tensor of fixed point equation
`X::Union{LocalTensor{R}, AdjointLocalTensor{R}}`: solution tensor of fixed point equation
`info`: information of the algorithm

# Fixed point equation & Maximum eigenequation
    if X₀ isa LocalTensor{R}
               (a)
         __     |        __
        |  | -- X ----- |  |
        |  |    |       |  |
        |  | —— O[1] —— |  |           (a)
        |  |            |  |            |
        |FL|    ⋮       |FR|  =  λ * -- X -----
        |  |    |       |  |            |
        |  | —— O[W] —— |  |
        |  |    |       |  |
        |  | --   ----- |  |
         ‾‾              ‾‾
    else X₀ isa AdjointLocalTensor{R}
         __              __
        |  | --   ----- |  |
        |  |    |       |  |
        |  | —— O[1] —— |  |
        |  |            |  |            |
        |FL|    ⋮       |FR|  =  λ * -- X -----
        |  |    |       |  |            |
        |  | —— O[W] —— |  |           (a)
        |  |    |       |  |
        |  | -- X ----- |  |
         ‾‾     |        ‾‾
               (a)
"""
function midFixedPoint(env::CenterEnvironment{N},
            X₀::Union{LocalTensor{R}, AdjointLocalTensor{R}}=_default_X₀_midFixedPoint(env, false),
            alg::EigenAlgorithm=Defaults.alg_eig;
            which::Symbol=:LM) where {N, R}
    if alg isa SimpleIterator
        func = x -> pushmid(x, env.FL, env.O, env.FR)
        λ, X, info = iterate(func, X₀, alg)
        return λ, X, info
    elseif alg isa KrylovKit.KrylovAlgorithm
        func = x -> pushmid(x, env.FL, env.O, env.FR)
        vals, vecs, info = eigsolve(func, X₀, 1, which, alg)

        λ = vals[1]
        X = vecs[1] / lambda(vecs[1])

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
function midFixedPoint(env::CenterEnvironment{N}, alg::EigenAlgorithm; which::Symbol=:LM) where N
    return midFixedPoint(env, _default_X₀_midFixedPoint(env, false), alg; which=which)
end


# ========== Auxiliary functions generating default X₀ ==========
function _default_X₀_midFixedPoint(env::BondEnvironment{N}, isadjoint::Bool) where N
    datatype = codomain(env.FL, 1) isa CartesianSpace ? Float64 : ComplexF64
    if isadjoint
        return AdjointBondTensor(TensorMap(rand, datatype, domain(env.FR, 1), codomain(env.FL, 1)))
    else
        return BondTensor(TensorMap(rand, datatype, domain(env.FL, N-1), codomain(env.FR, 1)))
    end
end
function _default_X₀_midFixedPoint(env::CenterEnvironment{N}, isadjoint::Bool) where N
    datatype = codomain(env.FL, 1) isa CartesianSpace ? Float64 : ComplexF64
    if isadjoint
        return AdjointMPSTensor(TensorMap(rand, datatype, domain(env.FR, 1), codomain(env.FL, 1) ⊗ codomain(env.O[end], 2)))
    else
        return MPSTensor(TensorMap(rand, datatype, domain(env.FL, N-1) ⊗ domain(env.O[1], 1), codomain(env.FR, 1)))
    end
end
