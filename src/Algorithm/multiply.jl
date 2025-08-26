"""
    multiply(ψ₀::DenseInfiniteMPS{L},
        ρ::InfiniteMPO{L},
        ψ::DenseInfiniteMPS{L},
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig;
        FL₀::AbstractVector{<:LeftEnvironmentTensor{3}}=_default_X₀_leftFixedPoint(environment(getA(ψ), getA(ρ), adjoint.(getA(ψ₀)))),
        FR₀::AbstractVector{<:RightEnvironmentTensor{3}}=_default_X₀_rightFixedPoint(environment(getA(ψ), getA(ρ), adjoint.(getA(ψ₀))))) -> ψ′::DenseInfiniteMPS{L}
    multiply(ψ₀::DenseInfiniteMPS{L},
        ρ::SparseUMPO{W, L},
        ψ::DenseInfiniteMPS{L},
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig;
        FL₀::AbstractVector{<:LeftEnvironmentTensor{3}}=_default_X₀_leftFixedPoint(environment(getA(ψ), getA(ρ), adjoint.(getA(ψ₀)))),
        FR₀::AbstractVector{<:RightEnvironmentTensor{3}}=_default_X₀_rightFixedPoint(environment(getA(ψ), getA(ρ), adjoint.(getA(ψ₀))))) -> ψ′::DenseInfiniteMPS{L}

Apply iMPO `ρ` to iMPS/iMPO `ψ`, return `ψ′` with maximum overlap with exact result and with D-manifold of `ψ₀` by variational methods.

# Note
For `ρ::SparseUMPO{W, L}` where `W > 1`, `multiply` executes `W` times, contracting only one layer at a time.

# Arguments
`ψ₀::DenseInfiniteMPS{L}`: the initial iMPS/iMPO, with target D-manifold
`ρ::Union{InfiniteMPO{L}, SparseUMPO{W, L}}`: iMPO applying to `ψ`
`ψ::DenseInfiniteMPS{L}`: iMPS/iMPO to find a `ψ′` with certain D-manifold maximizing the overlap with Hψ
`alg_grad::GradientAlgorithm=Defaults.alg_grad`: algorithm for solving gradient problem
`alg_eig::EigenAlgorithm=Defaults.alg_eig`: algorithm for solving eigenvalue problem

# Return
`ψ′::CanonicalMPS`: optimal approximation iMPS with D-manifold of `ψ₀`, maximizing overlap with `ψ`.
`FL::Vector{<:LeftEnvironmentTensor{3}}`: fixed point left environment tensor
`FR::Vector{<:RighttEnvironmentTensor{3}}`: fixed point right environment tensor
`info`: information of the algorithm.
"""
function multiply(ψ₀::DenseInfiniteMPS{L, T},
            ρ::InfiniteMPO{L, T},
            ψ::DenseInfiniteMPS{L, T},
            alg_grad::GradientAlgorithm=Defaults.alg_grad,
            alg_eig::EigenAlgorithm=Defaults.alg_eig;
            FL₀::AbstractVector{<:LeftEnvironmentTensor{3}}=_default_X₀_leftFixedPoint(environment(getA(ψ), getA(ρ), adjoint.(getA(ψ₀)))),
            FR₀::AbstractVector{<:RightEnvironmentTensor{3}}=_default_X₀_rightFixedPoint(environment(getA(ψ), getA(ρ), adjoint.(getA(ψ₀))))) where {L, T}

    # Canonicalize and normalize ψ₀ and ψ
    ψ′ = normalize!(canonicalize(ψ₀))
    ψ_c = normalize!(canonicalize(ψ))

    if alg_grad isa SimpleIterator
        # Define the update process
        func = x -> begin
            ψ_new′ = deepcopy(x)
            # Generate left and right environment
            _, FL, _ = leftFixedPoint(environment(ψ_c.AL, getAL(ρ), adjoint.(ψ_new′.AL)), FL₀, alg_eig)
            _, FR, _ = rightFixedPoint(environment(ψ_c.AR, getAR(ρ), adjoint.(ψ_new′.AR)), FR₀, alg_eig)
            # Update AC and C
            ψ_new′.AC[:] = map(l -> pushmid(ψ_c.AC[l], FL[l], getAC(ρ)[l], FR[l]), 1:L)
            ψ_new′.C[:] = map(l -> pushmid(ψ_c.C[l], FL[mod(l, L) + 1], FR[l]), 1:L)
            # Update AL and AR
            ψ_new′.AL[:] = map(l -> getAL(ψ_new′.AC[l], ψ_new′.C[l]), 1:L)
            ψ_new′.AR[:] = map(l -> getAR(ψ_new′.AC[l], ψ_new′.C[mod(l-2, L) + 1]), 1:L)
            # Update env for next iteration
            FL₀[:] = FL
            FR₀[:] = FR

            ψ_new′
        end
        _, ψ′, info = iterate(func, ψ′, alg_grad; which=:approximate)
        return ψ′, FL₀, FR₀, info
    elseif alg_grad isa OptimKit.OptimizationAlgorithm
        # TODO
        throw(ArgumentError("Undefined behavior"))
    end
end

function multiply(ψ₀::DenseInfiniteMPS{L, T},
            ρ::SparseUMPO{W, L},
            ψ::DenseInfiniteMPS{L, T},
            alg_grad::GradientAlgorithm=Defaults.alg_grad,
            alg_eig::EigenAlgorithm=Defaults.alg_eig;
            FL₀::AbstractVector{<:LeftEnvironmentTensor{3}}=_default_X₀_leftFixedPoint(environment(getA(ψ), getA(ρ)[1:1, :], adjoint.(getA(ψ₀)))),
            FR₀::AbstractVector{<:RightEnvironmentTensor{3}}=_default_X₀_rightFixedPoint(environment(getA(ψ), getA(ρ)[1:1, :], adjoint.(getA(ψ₀))))) where {W, L, T}

    # Canonicalize and normalize ψ₀ and ψ
    ψ′ = normalize!(canonicalize(ψ₀))
    ψ_c = normalize!(canonicalize(ψ))

    if alg_grad isa SimpleIterator
        # Define the update process
        func = [x -> begin
            ψ_new′ = deepcopy(x)
            # Generate left and right environment
            _, FL, _ = leftFixedPoint(environment(ψ_c.AL, getAL(ρ)[w:w, :], adjoint.(ψ_new′.AL)), FL₀, alg_eig)
            _, FR, _ = rightFixedPoint(environment(ψ_c.AR, getAR(ρ)[w:w, :], adjoint.(ψ_new′.AR)), FR₀, alg_eig)
            # Update AC and C
            ψ_new′.AC[:] = map(l -> pushmid(ψ_c.AC[l], FL[l], getAC(ρ)[w, l], FR[l]), 1:L)
            ψ_new′.C[:] = map(l -> pushmid(ψ_c.C[l], FL[mod(l, L) + 1], FR[l]), 1:L)
            # Update AL and AR
            ψ_new′.AL[:] = map(l -> getAL(ψ_new′.AC[l], ψ_new′.C[l]), 1:L)
            ψ_new′.AR[:] = map(l -> getAR(ψ_new′.AC[l], ψ_new′.C[mod(l-2, L) + 1]), 1:L)
            # Update env for next iteration
            FL₀[:] = FL
            FR₀[:] = FR

            ψ_new′
        end for w in 1:W]

        info = SimpleIteratorInfo[]
        for w in 1:W
            _, ψ′, info_w = iterate(func[w], ψ′, alg_grad; which=:approximate)

            ψ_c.AL[:] = ψ′.AL
            ψ_c.AR[:] = ψ′.AR
            ψ_c.AC[:] = ψ′.AC
            ψ_c.C[:] = ψ′.C

            push!(info, info_w)
        end

        return ψ′, FL₀, FR₀, info
    elseif alg_grad isa OptimKit.OptimizationAlgorithm
        # TODO
        throw(ArgumentError("Undefined behavior"))
    end
end
