"""
    approximate(ψ₀::DenseInfiniteMPS{L, T},
        ψ::DenseInfiniteMPS{L, T},
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig) -> ψ′::DenseInfiniteMPS{L}

Generate the optimal approximation `ψ′` with D-manifold of `ψ₀`, maximizing overlap with `ψ` by variational methods.

# Arguments
`ψ₀::DenseInfiniteMPS{L}`: the initial iMPS/iMPO with target D-manifold
`ψ::DenseInfiniteMPS{L}`: iMPS/iMPO to find a `ψ′` with certain D-manifold
`alg_grad::GradientAlgorithm=Defaults.alg_grad`: algorithm for solving gradient problem
`alg_eig::EigenAlgorithm=Defaults.alg_eig`: algorithm for solving eigenvalue problem

# Return
`ψ′::CanonicalMPS`: optimal approximation iMPS with D-manifold of `ψ₀`, maximizing overlap with `ψ`.
`info`: information of the algorithm.
"""
function approximate(ψ₀::DenseInfiniteMPS{L, T},
            ψ::DenseInfiniteMPS{L, T},
            alg_grad::GradientAlgorithm=Defaults.alg_grad,
            alg_eig::EigenAlgorithm=Defaults.alg_eig) where {L, T}

    # Canonicalize and normalize ψ₀ and ψ
    ψ′ = normalize!(canonicalize(ψ₀))
    ψ_c = normalize!(canonicalize(ψ))

    if alg_grad isa SimpleIterator
        # Define the update process
        func = x -> begin
            ψ_new′ = deepcopy(x)
            # Generate left and right environment
            _, FL, _ = leftFixedPoint(TransferMatrix(ψ_c.AL, adjoint.(ψ_new′.AL)), alg_eig)
            _, FR, _ = rightFixedPoint(TransferMatrix(ψ_c.AR, adjoint.(ψ_new′.AR)), alg_eig)
            # Update AC and C
            ψ_new′.AC[:] = map(l -> pushmid(ψ_c.AC[l], FL[l], FR[l]), 1:L)
            ψ_new′.C[:] = map(l -> pushmid(ψ_c.C[l], FL[mod(l, L) + 1], FR[l]), 1:L)
            # Update AL and AR
            ψ_new′.AL[:] = map(l -> getAL(ψ_new′.AC[l], ψ_new′.C[l]), 1:L)
            ψ_new′.AR[:] = map(l -> getAR(ψ_new′.AC[l], ψ_new′.C[mod(l-2, L) + 1]), 1:L)

            ψ_new′
        end
        _, ψ′, info = iterate(func, ψ′, alg_grad)
        return ψ′, info
    elseif alg_grad isa OptimKit.OptimizationAlgorithm
        # TODO
        throw(ArgumentError("Undefined behavior"))
    end
end
