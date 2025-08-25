"""
    approximate!(ψ₀::DenseInfiniteMPS{L},
        ψ::DenseInfiniteMPS{L},
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
"""
function approximate!(ψ₀::DenseInfiniteMPS{L},
            ψ::DenseInfiniteMPS{L},
            alg_grad::GradientAlgorithm=Defaults.alg_grad,
            alg_eig::EigenAlgorithm=Defaults.alg_eig) where L

    # Canonicalize and normalize ψ₀ and ψ
    ψ′ = normalize!(canonicalize(ψ₀))
    ψ_c = normalize!(canonicalize(ψ))

    if alg_grad isa SimpleIterator
        # Define the update process
        func = [x -> begin
            # Generate left and right environment
            
            # Update AC and C
            # Update AL and AR
        end for l in 1:L]
    elseif alg_grad isa OptimKit.OptimizationAlgorithm
        # TODO
        throw(ArgumentError("Undefined behavior"))
    end
end
