"""
    approximate(ψ::DenseInfiniteMPS{L},
        ψ₀::DenseInfiniteMPS{L},
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig) -> ψ′::DenseInfiniteMPS{L}

Generate the optimal approximation `ψ′` with maximum overlap with `ψ` and with D-manifold of `ψ₀` by variational methods.

`ψ::DenseInfiniteMPS{L}`: iMPS/iMPO to find a `ψ′` with certain D-manifold maximizing |⟨ψ′|ψ⟩|^2 / ⟨ψ′|ψ′⟩ or |Tr(ψ′†ψ)|^2 / Tr(ψ′†ψ′)
`ψ₀::DenseInfiniteMPS{L}`: the initial iMPS/iMPO, with target D-manifold
`alg_grad::GradientAlgorithm=Defaults.alg_grad`: algorithm for solving gradient problem
`alg_eig::EigenAlgorithm=Defaults.alg_eig`: algorithm for solving eigenvalue problem
"""
function approximate(ψ::DenseInfiniteMPS{L},
        ψ₀::DenseInfiniteMPS{L},
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig) where L
    
end
