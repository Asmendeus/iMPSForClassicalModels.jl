"""
    multiply(ρ::InfiniteMPO{L},
        ψ::DenseInfiniteMPS{L},
        ψ₀::DenseInfiniteMPS{L},
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig) -> ψ′::DenseInfiniteMPS{L}
    multiply(ρ::SparseUMPO{W, L},
        ψ::DenseInfiniteMPS{L},
        ψ₀::DenseInfiniteMPS{L},
        maxlayer::Int64=1,
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig) -> ψ′::DenseInfiniteMPS{L}

Apply iMPO `ρ` to iMPS/iMPO `ψ`, return `ψ′` with maximum overlap with exact result and with D-manifold of `ψ₀` by variational methods.

# Arguments
`ρ::Union{InfiniteMPO{L}, SparseUMPO{W, L}}`: iMPO applying to `ψ`
`ψ::DenseInfiniteMPS{L}`: iMPS/iMPO to find a `ψ′` with certain D-manifold maximizing the overlap with Hψ
`ψ₀::DenseInfiniteMPS{L}`: the initial iMPS/iMPO, with target D-manifold
`alg_grad::GradientAlgorithm=Defaults.alg_grad`: algorithm for solving gradient problem
`alg_eig::EigenAlgorithm=Defaults.alg_eig`: algorithm for solving eigenvalue problem
"""
function multiply(ρ::InfiniteMPO{L},
            ψ::DenseInfiniteMPS{L},
            ψ₀::DenseInfiniteMPS{L},
            alg_grad::GradientAlgorithm=Defaults.alg_grad,
            alg_eig::EigenAlgorithm=Defaults.alg_eig) where L
end
function multiply(ρ::SparseUMPO{W, L},
            ψ::DenseInfiniteMPS{L},
            ψ₀::DenseInfiniteMPS{L},
            alg_grad::GradientAlgorithm=Defaults.alg_grad,
            alg_eig::EigenAlgorithm=Defaults.alg_eig) where {W, L}
end
