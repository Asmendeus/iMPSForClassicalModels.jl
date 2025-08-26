"""
    approximate(ψ₀::DenseInfiniteMPS{L, T},
        ψ::DenseInfiniteMPS{L, T},
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig;
        FL₀::AbstractVector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(TransferMatrix(getA(ψ), adjoint.(getA(ψ₀)))),
        FR₀::AbstractVector{<:RightEnvironmentTensor{2}}=_default_X₀_rightFixedPoint(TransferMatrix(getA(ψ), adjoint.(getA(ψ₀)))))) -> ψ′::DenseInfiniteMPS{L}

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
            alg_eig::EigenAlgorithm=Defaults.alg_eig;
            FL₀::AbstractVector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(TransferMatrix(getAL(ψ), adjoint.(getAL(ψ₀)))),
            FR₀::AbstractVector{<:RightEnvironmentTensor{2}}=_default_X₀_rightFixedPoint(TransferMatrix(getAR(ψ), adjoint.(getAR(ψ₀))))) where {L, T}

    # Canonicalize and normalize ψ₀ and ψ
    ψ′ = normalize!(canonicalize(ψ₀))
    ψ_c = normalize!(canonicalize(ψ))

    if alg_grad isa SimpleIterator
        # Define the update process
        func = x -> begin
            ψ_new′ = deepcopy(x)
            # Generate left and right environment
            _, FL, _ = leftFixedPoint(TransferMatrix(ψ_c.AL, adjoint.(ψ_new′.AL)), FL₀, alg_eig)
            _, FR, _ = rightFixedPoint(TransferMatrix(ψ_c.AR, adjoint.(ψ_new′.AR)), FR₀, alg_eig)
            # Update AC and C
            ψ_new′.AC[:] = map(l -> pushmid(ψ_c.AC[l], FL[l], FR[l]), 1:L)
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
        return ψ′, info
    elseif alg_grad isa OptimKit.OptimizationAlgorithm
        # TODO
        throw(ArgumentError("Undefined behavior"))
    end
end
