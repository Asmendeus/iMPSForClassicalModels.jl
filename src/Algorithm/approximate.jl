"""
    approximate(ψ₀::DenseInfiniteMPS{L, T},
        ψ::DenseInfiniteMPS{L, T},
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig;
        FL₀::AbstractVector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(TransferMatrix(ψ.AL, adjoint.(ψ₀.AL))),
        FR₀::AbstractVector{<:RightEnvironmentTensor{2}}=_default_X₀_rightFixedPoint(TransferMatrix(ψ.AR, adjoint.(ψ₀.AR))))) -> ψ′::DenseInfiniteMPS{L}

Generate the optimal approximation `ψ′` with D-manifold of `ψ₀`, maximizing overlap with `ψ` by variational methods.

# Arguments
`ψ₀::DenseInfiniteMPS{L}`: the initial iMPS/iMPO with target D-manifold
`ψ::DenseInfiniteMPS{L}`: iMPS/iMPO to find a `ψ′` with certain D-manifold
`alg_grad::GradientAlgorithm=Defaults.alg_grad`: algorithm for solving gradient problem
`alg_eig::EigenAlgorithm=Defaults.alg_eig`: algorithm for solving eigenvalue problem

# Keyword arguments
`FL₀::AbstractVector{<:LeftEnvironmentTensor{2}}`: left environment
`FR₀::AbstractVector{<:RightEnvironmentTensor{2}}`: right environment

# Return
`ψ′::DenseInfiniteMPS{L}`: optimal approximation iMPS with D-manifold of `ψ₀`, maximizing overlap with `ψ`.
`info`: information of the algorithm.
"""
function approximate(ψ₀::DenseInfiniteMPS{L, T},
            ψ::DenseInfiniteMPS{L, T},
            alg_grad::GradientAlgorithm=Defaults.alg_grad,
            alg_eig::EigenAlgorithm=Defaults.alg_eig;
            FL₀::AbstractVector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(TransferMatrix(ψ.AL, adjoint.(ψ₀.AL))),
            FR₀::AbstractVector{<:RightEnvironmentTensor{2}}=_default_X₀_rightFixedPoint(TransferMatrix(ψ.AR, adjoint.(ψ₀.AR)))) where {L, T}

    # Canonicalize and normalize ψ₀ and ψ
    if alg_grad isa SimpleIterator
        # Define the update process
        func = x -> begin
            ψ_new = deepcopy(x)
            # Generate left and right environment
            _, FL, _ = leftFixedPoint(TransferMatrix(ψ.AL, adjoint.(ψ_new.AL)), FL₀, alg_eig)
            _, FR, _ = rightFixedPoint(TransferMatrix(ψ.AR, adjoint.(ψ_new.AR)), FR₀, alg_eig)
            # Update AC and C
            ψ_new.AC[:] = map(l -> pushmid(ψ.AC[l], FL[l], FR[l]), 1:L)
            ψ_new.C[:] = map(l -> pushmid(ψ.C[l], FL[mod(l, L) + 1], FR[l]), 1:L)
            # Update AL and AR
            ψ_new.AL[:] = map(l -> getAL(ψ_new.AC[l], ψ_new.C[l]), 1:L)
            ψ_new.AR[:] = map(l -> getAR(ψ_new.AC[l], ψ_new.C[mod(l-2, L) + 1]), 1:L)
            # Update env for next iteration
            FL₀[:] = FL
            FR₀[:] = FR

            ψ_new
        end
        _, ψ′, info = iterate(func, ψ₀, alg_grad; which=:approximate)
        return ψ′, info
    elseif alg_grad isa OptimKit.OptimizationAlgorithm
        # TODO
        throw(ArgumentError("Undefined behavior"))
    end
end
