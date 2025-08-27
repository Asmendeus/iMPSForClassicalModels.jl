"""
    multiply(ψ₀::DenseInfiniteMPS{L},
        ρ::InfiniteMPO{L},
        ψ::DenseInfiniteMPS{L},
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig;
        FL₀::AbstractVector{<:LeftEnvironmentTensor{3}}=_default_X₀_leftFixedPoint(environment(ψ.AL, ρ.AL, adjoint.(ψ₀.AL))),
        FR₀::AbstractVector{<:RightEnvironmentTensor{3}}=_default_X₀_rightFixedPoint(environment(ψ.AR, ρ.AR, adjoint.(ψ₀.AR)))) -> ψ′::DenseInfiniteMPS{L}
    multiply(ψ₀::DenseInfiniteMPS{L},
        ρ::SparseGeneralMPO{W, L},
        ψ::DenseInfiniteMPS{L},
        alg_grad::GradientAlgorithm=Defaults.alg_grad,
        alg_eig::EigenAlgorithm=Defaults.alg_eig;
        FL₀::AbstractVector{<:LeftEnvironmentTensor{3}}=_default_X₀_leftFixedPoint(environment(ψ.AL, ρ.A[1:1, :], adjoint.(ψ₀.AL))),
        FR₀::AbstractVector{<:RightEnvironmentTensor{3}}=_default_X₀_rightFixedPoint(environment(ψ.AR, ρ.A[1:1, :], adjoint.(ψ₀.AR)))) -> ψ′::DenseInfiniteMPS{L}

Apply iMPO `ρ` to iMPS/iMPO `ψ`, return `ψ′` with maximum overlap with exact result and with D-manifold of `ψ₀` by variational methods.

# Note
For `ρ::SparseGeneralMPO{W, L}` where `W > 1`, `multiply` executes `W` times, contracting only one layer at a time.

# Arguments
`ψ₀::DenseInfiniteMPS{L}`: the initial iMPS/iMPO, with target D-manifold
`ρ::SparseGeneralMPO{L}`: iMPO applying to `ψ`
`ψ::DenseInfiniteMPS{L}`: iMPS/iMPO to find a `ψ′` with certain D-manifold maximizing the overlap with Hψ
`alg_grad::GradientAlgorithm=Defaults.alg_grad`: algorithm for solving gradient problem
`alg_eig::EigenAlgorithm=Defaults.alg_eig`: algorithm for solving eigenvalue problem

# Return
`ψ′::CanonicalMPS`: optimal approximation iMPS with D-manifold of `ψ₀`, maximizing overlap with `ψ`.
`info`: information of the algorithm.
"""
function multiply(ψ₀::DenseInfiniteMPS{L, T},
            ρ::InfiniteMPO{L, T},
            ψ::DenseInfiniteMPS{L, T},
            alg_grad::GradientAlgorithm=Defaults.alg_grad,
            alg_eig::EigenAlgorithm=Defaults.alg_eig;
            FL₀::AbstractVector{<:LeftEnvironmentTensor{3}}=_default_X₀_leftFixedPoint(environment(ψ.AL, reshape(ρ.AL, 1, L), adjoint.(ψ₀.AL))),
            FR₀::AbstractVector{<:RightEnvironmentTensor{3}}=_default_X₀_rightFixedPoint(environment(ψ.AR, reshape(ρ.AR, 1, L), adjoint.(ψ₀.AR)))) where {L, T}

    if alg_grad isa SimpleIterator
        # Define the update process
        func = x -> begin
            ψ_new = deepcopy(x)
            # Generate left and right environment
            _, FL, _ = leftFixedPoint(environment(ψ.AL, reshape(ρ.AL, 1, L), adjoint.(ψ_new.AL)), FL₀, alg_eig)
            _, FR, _ = rightFixedPoint(environment(ψ.AR, reshape(ρ.AR, 1, L), adjoint.(ψ_new.AR)), FR₀, alg_eig)
            # Update AC and C
            ψ_new.AC[:] = map(l -> pushmid(ψ.AC[l], FL[l], ρ.AC[l], FR[l]), 1:L)
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

function multiply(ψ₀::DenseInfiniteMPS{L, T},
            ρ::SparseGeneralMPO{W, L},
            ψ::DenseInfiniteMPS{L, T},
            alg_grad::GradientAlgorithm=Defaults.alg_grad,
            alg_eig::EigenAlgorithm=Defaults.alg_eig;
            FL₀::AbstractVector{<:LeftEnvironmentTensor{3}}=_default_X₀_leftFixedPoint(environment(ψ.AL, ρ.A[1:1, :], adjoint.(ψ₀.AL))),
            FR₀::AbstractVector{<:RightEnvironmentTensor{3}}=_default_X₀_rightFixedPoint(environment(ψ.AR, ρ.A[1:1, :], adjoint.(ψ₀.AR)))) where {W, L, T}

    if alg_grad isa SimpleIterator
        ψ′ = deepcopy(ψ₀)
        ψc = deepcopy(ψ)
        # Define the update process
        func = [x -> begin
            ψ_new = deepcopy(x)
            # Generate left and right environment
            _, FL, _ = leftFixedPoint(environment(ψc.AL, ρ.A[w:w, :], adjoint.(ψ_new.AL)), FL₀, alg_eig)
            _, FR, _ = rightFixedPoint(environment(ψc.AR, ρ.A[w:w, :], adjoint.(ψ_new.AR)), FR₀, alg_eig)
            # Update AC and C
            ψ_new.AC[:] = map(l -> pushmid(ψc.AC[l], FL[l], ρ.A[w, l], FR[l]), 1:L)
            ψ_new.C[:] = map(l -> pushmid(ψc.C[l], FL[mod(l, L) + 1], FR[l]), 1:L)
            # Update AL and AR
            ψ_new.AL[:] = map(l -> getAL(ψ_new.AC[l], ψ_new.C[l]), 1:L)
            ψ_new.AR[:] = map(l -> getAR(ψ_new.AC[l], ψ_new.C[mod(l-2, L) + 1]), 1:L)
            # Update env for next iteration
            FL₀[:] = FL
            FR₀[:] = FR

            ψ_new
        end for w in 1:W]

        info = SimpleIteratorInfo[]
        for w in 1:W
            _, ψ′, info_w = iterate(func[w], ψ′, alg_grad; which=:approximate)
            push!(info, info_w)

            ψc.AL[:] = ψ′.AL
            ψc.AR[:] = ψ′.AR
            ψc.AC[:] = ψ′.AC
            ψc.C[:] = ψ′.C
        end

        return ψ′, info
    elseif alg_grad isa OptimKit.OptimizationAlgorithm
        # TODO
        throw(ArgumentError("Undefined behavior"))
    end
end
