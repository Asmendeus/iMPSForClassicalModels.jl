"""
    fixedBondary(ψ::DenseInfiniteMPS{L, T}, ρ::SparseGeneralMPO{W, L}, alg::ViTEBD)

Slove the following maximum eigenvalue equation by VUMPS algorithm

    ... —— A[1] —— ... —— A[L] —— ...  =  λ * ... —— A[1] —— ... —— A[L] —— ...
           |              |                          |              |
    ... —— O[1] —— ... —— O[L] —— ...
           |              |

# Arguments
`ψ::DenseInfiniteMPS{L, T}`: initial variational boundary iMPS
`ρ::SparseGeneralMPO{W, L}`: partition function of tensor network representation
"""
function fixedBondary(ψ::DenseInfiniteMPS{L, T}, ρ::SparseGeneralMPO{W, L}, alg::VUMPS;
            FL₀::AbstractVector{<:LeftEnvironmentTensor{N}}=_default_X₀_leftFixedPoint(environment(ψ.AL, ρ.A, adjoint.(ψ.AL))),
            FR₀::AbstractVector{<:RightEnvironmentTensor{N}}=_default_X₀_rightFixedPoint(environment(ψ.AR, ρ.A, adjoint.(ψ.AR)))) where {W, L, N, T}

    N == W + 2 || throw(ArgumentError("Dismatched `N` and `W + 2`: $N ≠ $(W+2)"))

    if alg.alg_grad isa SimpleIterator
        # Define the update process
        func = x -> begin
            ψ_new = deepcopy(x)
            # Generate left and right environment
            _, FL, _ = leftFixedPoint(environment(ψ_new.AL, ρ.A, adjoint.(ψ_new.AL)), FL₀, alg.alg_eig_environment)
            _, FR, _ = rightFixedPoint(environment(ψ_new.AR, ρ.A, adjoint.(ψ_new.AR)), FR₀, alg.alg_eig_environment)
            # Update AC and C
            tmp_AC = map(l -> midFixedPoint(environment(FL[l], ρ.A[:, l], FR[l]))[2], 1:L)
            ψ_new.AC[:] = tmp_AC
            tmp_C = map(l -> midFixedPoint(environment(FL[l], FR[l]))[2], 1:L)
            ψ_new.C[:] = tmp_C
            # Update AL and AR
            ψ_new.AL[:] = map(l -> getAL(ψ_new.AC[l], ψ_new.C[l]), 1:L)
            ψ_new.AR[:] = map(l -> getAR(ψ_new.AC[l], ψ_new.C[mod(l-2, L) + 1]), 1:L)
            # Update env for next iteration
            FL₀[:] = FL
            FR₀[:] = FR

            ψ_new
        end
        _, ψ′, info = iterate(func, ψ, alg.alg_grad; which=:VUMPS)
        return ψ′, info
    elseif alg.alg_grad isa OptimKit.OptimizationAlgorithm
        # TODO
        throw(ArgumentError("Undefined behavior"))
    end
end

"""
    fixedBondary(ψ::DenseInfiniteMPS{L, T}, ρ::SparseGeneralMPO{W, L}, alg::ViTEBD)

Slove the following maximum eigenvalue equation by ViTEBD algorithm

    ... —— A[1] —— ... —— A[L] —— ...  =  λ * ... —— A[1] —— ... —— A[L] —— ...
           |              |                          |              |
    ... —— O[1] —— ... —— O[L] —— ...
           |              |

# Arguments
`ψ::DenseInfiniteMPS{L, T}`: initial variational boundary iMPS
`ρ::SparseGeneralMPO{W, L}`: partition function of tensor network representation
"""
function fixedBondary(ψ::DenseInfiniteMPS{L, T}, ρ::SparseGeneralMPO{W, L}, alg::ViTEBD) where {W, L, T}
    if alg.alg_grad isa SimpleIterator
        func = x -> multiply(x, ρ, x, alg.alg_grad_mutiply, alg.alg_eig_mutiply)[1]
        _, ψ′, info = iterate(func, ψ, alg.alg_grad; which=:ViTEBD)
        return ψ′, info
    elseif alg.alg_grad isa OptimKit.OptimizationAlgorithm
        # TODO
        throw(ArgumentError("Undefined behavior"))
    end
end
