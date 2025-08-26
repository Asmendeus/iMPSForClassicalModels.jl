function fixedBondary(ψ::DenseInfiniteMPS{L, T}, ρ::SparseUMPO{W, L}, alg::VUMPS) where {W, L, T}
end

"""
    fixedBondary(ψ::DenseInfiniteMPS{L, T}, ρ::SparseUMPO{W, L}, alg::ViTEBD)

Slove the following maximum eigenvalue equation by ViTEBD algorithm

    ... —— A[1] —— ... —— A[L] —— ...  =  λ * ... —— A[1] —— ... —— A[L] —— ...
           |              |                          |              |
    ... —— O[1] —— ... —— O[L] —— ...
           |              |

# Arguments
`ψ::DenseInfiniteMPS{L, T}`: initial variational boundary iMPS
`ρ::SparseUMPO{W, L}`: partition function of tensor network representation
"""
function fixedBondary(ψ::DenseInfiniteMPS{L, T}, ρ::SparseUMPO{W, L}, alg::ViTEBD) where {W, L, T}
    if alg.alg_grad isa SimpleIterator
        func = x -> multiply(x, ρ, x, alg.alg_grad_mutiply, alg_eig_mutiply)[1]
        _, ψ′, info = iterate(func, ψ′, alg_grad; which=:ViTEBD)
        return ψ′, info
    elseif alg.alg_grad isa OptimKit.OptimizationAlgorithm
        # TODO
        throw(ArgumentError("Undefined behavior"))
    end
end
