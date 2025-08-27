using TensorKit, iMPSForClassicalModels
using Plots

# 2D Ising: H = ∑_{ij} - J * s_i s_j
# local tensor of partition function:
#          |
#       —— O ——
#          |

# 2D Ising transformation temperature: Tc = 2.269 J
function getLocalTensor(β::Number, J::Number)
    m = sqrt(TensorMap([exp(β * J) exp(- β * J); exp(- β * J) exp(β * J)], ℂ^2, ℂ^2))
    Id = TensorMap(zeros, ℂ^2⊗ℂ^2, ℂ^2⊗ℂ^2)
    Id[1,1,1,1] = Id[2,2,2,2] = 1.0
    @tensor O[-1 -2; -3 -4] := Id[1 2 3 4] * m[-1 1] * m[-2 2] * m[3 -3] * m[4 -4]
    return MPOTensor(O)
end

T_list = vcat(0.4:0.4:1.2, 1.4:0.2:2, 2.1:0.1:2.6, 2.8:0.2:4)
β_list = 1 ./ T_list

EE = Float64[]

ψ = randInfiniteMPS(ComplexF64, 1, 2, 100)

for β in β_list
    O = getLocalTensor(β, 1.0)
    ρ = SparseGeneralMPO(O)
    ψ′, _ = fixedBondary(ψ, ρ, VUMPS())
    push!(EE, entropy(ψ′)[1])

    ψ.AL[:] = ψ′.AL[:]
    ψ.AR[:] = ψ′.AR[:]
    ψ.AC[:] = ψ′.AC[:]
    ψ.C[:] = ψ′.C[:]
end

plot(T_list, EE,
    title="Entropy of entanglement of 2D classical Ising model",
    xlabel="temperature: T/J",
    ylabel="entropy of entanglement: EE",
    label="")
savefig("iMPSForClassicalModels/examples/2DIsing.png")
