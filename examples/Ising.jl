using TensorKit, iMPSForClassicalModels

# 2D Ising: H = ∑_{ij} - J * s_i s_j
# local tensor of partition function:
#          |
#       —— O ——
#          |
function getLocalTensor(β::Number, J::Number)
    m = sqrt(TensorMap([exp(β * J) exp(- β * J); exp(- β * J) exp(β * J)], ℂ^2, ℂ^2))
    Id = TensorMap(zeros, ℂ^2⊗ℂ^2, ℂ^2⊗ℂ^2)
    @tensor O[-1 -2; -3 -4] := Id[1 2 3 4] * m[-1 1] * m[-2 2] * m[-3 3] * m[-4 4]
    return O
end

T = vcat(0:0.4:1.2, 1.4:0.2:2, 2.1:0.1:2.6, 2.8:0.2:4)
