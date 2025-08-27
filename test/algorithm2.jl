using TensorKit, Test, iMPSForClassicalModels

tol = 1e-8

function getLocalTensor(β::Number, J::Number)
    m = sqrt(TensorMap([exp(β * J) exp(- β * J); exp(- β * J) exp(β * J)], ℂ^2, ℂ^2))
    Id = TensorMap(zeros, ℂ^2⊗ℂ^2, ℂ^2⊗ℂ^2)
    Id[1,1,1,1] = Id[2,2,2,2] = 1.0
    @tensor O[-1 -2; -3 -4] := Id[1 2 3 4] * m[-1 1] * m[-2 2] * m[3 -3] * m[4 -4]
    return MPOTensor(O)
end

O = getLocalTensor(0.5, 1.0)

ρ1 = SparseGeneralMPO(O)
ρ2 = SparseGeneralMPO([O; O;;])

ψ₀ = randInfiniteMPS(ComplexF64, 1, 2, 100)

ψ1, _ = fixedBondary(ψ₀, ρ1, VUMPS())
ψ2, _ = fixedBondary(ψ₀, ρ1, ViTEBD())
ψ3, _ = fixedBondary(ψ₀, ρ2, ViTEBD())

@testset "2D Ising (T = 2 * J)" begin
    @test abs(log(overlap(ψ1, ψ2))) < tol
    @test abs(log(overlap(ψ2, ψ3))) < tol
    @test abs(log(overlap(ψ3, ψ1))) < tol

    @test abs(entropy(ψ1)[1] - entropy(ψ2)[1]) < tol
end
