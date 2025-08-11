using Test, iMPSForClassicalModels

tol1 = 1e-10
tol2 = 1e-8

D = 10
d = 4

@testset "Real MPSTransferMatrix{1}" begin
    A = MPSTensor(TensorMap(rand, ℝ^D⊗ℝ^d, ℝ^D))
    B = AdjointMPSTensor(TensorMap(rand, ℝ^D, ℝ^D⊗ℝ^d))
    t = TransferMatrix(A, B)

    λ, L, _ = leftFixedPoint(t)

    x₀ = LeftEnvironmentTensor(TensorMap(rand, ℝ^D, ℝ^D))
    λm, Lm = eigsolve(x->pushleft(x, A, B), x₀, 1, :LM)

    @test abs(λ[1] - λm[1]) < tol
    @test norm(L[1] - Lm[1] / sign_first_element(Lm[1])) < tol1
end

@testset "Complex MPOTransferMatrix{2}" begin
    A1 = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))
    A2 = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))
    B1 = AdjointMPOTensor(TensorMap(rand, ℂ^d⊗ℂ^D, ℂ^D⊗ℂ^d))
    B2 = AdjointMPOTensor(TensorMap(rand, ℂ^d⊗ℂ^D, ℂ^D⊗ℂ^d))
    t = TransferMatrix([A1, A2], [B1, B2])

    λ1, L1, _ = leftFixedPoint(t)
    λ2, L2, _ = leftFixedPoint(t, Arnoldi())

    x₀ = LeftEnvironmentTensor(TensorMap(rand, ℂ^D, ℂ^D))
    λm1, Lm1 = eigsolve(x->pushleft(pushleft(x, A1, B1), A2, B2), x₀, 1, :LM)
    λm2, Lm2 = eigsolve(x->pushleft(pushleft(x, A2, B2), A1, B1), x₀, 1, :LM)

    @test abs(prod(λ1) - prod(λ2)) < tol2
    @test abs(prod(λ1) - λm1[1]) < tol2
    @test abs(prod(λ1) - λm2[1]) < tol2
    @test norm(L1[1] - L2[1]) < tol2
    @test norm(L1[2] - L2[2]) < tol2
    @test norm(L1[1] - Lm1[1] / sign_first_element(Lm1[1])) < tol2
    @test norm(L1[2] - Lm2[1] / sign_first_element(Lm2[1])) < tol2
end
