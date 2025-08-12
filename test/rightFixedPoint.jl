using Test, iMPSForClassicalModels

tol = 1e-14
tol1 = 1e-10
tol2 = 1e-8

D = 10
d = 4

@testset "Real MPSTransferMatrix{1}" begin
    A = MPSTensor(TensorMap(rand, ℝ^D⊗ℝ^d, ℝ^D))
    B = AdjointMPSTensor(TensorMap(rand, ℝ^D, ℝ^D⊗ℝ^d))
    t = TransferMatrix(A, B)

    λ1, L1, _ = rightFixedPoint(t, SimpleIteration(; tol=[tol,]))
    λ2, L2, _ = rightFixedPoint(t, Arnoldi())

    x₀ = RightEnvironmentTensor(TensorMap(rand, ℝ^D, ℝ^D))
    λm, Lm = eigsolve(x->pushright(x, A, B), x₀, 1, :LM)

    @test abs(λ1[1] - λ2[1]) < tol1
    @test norm(L1[1] - L2[1]) < tol1
    @test abs(λ1[1] - λm[1]) < tol1
    @test norm(L1[1] - Lm[1] / sign_first_element(Lm[1])) < tol1
end

@testset "Complex MPOTransferMatrix{2}" begin
    A1 = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))
    A2 = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))
    B1 = AdjointMPOTensor(TensorMap(rand, ℂ^d⊗ℂ^D, ℂ^D⊗ℂ^d))
    B2 = AdjointMPOTensor(TensorMap(rand, ℂ^d⊗ℂ^D, ℂ^D⊗ℂ^d))
    t = TransferMatrix([A1, A2], [B1, B2])

    λ1, L1, _ = rightFixedPoint(t, SimpleIteration(; tol=[tol, tol]))
    λ2, L2, _ = rightFixedPoint(t, Arnoldi())

    x₀ = RightEnvironmentTensor(TensorMap(rand, ℂ^D, ℂ^D))
    λm1, Lm1 = eigsolve(x->pushright(pushright(x, A1, B1), A2, B2), x₀, 1, :LM)
    λm2, Lm2 = eigsolve(x->pushright(pushright(x, A2, B2), A1, B1), x₀, 1, :LM)

    @test abs(prod(λ1) - prod(λ2)) < tol2
    @test abs(prod(λ1) - λm1[1]) < tol2
    @test abs(prod(λ1) - λm2[1]) < tol2
    @test norm(L1[1] - L2[1]) < tol2
    @test norm(L1[2] - L2[2]) < tol2
    @test norm(L1[1] - Lm1[1] / sign_first_element(Lm1[1])) < tol2
    @test norm(L1[2] - Lm2[1] / sign_first_element(Lm2[1])) < tol2
end

@testset "Real right canonicalize (L = 1, R = 3)" begin
    A = MPSTensor(TensorMap(rand, ℝ^D⊗ℝ^d, ℝ^D))
    B = A'

    λ1, L1, AR1, _ = rightFixedPoint([A,], SimpleIteration(; tol=[tol,]))
    λ2, L2, AR2, _ = rightFixedPoint([A,], Arnoldi())
    λ1′, L1′, BR1, _ = rightFixedPoint([B,], SimpleIteration(; tol=[tol,]))
    λ2′, L2′, BR2, _ = rightFixedPoint([B,], Arnoldi())

    @test abs(λ1[1] - λ2[1]) < tol1
    @test norm(L1[1] - L2[1]) < tol1
    @test norm(AR1[1] - AR2[1]) < tol1

    @test abs(λ1[1] - λ1′[1]) < tol1
    @test norm(L1[1] - L1′[1]') < tol1
    @test norm(AR1[1] - BR1[1]') < tol1

    @test abs(λ1′[1] - λ2′[1]) < tol1
    @test norm(L1′[1] - L2′[1]) < tol1
    @test norm(BR1[1] - BR2[1]) < tol1
end

@testset "Complex right canonicalize (L = 2, R = 4)" begin
    A1 = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))
    A2 = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))
    B1 = A1'
    B2 = A2'

    λ1, L1, AR1, _ = rightFixedPoint([A1, A2], SimpleIteration(; tol=[tol, tol]))
    λ2, L2, AR2, _ = rightFixedPoint([A1, A2], Arnoldi())
    λ1′, L1′, BR1, _ = rightFixedPoint([B1, B2], SimpleIteration(; tol=[tol, tol]))
    λ2′, L2′, BR2, _ = rightFixedPoint([B1, B2], Arnoldi())

    @test abs(λ1[1] - λ2[1]) < tol2
    @test abs(λ1[2] - λ2[2]) < tol2
    @test norm(L1[1] - L2[1]) < tol2
    @test norm(L1[2] - L2[2]) < tol2
    @test norm(AR1[1] - AR2[1]) < tol2
    @test norm(AR1[1] - AR2[1]) < tol2

    @test abs(λ1[1] - λ1[1]) < tol2
    @test abs(λ1[2] - λ1[2]) < tol2
    @test norm(L1[1] - L1′[1]') < tol2
    @test norm(L1[2] - L1′[2]') < tol2
    @test norm(AR1[1] - BR1[1]') < tol2
    @test norm(AR1[1] - BR1[1]') < tol2

    @test abs(λ1′[1] - λ2′[1]) < tol2
    @test abs(λ1′[2] - λ2′[2]) < tol2
    @test norm(L1′[1] - L2′[1]) < tol2
    @test norm(L1′[2] - L2′[2]) < tol2
    @test norm(BR1[1] - BR2[1]) < tol2
    @test norm(BR1[1] - BR2[1]) < tol2
end
