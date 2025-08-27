using Test, iMPSForClassicalModels

tol = 1e-8

D = 10
d = 4

@testset "Real left canonicalize (L = 1, R = 3)" begin
    A = MPSTensor(TensorMap(rand, ℝ^D⊗ℝ^d, ℝ^D))
    B = A'

    λ1, L1, AL1, _ = leftFixedPoint([A,], SimpleIterator())
    λ2, L2, AL2, _ = leftFixedPoint([A,], Arnoldi())
    λ1′, L1′, BL1, _ = leftFixedPoint([B,], SimpleIterator())
    λ2′, L2′, BL2, _ = leftFixedPoint([B,], Arnoldi())

    @test abs(λ1[1] - λ2[1]) < tol
    @test norm(L1[1] - L2[1]) < tol
    @test norm(AL1[1] - AL2[1]) < tol

    @test abs(λ1[1] - λ1′[1]) < tol
    @test norm(L1[1] - L1′[1]') < tol
    @test norm(AL1[1] - BL1[1]') < tol

    @test abs(λ1′[1] - λ2′[1]) < tol
    @test norm(L1′[1] - L2′[1]) < tol
    @test norm(BL1[1] - BL2[1]) < tol
end

@testset "Complex left canonicalize (L = 2, R = 4)" begin
    A1 = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))
    A2 = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))
    B1 = A1'
    B2 = A2'

    λ1, L1, AL1, _ = leftFixedPoint([A1, A2], SimpleIterator(; tol=1e-14))
    λ2, L2, AL2, _ = leftFixedPoint([A1, A2], Arnoldi())
    λ1′, L1′, BL1, _ = leftFixedPoint([B1, B2], SimpleIterator(; tol=1e-14))
    λ2′, L2′, BL2, _ = leftFixedPoint([B1, B2], Arnoldi())

    @test abs(λ1[1] - λ2[1]) < tol
    @test abs(λ1[2] - λ2[2]) < tol
    @test norm(L1[1] - L2[1]) < tol
    @test norm(L1[2] - L2[2]) < tol
    @test norm(AL1[1] - AL2[1]) < tol
    @test norm(AL1[2] - AL2[2]) < tol

    @test abs(λ1[1] - λ2′[1]) < tol
    @test abs(λ1[2] - λ2′[2]) < tol
    @test norm(L1[1] - L2′[1]') < tol
    @test norm(L1[2] - L2′[2]') < tol
    @test norm(AL1[1] - BL2[1]') < tol
    @test norm(AL1[2] - BL2[2]') < tol

    @test abs(λ1′[1] - λ2′[1]) < tol
    @test abs(λ1′[2] - λ2′[2]) < tol
    @test norm(L1′[1]' - L2′[1]') < tol
    @test norm(L1′[2]' - L2′[2]') < tol
    @test norm(BL1[1] - BL2[1]) < tol
    @test norm(BL1[2] - BL2[2]) < tol
end

@testset "Real MPSTransferMatrix{1}" begin
    A = MPSTensor(TensorMap(rand, ℝ^D⊗ℝ^d, ℝ^D))
    B = AdjointMPSTensor(TensorMap(rand, ℝ^D, ℝ^D⊗ℝ^d))
    t = TransferMatrix(A, B)

    λ1, L1, _ = leftFixedPoint(t, SimpleIterator(; tol=1e-14))
    λ2, L2, _ = leftFixedPoint(t, Arnoldi())

    @test abs(λ1[1] - λ2[1]) < tol
    @test norm(L1[1] - L2[1]) < tol

    x₀ = LeftEnvironmentTensor(TensorMap(rand, ℝ^D, ℝ^D))
    λm, Lm = eigsolve(x->pushleft(x, A, B), x₀, 1, :LM)
    @test abs(λ1[1] - λm[1]) < tol
    @test norm(L1[1] - Lm[1]/sign(Lm[1].A[1])) < tol
end

@testset "Complex MPOTransferMatrix{2}" begin
    A1 = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))
    A2 = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))
    B1 = AdjointMPOTensor(TensorMap(rand, ℂ^d⊗ℂ^D, ℂ^D⊗ℂ^d))
    B2 = AdjointMPOTensor(TensorMap(rand, ℂ^d⊗ℂ^D, ℂ^D⊗ℂ^d))
    t = TransferMatrix([A1, A2], [B1, B2])

    λ1, L1, _ = leftFixedPoint(t, SimpleIterator(; tol=1e-14))
    λ2, L2, _ = leftFixedPoint(t, Arnoldi())

    @test abs(λ1[1] - λ2[1]) < tol
    @test abs(λ1[2] - λ2[2]) < tol
    @test norm(L1[1] - L2[1]) < tol
    @test norm(L1[2] - L2[2]) < tol

    x₀ = LeftEnvironmentTensor(TensorMap(rand, ℂ^D, ℂ^D))
    λm1, Lm1 = eigsolve(x->pushleft(pushleft(x, A1, B1), A2, B2), x₀, 1, :LM)
    λm2, Lm2 = eigsolve(x->pushleft(pushleft(x, A2, B2), A1, B1), x₀, 1, :LM)
    @test abs(prod(λ1) - λm1[1]) < tol
    @test abs(prod(λ1) - λm2[1]) < tol
    @test norm(L1[1] - Lm1[1]/sign(Lm1[1].A[1])) < tol
    @test norm(L1[2] - Lm2[1]/sign(Lm2[1].A[1])) < tol
end

@testset "Real left ChannelEnvironment{3, 1, 3}" begin
    A = MPSTensor(TensorMap(rand, ℝ^D⊗ℝ^d, ℝ^D))
    O = MPOTensor(TensorMap(rand, ℝ^d⊗ℝ^d, ℝ^d⊗ℝ^d))
    B = AdjointMPSTensor(TensorMap(rand, ℝ^D, ℝ^D⊗ℝ^d))

    env = ChannelEnvironment(A, O, B)

    λ1, FL1, _ = leftFixedPoint(env, SimpleIterator(; tol=1e-14))
    λ2, FL2, _ = leftFixedPoint(env, Arnoldi())

    @test abs(λ1[1] - λ2[1]) < tol
    @test norm(FL1[1] - FL2[1]) < tol

    @test norm(pushleft(FL1[1], A, O, B) - λ1[1] * FL1[1]) < tol
end

@testset "Complex left ChannelEnvironment{4, 2, 4}" begin
    A1 = MPSTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^D))
    A2 = MPSTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^D))
    O11 = MPOTensor(TensorMap(rand, ℂ^d⊗ℂ^d, ℂ^d⊗ℂ^d))
    O12 = MPOTensor(TensorMap(rand, ℂ^d⊗ℂ^d, ℂ^d⊗ℂ^d))
    O21 = MPOTensor(TensorMap(rand, ℂ^d⊗ℂ^d, ℂ^d⊗ℂ^d))
    O22 = MPOTensor(TensorMap(rand, ℂ^d⊗ℂ^d, ℂ^d⊗ℂ^d))
    B1 = AdjointMPSTensor(TensorMap(rand, ℂ^D, ℂ^D⊗ℂ^d))
    B2 = AdjointMPSTensor(TensorMap(rand, ℂ^D, ℂ^D⊗ℂ^d))

    env = ChannelEnvironment([A1, A2], [O11 O12; O21 O22], [B1, B2])

    λ1, FL1, _ = leftFixedPoint(env, SimpleIterator(; tol=1e-14))
    λ2, FL2, _ = leftFixedPoint(env, Arnoldi())

    @test abs(λ1[1] - λ2[1]) < tol
    @test abs(λ1[2] - λ2[2]) < tol
    @test norm(FL1[1] - FL2[1]) < tol
    @test norm(FL1[2] - FL2[2]) < tol

    @test norm(pushleft(FL2[1], A1, [O11, O21], B1) - λ2[1] * FL2[2]) < tol
    @test norm(pushleft(FL2[2], A2, [O12, O22], B2) - λ2[2] * FL2[1]) < tol
    @test norm(pushleft(FL2[1], env) - λ2[1] * λ2[2] * FL2[1]) < tol
end
