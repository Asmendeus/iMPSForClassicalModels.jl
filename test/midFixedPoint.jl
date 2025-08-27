using Test, iMPSForClassicalModels

tol = 1e-14
tol1 = 1e-10
tol2 = 1e-8

D = 10
d = 4

@testset "Real BondEnvironment{3}" begin
    FL = LeftEnvironmentTensor(TensorMap(rand, ℝ^D, ℝ^d⊗ℝ^D))
    FR = RightEnvironmentTensor(TensorMap(rand, ℝ^D⊗ℝ^d, ℝ^D))
    env = BondEnvironment(FL, FR)

    λ1, C1, _ = midFixedPoint(env, SimpleIterator(; tol=1e-14))
    λ2, C2, _ = midFixedPoint(env, Arnoldi())
    λ1′, C1′, _ = midFixedPoint(env, C1', SimpleIterator(; tol=1e-14))
    λ2′, C2′, _ = midFixedPoint(env, C2', Arnoldi())

    @test abs(λ1 - λ2) < tol1
    @test norm(C1 - C2) < tol1

    @test abs(λ1′ - λ2′) < tol1
    @test norm(C1′ - C2′) < tol1

    @test abs(λ1 - λ1′) < tol1
end

@testset "Complex BondEnvironment{4}" begin
    FL = LeftEnvironmentTensor(TensorMap(rand, ℂ^D, ℂ^d⊗ℂ^d⊗ℂ^D))
    FR = RightEnvironmentTensor(TensorMap(rand, ℂ^D⊗ℂ^d⊗ℂ^d, ℂ^D))
    env = BondEnvironment(FL, FR)

    λ1, C1, _ = midFixedPoint(env, SimpleIterator(; tol=1e-14))
    λ2, C2, _ = midFixedPoint(env, Arnoldi())
    λ1′, C1′, _ = midFixedPoint(env, C1', SimpleIterator(; tol=1e-14))
    λ2′, C2′, _ = midFixedPoint(env, C2', Arnoldi())

    @test abs(λ1 - λ2) < tol1
    @test norm(C1 - C2) < tol1

    @test abs(λ1′ - λ2′) < tol1
    @test norm(C1′ - C2′) < tol1

    @test abs(λ1 - λ1′) < tol1
end

@testset "Real CenterEnvironment{3}" begin
    FL = LeftEnvironmentTensor(TensorMap(rand, ℝ^D, ℝ^d⊗ℝ^D))
    O = MPOTensor(TensorMap(rand, ℝ^d⊗ℝ^d, ℝ^d⊗ℝ^d))
    FR = RightEnvironmentTensor(TensorMap(rand, ℝ^D⊗ℝ^d, ℝ^D))
    env = CenterEnvironment(FL, O, FR)

    λ1, AC1, _ = midFixedPoint(env, SimpleIterator(; tol=1e-14))
    λ2, AC2, _ = midFixedPoint(env, Arnoldi())
    λ1′, AC1′, _ = midFixedPoint(env, AC1', SimpleIterator(; tol=1e-14))
    λ2′, AC2′, _ = midFixedPoint(env, AC2', Arnoldi())

    @test abs(λ1 - λ2) < tol1
    @test norm(AC1 - AC2) < tol1

    @test abs(λ1′ - λ2′) < tol1
    @test norm(AC1′ - AC2′) < tol1

    @test abs(λ1 - λ1′) < tol1
end

@testset "Complex CenterEnvironment{4}" begin
    FL = LeftEnvironmentTensor(TensorMap(rand, ℂ^D, ℂ^d⊗ℂ^d⊗ℂ^D))
    O1 = MPOTensor(TensorMap(rand, ℂ^d⊗ℂ^d, ℂ^d⊗ℂ^d))
    O2 = MPOTensor(TensorMap(rand, ℂ^d⊗ℂ^d, ℂ^d⊗ℂ^d))
    FR = RightEnvironmentTensor(TensorMap(rand, ℂ^D⊗ℂ^d⊗ℂ^d, ℂ^D))
    env = CenterEnvironment(FL, [O1, O2], FR)

    X₀ = MPOTensor(TensorMap(rand, ℂ^D⊗ℂ^d, ℂ^d⊗ℂ^D))

    λ1, AC1, _ = midFixedPoint(env, X₀, SimpleIterator(; tol=1e-14))
    λ2, AC2, _ = midFixedPoint(env, X₀, Arnoldi())
    λ1′, AC1′, _ = midFixedPoint(env, AC1', SimpleIterator(; tol=1e-14))
    λ2′, AC2′, _ = midFixedPoint(env, AC2', Arnoldi())

    @test abs(λ1 - λ2) < tol1
    # @test norm(AC1 - AC2) < tol1    # Fail easily here because of degeneracy

    @test abs(λ1′ - λ2′) < tol1
    # @test norm(AC1′ - AC2′) < tol1  # Fail easily here because of degeneracy

    @test abs(λ1 - λ1′) < tol1

    @test norm(pushmid(AC1′, env) / λ1′ - AC1′) < tol1
    @test norm(pushmid(AC2′, env) / λ2′ - AC2′) < tol1
end
