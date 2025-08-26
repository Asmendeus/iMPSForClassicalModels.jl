using Test, iMPSForClassicalModels

D1 = 80
D2 = 100
d = 4

tol = 1e-8
tol_low = 1e-6

@testset "approximate(iMPS{1, ComplexF64})" begin
    ψ = randCMPS(ComplexF64, 1, d, D2)

    ψ1₀ = randUMPS(ComplexF64, 1, d, D1)
    ψ1′, info = approximate(ψ1₀, ψ)
    @test abs(log(overlap(ψ, ψ1′))) < tol
    if !info.converged
        @show info.normres
    end

    ψ2₀ = randUMPS(ComplexF64, 1, d, D2)
    ψ2′, info = approximate(ψ2₀, ψ)
    @test abs(log(overlap(ψ, ψ2′))) < tol
    if !info.converged
        @show info.normres
    end
end

@testset "approximate(iMPS{2, ComplexF64})" begin
    ψ = randCMPS(ComplexF64, 2, d, D2)

    ψ1₀ = randUMPS(ComplexF64, 2, d, D1)
    ψ1′, info = approximate(ψ1₀, ψ)
    @test abs(log(overlap(ψ, ψ1′))) < tol
    if !info.converged
        @show info.normres
    end

    ψ2₀ = randUMPS(ComplexF64, 2, d, D2)
    ψ2′, info = approximate(ψ2₀, ψ)
    @test abs(log(overlap(ψ, ψ2′))) < tol
    if !info.converged
        @show info.normres
    end
end

@testset "approximate(iMPO{1, ComplexF64})" begin
    O1 = MPOTensor(TensorMap(rand, ℂ^D1⊗ℂ^d, ℂ^d⊗ℂ^D1))
    O2 = MPOTensor(TensorMap(rand, ℂ^D2⊗ℂ^d, ℂ^d⊗ℂ^D2))
    O2′ = MPOTensor(TensorMap(rand, ℂ^D2⊗ℂ^d, ℂ^d⊗ℂ^D2))

    ρ = UMPO(O2)

    ρ1₀ = UMPO(O1)
    ρ1′, info = approximate(ρ1₀, ρ)
    @test abs(log(overlap(ρ, ρ1′))) < tol
    @show overlap(ρ, ρ1′)
    if !info.converged
        @show info.normres
    end

    ρ2₀ = UMPO(O2′)
    ρ2′, info = approximate(ρ2₀, ρ)
    @test abs(log(overlap(ρ, ρ2′))) < tol
    if !info.converged
        @show info.normres
    end
end

@testset "multiply(iMPS{1, ComplexF64})" begin
    ψ = randCMPS(ComplexF64, 1, d, D1)
    ψ₀ = randCMPS(ComplexF64, 1, d, D1)
    O = MPOTensor(TensorMap(rand, ℂ^D1⊗ℂ^d, ℂ^d⊗ℂ^D1))
    ρ = SparseUMPO(O)

    ψ′, _, _, _ = multiply(ψ₀, ρ, ψ)

    λ1, _, _ = leftFixedPoint(environment(ψ.AL, [O;;], adjoint.(ψ′.AL)))
    λ2, _, _ = leftFixedPoint(environment(ψ.AL, [O; O;;], adjoint.(ψ.AL)))
    @test abs(λ1[1] - sqrt(λ2[1])) < tol_low
end
