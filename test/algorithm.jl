using Test, iMPSForClassicalModels

D1 = 10
D2 = 200
d = 4

tol = 1e-12

@testset "approximate(iMPS{1, ComplexF64})" begin
    ψ = randCMPS(ComplexF64, 1, d, D2)

    ψ1₀ = randUMPS(ComplexF64, 1, d, D1)
    ψ1′, info = approximate(ψ1₀, ψ)
    @test abs(overlap(ψ, ψ1′)) ≈ 1
    @show overlap(ψ, ψ1′)
    if !info.converged
        @show info.normres
    end

    ψ2₀ = randUMPS(ComplexF64, 1, d, D2)
    ψ2′, info = approximate(ψ2₀, ψ)
    @test abs(overlap(ψ, ψ2′) - 1) < tol
    if !info.converged
        @show info.normres
    end
end

@testset "approximate(iMPS{2, ComplexF64})" begin
    ψ = randCMPS(ComplexF64, 2, d, D2)

    ψ1₀ = randUMPS(ComplexF64, 2, d, D1)
    ψ1′, info = approximate(ψ1₀, ψ)
    @test abs(overlap(ψ, ψ1′)) ≈ 1
    @show overlap(ψ, ψ1′)
    if !info.converged
        @show info.normres
    end

    ψ2₀ = randUMPS(ComplexF64, 2, d, D2)
    ψ2′, info = approximate(ψ2₀, ψ)
    @test abs(overlap(ψ, ψ2′) - 1) < tol
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
    @test abs(overlap(ρ, ρ1′)) ≈ 1
    @show overlap(ρ, ρ1′)
    if !info.converged
        @show info.normres
    end

    ρ2₀ = UMPO(O2′)
    ρ2′, info = approximate(ρ2₀, ρ)
    @test abs(overlap(ρ, ρ2′) - 1) < tol
    if !info.converged
        @show info.normres
    end
end
