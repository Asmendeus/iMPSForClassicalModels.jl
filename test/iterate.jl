using Test, iMPSForClassicalModels
import iMPSForClassicalModels: lambda, division, minus, norm_max
using LinearAlgebra

tol = 1e-10

lambda(n::Number) = norm(n) * sign(n)
division(n1::Number, n2::Number) = n1 / n2
minus(n1::Number, n2::Number) = n1 - n2
norm_max(n::Number) = norm(n)
@testset "λx = x^2 - 2*x + 3" begin

    f = x -> x^2 - 2*x + 3      # x = f(x) / λ  => λ = 2, x = 1
    v = 2.0
    λ, v, _ = iterate(f, v)

    @test λ == 2
    @test v == 1
end

lambda(v::Vector{<:Number}) = norm(v) * sign(v[1])
division(v::Vector{<:Number}, n::Number) = v ./ n
minus(v1::Vector{<:Number}, v2::Vector{<:Number}) = v1 - v2
norm_max(v::Vector{<:Number}) = norm(v)
@testset "Av = λv" begin
    A = Hermitian(rand(4, 4))
    f = x -> A * x
    v = rand(4)
    λ, v, _ = iterate(f, v)

    λm, vm = eigsolve(A, rand(4), 1, :LM)
    @test abs(λm[1] - λ) < tol
    @test norm(vm[1] - v / sign(v[1])) < tol
end

lambda(v::Vector{<:Matrix{<:Number}}) = norm.(v) .* map(x->sign(x[1]), v)
division(v::Vector{<:Matrix{<:Number}}, n::Vector{<:Number}) = v ./ n
minus(v1::Vector{<:Matrix{<:Number}}, v2::Vector{<:Matrix{<:Number}}) = v1 .- v2
norm_max(v::Vector{<:Matrix{<:Number}}) = maximum(norm.(v))
@testset "ABCv = ∏λ v & vABC = ∏λ v" begin
    M = [Hermitian(rand(4, 4)) for _ in 1:3]

    f₊ = [x -> x * M[i] for i in 1:3]
    f₋ = [x -> M[i] * x for i in 1:3]

    vL = [rand(1, 4) for _ in 1:3]
    vR = [rand(4, 1) for _ in 1:3]

    λL, vL, _ = iterate(f₊, vL, true)
    λR, vR, _ = iterate(f₋, vR, false)

    H123 = M[1] * M[2] * M[3]
    H231 = M[2] * M[3] * M[1]
    H312 = M[3] * M[1] * M[2]
    H321 = transpose(M[3]) * transpose(M[2]) * transpose(M[1])
    H132 = transpose(M[1]) * transpose(M[3]) * transpose(M[2])
    H213 = transpose(M[2]) * transpose(M[1]) * transpose(M[3])

    val₊1, vecs123 = eigsolve(H123, rand(4), 1, :LM)
    val₊2, vecs231 = eigsolve(H231, rand(4), 1, :LM)
    val₊3, vecs312 = eigsolve(H312, rand(4), 1, :LM)
    val₋1, vecs321 = eigsolve(H321, rand(4), 1, :LM)
    val₋2, vecs132 = eigsolve(H132, rand(4), 1, :LM)
    val₋3, vecs213 = eigsolve(H213, rand(4), 1, :LM)

    @test abs(prod(λL) - val₊1[1]) < tol && abs(prod(λL) - val₊2[1]) < tol && abs(prod(λL) - val₊3[1]) < tol
    @test norm(vL[1][:] - vecs321[1][:, 1] / sign(vecs321[1][:, 1][1])) < tol
    @test norm(vL[2][:] - vecs132[1][:, 1] / sign(vecs132[1][:, 1][1])) < tol
    @test norm(vL[3][:] - vecs213[1][:, 1] / sign(vecs213[1][:, 1][1])) < tol
    @test abs(prod(λR) - val₋1[1]) < tol && abs(prod(λR) - val₋2[1]) < tol && abs(prod(λR) - val₋3[1]) < tol
    @test norm(vR[1][:] - vecs231[1][:, 1] / sign(vecs231[1][:, 1][1])) < tol
    @test norm(vR[2][:] - vecs312[1][:, 1] / sign(vecs312[1][:, 1][1])) < tol
    @test norm(vR[3][:] - vecs123[1][:, 1] / sign(vecs123[1][:, 1][1])) < tol
end
