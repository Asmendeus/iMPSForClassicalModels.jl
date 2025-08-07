using Test, TensorKit
using iMPSForClassicalModels

tol = 1e-8

DL = 6
d = 3
DR = 5

B = AdjointMPSTensor(TensorMap(rand, ℝ^DR, ℝ^DL⊗ℝ^d))
O = MPOTensor(TensorMap(rand, ℝ^d⊗ℝ^d, ℝ^d⊗ℝ^d))
A = MPSTensor(TensorMap(rand, ℝ^DL⊗ℝ^d, ℝ^DR))

AL, RA = leftorth(A)
LA, AR = rightorth(A)
BL, RB = leftorth(B)
LB, BR = rightorth(B)

env = Environment(A, O, B)

@testset "leftorth & rightorth" begin
    @test isLeftIsometric(AL)
    @test norm(AL.A * RA - A.A) < tol

    @test isRightIsometric(AR)
    @test norm(LA * permute(AR.A, (1,), (2, 3)) - permute(A.A, (1,), (2, 3))) < tol

    @test isLeftIsometric(BL)
    @test norm(permute(RB, (2,), (1,)) * BL.A - B.A) < tol

    @test isRightIsometric(BR)
    @test norm(LB * permute(BR.A, (2,), (3, 1)) - permute(B.A, (2,), (3, 1))) < tol
end
