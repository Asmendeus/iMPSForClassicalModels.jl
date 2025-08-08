using Test, TensorKit
using iMPSForClassicalModels

tol = 1e-8

DL = 6
d = 3
DR = 5

A = MPSTensor(TensorMap(rand, ℝ^DL⊗ℝ^d, ℝ^DR))           # MPS tensor
O = MPSTensor(TensorMap(rand, ℝ^d⊗ℝ^d, ℝ^d⊗ℝ^d))        # MPO tensor
B = AdjointMPSTensor(TensorMap(rand, ℝ^DR, ℝ^DL⊗ℝ^d))    # adjoint MPS tensor
C = MPSTensor(TensorMap(rand, ℝ^DL, ℝ^DR))                # bond tensor

AL, RA, _ = leftorth(A)
LA, AR, _ = rightorth(A)
OL, RO, _ = leftorth(O)
LO, OR, _ = rightorth(O)
BL, RB, _ = leftorth(B)
LB, BR, _ = rightorth(B)

trans_AB = TransferMatrix(A, B)

env = Environment(A, O, B)

@testset "leftorth & rightorth" begin
    @test isLeftIsometric(AL)
    @test norm(AL.A * RA.A - A.A) < tol

    @test isRightIsometric(AR)
    @test norm(LA.A * permute(AR.A, (1,), (2, 3)) - permute(A.A, (1,), (2, 3))) < tol

    @test isLeftIsometric(OL)
    @test norm(LO.A * permute(OR.A, (1,), (2, 3, 4)) - permute(O.A, (1,), (2, 3, 4))) < tol

    @test isRightIsometric(OR)
    @test norm(permute(OL.A, (1, 2, 3), (4,)) * RO.A - permute(O.A, (1, 2, 3), (4,))) < tol

    @test isLeftIsometric(BL)
    @test norm(BL.A' * RB.A' - B.A') < tol

    @test isRightIsometric(BR)
    @test norm(LB.A' * permute(BR.A', (1,), (2, 3)) - permute(B.A', (1,), (2, 3))) < tol
end