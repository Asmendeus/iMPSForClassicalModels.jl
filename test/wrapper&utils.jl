using Test, iMPSForClassicalModels

tol = 1e-8

DL = 6  # left virtual bond dimension
DR = 8  # right virtual bond dimension

dp = 3    # physical bond dimension
d1 = 4    # first layer MPO tensor's auxiliary bond dimension
d2 = 5    # second layer  MPO tensor's auxiliary bond dimension

# MPS: -- A --
#         |
A = MPSTensor(TensorMap(rand, ℝ^DL⊗ℝ^dp, ℝ^DL))
# MPS: -- A1 -- A2 --
#         |     |
A1 = MPSTensor(TensorMap(rand, ℝ^DL⊗ℝ^dp, ℝ^DR))
A2 = MPSTensor(TensorMap(rand, ℝ^DR⊗ℝ^dp, ℝ^DL))

# MPO:    |
#      -- O --
#         |
O = MPOTensor(TensorMap(rand, ℝ^d1⊗ℝ^dp, ℝ^dp⊗ℝ^d1))
# MultirowMPO:    |     |
#              —— O1 —— O1 ——
#                 |     |
#              —— O2 —— O2 ——
#                 |     |
O1 = MPOTensor(TensorMap(rand, ℝ^d1⊗ℝ^dp, ℝ^dp⊗ℝ^d1))
O2 = MPOTensor(TensorMap(rand, ℝ^d2⊗ℝ^dp, ℝ^dp⊗ℝ^d2))

# AdjointMPS:    |
#             —— B ——
B = AdjointMPSTensor(TensorMap(rand, ℝ^DL, ℝ^DL⊗ℝ^dp))
# AdjointMPS:    |     |
#             —— B1 —— B2 ——
B1 = AdjointMPSTensor(TensorMap(rand, ℝ^DR, ℝ^DL⊗ℝ^dp))
B2 = AdjointMPSTensor(TensorMap(rand, ℝ^DL, ℝ^DR⊗ℝ^dp))

# BondTensor: -- C --
C = BondTensor(TensorMap(rand, ℝ^DL, ℝ^DR))

# single layer EnvironmentTensor:
#    ___             ___
#   |   | --     -- |   |
#   |   |           |   |
#   |FL1| ——     —— |FR1|
#   |   |           |   |
#   |   | --     -- |   |
#    ‾‾‾             ‾‾‾
FL1 = LeftEnvironmentTensor(TensorMap(rand, ℝ^DL, ℝ^d1⊗ℝ^DL))
FR1 = RightEnvironmentTensor(TensorMap(rand, ℝ^DR⊗ℝ^d1, ℝ^DR))
# bilayer EnvironmentTensor:
#    ___             ___
#   |   | --     -- |   |
#   |   |           |   |
#   |   | ——     —— |   |
#   |FL2|           |FR2|
#   |   | ——     —— |   |
#   |   |           |   |
#   |   | --     -- |   |
#    ‾‾‾             ‾‾‾
FL2 = LeftEnvironmentTensor(TensorMap(rand, ℝ^DL, ℝ^d2⊗ℝ^d1⊗ℝ^DL))
FR2 = RightEnvironmentTensor(TensorMap(rand, ℝ^DL⊗ℝ^d1⊗ℝ^d2, ℝ^DL))

# leftorth and rightorth
AL, RA, _ = leftorth(A)
LA, AR, _ = rightorth(A)
OL, RO, _ = leftorth(O)
LO, OR, _ = rightorth(O)
BL, RB, _ = leftorth(B)
LB, BR, _ = rightorth(B)
CL, RC, _ = leftorth(C)
LC, CR, _ = rightorth(C)

AL1, _, _ = leftorth(A1)
_, AR1, _ = rightorth(A1)
BL1, _, _ = leftorth(B1)
_, BR1, _ = rightorth(B1)

AL2, _, _ = leftorth(A2)
_, AR2, _ = rightorth(A2)
BL2, _, _ = leftorth(B2)
_, BR2, _ = rightorth(B2)

# TransferMatrix
trans_AB = TransferMatrix(A, B)
trans_ABL = TransferMatrix(A, BL)
trans_AsBs = TransferMatrix([A1, A2], [B1, B2])
trans_OO = TransferMatrix(O, O')

# Environment
bondEnv1 = Environment(FL1, FR1)
bondEnv2 = Environment(FL2, FR2)

centerEnv = Environment(FL1, O1, FR1)

leftEnv1 = Environment(AL, O, BL)
leftEnv2 = Environment([AL1, AL2], [O1, O1], [BL1, BL2])
leftEnv3 = Environment(OL, O, OL')

rightEnv1 = Environment(AR, O, BR)
rightEnv2 = Environment([AR1, AR2], [O1, O1], [BR1, BR2])
rightEnv3 = Environment(OR, O, OR')

multirowCenterEnv = Environment(FL2, [O1, O2], FR2)

multirowLeftEnv1 = Environment([AL1, AL2], [O1 O1; O2 O2], [BL1, BL2])
multirowLeftEnv2 = Environment(AL1, [O1; O2], BL1)
multirowLeftEnv3 = Environment(OL, [O1, O2], OL')

multirowRightEnv1 = Environment([AR1, AR2], [O1 O1; O2 O2], [BR1, BR2])
multirowRightEnv2 = Environment(AR1, [O1; O2], BR1)
multirowRightEnv3 = Environment(OR, [O1, O2], OR')

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

    @test isLeftIsometric(CL)
    @test norm(CL.A * RC.A - C.A) < tol

    @test isRightIsometric(CR)
    @test norm(LC.A * CR.A - C.A) < tol
end
