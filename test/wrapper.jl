using Test, iMPSForClassicalModels

tol = 1e-8

D = 6   # bond dimension
DL = 6  # left bond dimension
DR = 8  # right bond dimension

dp = 3    # physical dimension
d1 = 4    # first layer MPO tensor's auxiliary bond dimension
d2 = 5    # second layer  MPO tensor's auxiliary bond dimension

# MPS: -- A --
#         |
A = MPSTensor(TensorMap(rand, ℝ^D⊗ℝ^dp, ℝ^D))
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

# TransferMatrix
trans_AB = environment(A, B)
trans_AB = environment([A1, A2], [A1', A2'])
trans_OO = environment(O, O')
trans_OO = environment([O, O], [O', O'])

# Environment
channelEnv1 = environment(AL, O, BL)
channelEnv2 = environment(OR, O, OR')
channelEnv3 = environment([A1, A2], [O1 O1; O2 O2], [A1', A2'])

bondEnv1 = environment(FL1, FR1)
bondEnv2 = environment(FL2, FR2)

centerEnv1 = environment(FL1, O1, FR1)
centerEnv2 = environment(FL1, [O1, O2], FR1)

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

ψ1 = iMPS(2)
ψ2 = randInfiniteMPS(ComplexF64, [ℂ^4, ℂ^3], [ℂ^10, ℂ^10])
ψ3 = randInfiniteMPS(Float64, 1, 4, 8)'

ρ1 = iMPO(2)
ρ2 = identityInfiniteMPO(ComplexF64, [ℝ^4, ℝ^3])
ρ3 = identityInfiniteMPO(Float64, 3, ℂ^4)'

Z = SparseMPO(MPOTensor[O1 O1; O2 O2])
M = LocalImpurity(MPOTensor[O1 O1; O2 O2])
