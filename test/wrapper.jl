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
A = MPSTensor(TensorMap(rand, ℂ^D⊗ℂ^dp, ℂ^D))
# MPS: -- A1 -- A2 --
#         |     |
A1 = MPSTensor(TensorMap(rand, ℂ^DL⊗ℂ^dp, ℂ^DR))
A2 = MPSTensor(TensorMap(rand, ℂ^DR⊗ℂ^dp, ℂ^DL))

# MPO:    |
#      -- O --
#         |
O = MPOTensor(TensorMap(rand, ℂ^d1⊗ℂ^dp, ℂ^dp⊗ℂ^d1))

# MultirowMPO:    |     |
#              —— O1 —— O1 ——
#                 |     |
#              —— O2 —— O2 ——
#                 |     |
O1 = MPOTensor(TensorMap(rand, ℂ^d1⊗ℂ^dp, ℂ^dp⊗ℂ^d1))
O2 = MPOTensor(TensorMap(rand, ℂ^d2⊗ℂ^dp, ℂ^dp⊗ℂ^d2))

# AdjointMPS:    |
#             —— B ——
B = AdjointMPSTensor(TensorMap(rand, ℂ^DL, ℂ^DL⊗ℂ^dp))

# BondTensor: -- C --
C = BondTensor(TensorMap(rand, ℂ^DL, ℂ^DL))

# single layer EnvironmentTensor:
#    ___             ___
#   |   | --     -- |   |
#   |   |           |   |
#   |FL1| ——     —— |FR1|
#   |   |           |   |
#   |   | --     -- |   |
#    ‾‾‾             ‾‾‾
FL1 = LeftEnvironmentTensor(TensorMap(rand, ℂ^DL, ℂ^d1⊗ℂ^DL))
FR1 = RightEnvironmentTensor(TensorMap(rand, ℂ^DR⊗ℂ^d1, ℂ^DR))
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
FL2 = LeftEnvironmentTensor(TensorMap(rand, ℂ^DL, ℂ^d2⊗ℂ^d1⊗ℂ^DL))
FR2 = RightEnvironmentTensor(TensorMap(rand, ℂ^DL⊗ℂ^d1⊗ℂ^d2, ℂ^DL))

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

# MPS
ψ1 = iMPS(A)
ψ2 = randInfiniteMPS(ComplexF64, 2, 4, 10)

# MPO
ρ1 = iMPO(O)
ρ2 = identityInfiniteMPO(ComplexF64, 2, 4)

# SparseUMPO
Z1 = SparseGeneralMPO(O)
Z2 = SparseGeneralMPO([O1 O1; O2 O2])

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

@testset "iMPS (L = 1)" begin
    @test abs(overlap(ψ1, ψ1) - 1) < tol
    @test norm(ψ1.AL[1] * ψ1.C[1] - ψ1.C[1] * ψ1.AR[1]) < tol
end
@testset "iMPS (L = 2)" begin
    @test abs(overlap(ψ2, ψ2) - 1) < tol
    @test norm(ψ2.AL[1] * ψ2.C[1] - ψ2.C[2] * ψ2.AR[1]) < tol
    @test norm(ψ2.AL[2] * ψ2.C[2] - ψ2.C[1] * ψ2.AR[2]) < tol
end
@testset "iMPO (L = 1)" begin
    @test abs(overlap(ρ1, ρ1) - 1) < tol
    @test norm(ρ1.AL[1] * ρ1.C[1] - ρ1.C[1] * ρ1.AR[1]) < tol
end
