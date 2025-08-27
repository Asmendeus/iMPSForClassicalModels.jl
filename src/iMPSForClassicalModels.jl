module iMPSForClassicalModels

using Reexport
@reexport using KrylovKit, OptimKit, TensorKit, TensorKit.TensorOperations
@reexport import Base: +, -, *, /, ==, iterate, promote_rule, convert, length, size, show, getindex, setindex!, lastindex, keys, similar, merge, merge!, iterate, complex, sign
@reexport import TensorKit: ×, tr, one, zero, dim, inner, scalar, space, domain, codomain, eltype, scalartype, numin, numout, numind, leftorth, rightorth, leftnull, rightnull, tsvd, adjoint, AdjointTensorMap, norm, normalize!, normalize, axpy!, axpby!, add!, add!!, dot, mul!, rmul!, NoTruncation, fuse, zerovector!, zerovector, scale, scale!, scale!!, sqrt, inv
@reexport import LinearAlgebra: BLAS, rank, qr, diag, I, diagm, ishermitian

# utils
export SimpleIterator
# export EigenAlgorithm, GradientAlgorithm
export trivial, istrivial
export SimpleIteratorInfo, LanczosInfo, ArnoldiInfo, BondInfo

# Tensor wrapper
export AbstractTensorWrapper
export lambda, division, minus, norm_max
export AbstractLocalTensor, AbstractBondTensor, AbstractMPSTensor, AbstractMPOTensor
export LocalTensor, BondTensor, MPSTensor, MPOTensor
export AdjointLocalTensor, AdjointBondTensor, AdjointMPSTensor, AdjointMPOTensor
export AbstractEnvironmentTensor, LeftEnvironmentTensor, RightEnvironmentTensor
export leftVirtualSpace, rightVirtualSpace, physicalSpace, extraPhysicalSpace
export isadjoint, isLeftIsometric, isRightIsometric, leftorth, rightorth

# Environment
export AbstractEnvironment
export TransferMatrix, MPSTransferMatrix, MPOTransferMatrix
export ChannelEnvironment, MPSChannelEnvironment, MPOChannelEnvironment
export BondEnvironment, CenterEnvironment
export environment, contract

# Method
export pushleft, pushright, pushmid
export leftFixedPoint, rightFixedPoint, midFixedPoint

# iMPS
export AbstractGeneralInfiniteMPS, AbstractInfiniteMPS, DenseInfiniteMPS, getAL, getAR, getAR, getC, entropy
export InfiniteMPS, iMPS, randInfiniteMPS
export canonicalize
export AdjointInfiniteMPS
export overlap

# iMPO
export InfiniteMPO, iMPO, identityInfiniteMPO

# SparseMPO
export SparseGeneralMPO, getA

# Impurity tensor
export AbstractLocalImpurity, LocalImpurity, expectation

# Algorithm
export approximate, multiply
export VUMPS, ViTEBD, fixedBondary

include("utils/SimpleIterator.jl")
include("utils/iterate.jl")
include("utils/TensorMap.jl")
include("utils/trivial.jl")
include("utils/Info.jl")

include("TensorWrapper/TensorWrapper.jl")
include("TensorWrapper/LocalTensor.jl")
include("TensorWrapper/AdjointLocalTensor.jl")
include("TensorWrapper/EnvironmentTensor.jl")
include("TensorWrapper/utils.jl")

include("Environment/AbstractEnvironment.jl")
include("Environment/TransferMatrix.jl")
include("Environment/ChannelEnvironment.jl")
include("Environment/BondEnvironment.jl")
include("Environment/CenterEnvironment.jl")
include("Environment/environment.jl")
include("Environment/contract.jl")

include("Method/pushleft.jl")
include("Method/pushright.jl")
include("Method/pushmid.jl")
include("Method/leftFixedPoint.jl")
include("Method/rightFixedPoint.jl")
include("Method/midFixedPoint.jl")

include("MPS/AbstractInfiniteMPS.jl")
include("MPS/InfiniteMPS.jl")
include("MPS/canonicalize.jl")
include("MPS/AdjointInfiniteMPS.jl")
include("MPS/overlap.jl")

include("MPO/InfiniteMPO.jl")

include("SparseMPO/SparseGeneralMPO.jl")

include("LocalImpurity/LocalImpurity.jl")
include("LocalImpurity/expectation.jl")

include("Algorithm/getAL.jl")
include("Algorithm/getAR.jl")
include("Algorithm/approximate.jl")
include("Algorithm/multiply.jl")
include("Algorithm/fixedBondary/algorithm.jl")
include("Algorithm/fixedBondary/fixedBondary.jl")

include("utils/Defaults.jl")

end # module iMPSForClassicalModels
