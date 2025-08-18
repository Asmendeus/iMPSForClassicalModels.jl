module iMPSForClassicalModels

using Reexport
@reexport using KrylovKit, TensorKit, TensorKit.TensorOperations
@reexport import Base: +, -, *, /, ==, iterate, promote_rule, convert, length, size, show, getindex, setindex!, lastindex, keys, similar, merge, merge!, iterate, complex
@reexport import TensorKit: ×, one, zero, dim, inner, scalar, space, domain, codomain, eltype, scalartype, numin, numout, numind, leftorth, rightorth, leftnull, rightnull, tsvd, adjoint, AdjointTensorMap, normalize!, normalize, norm, axpy!, axpby!, add!, add!!, dot, mul!, rmul!, NoTruncation, fuse, zerovector!, zerovector, scale, scale!, scale!!, sqrt, inv
@reexport import LinearAlgebra: BLAS, rank, qr, diag, I, diagm, ishermitian

# utils
export sign_first_element, SimpleIteration
# export EigenAlgorithm, GradientAlgorithm
export trivial, istrivial
export SimpleIterationInfo, LanczosInfo, ArnoldiInfo
export BondInfo, FixedPointInfo

# Tensor wrapper
export AbstractTensorWrapper
export AbstractLocalTensor, AbstractBondTensor, AbstractMPSTensor, AbstractMPOTensor
export LocalTensor, BondTensor, MPSTensor, MPOTensor
export AdjointLocalTensor, AdjointBondTensor, AdjointMPSTensor, AdjointMPOTensor
export AbstractEnvironmentTensor, LeftEnvironmentTensor, RightEnvironmentTensor
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
export getAL, getAR

# iMPS
export AbstractInfiniteMPS, DenseInfiniteMPS, coef, Center, iscanonical, isuniform, getAllCanonicalFormTensors
export AdjointInfiniteMPS
export canonicalize!, canonicalize
export setCenter!, setCenter
export uniformize!, uniformize
export InfiniteMPS, iMPS, randInfiniteMPS

# iMPO
export InfiniteMPO, iMPO, identityInfiniteMPO

# SparseInfiniteMPO
export SparseInfiniteMPO

# Impurity tensor
export AbstractLocalImpurity, LocalImpurity, expectation

# Algorithm

include("utils/SimpleIteration.jl")
include("utils/iterate.jl")
include("utils/Defaults.jl")
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
include("Method/getAL.jl")
include("Method/getAR.jl")

include("MPS/AbstractInfiniteMPS.jl")
include("MPS/InfiniteMPS.jl")
include("MPS/AdjointInfiniteMPS.jl")
include("MPS/canonicalize.jl")
include("MPS/setCenter.jl")
include("MPS/uniformize.jl")

include("MPO/InfiniteMPO.jl")

include("SparseMPO/SparseInfiniteMPO.jl")

include("LocalImpurity/LocalImpurity.jl")
include("LocalImpurity/expectation.jl")


end # module iMPSForClassicalModels
