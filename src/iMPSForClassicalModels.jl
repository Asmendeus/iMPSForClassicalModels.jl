module iMPSForClassicalModels

using Reexport
@reexport using KrylovKit, TensorKit, TensorKit.TensorOperations
@reexport import Base: +, -, *, /, ==, iterate, promote_rule, convert, length, size, show, getindex, setindex!, lastindex, keys, similar, merge, merge!, iterate, complex
@reexport import TensorKit: ×, one, zero, dim, inner, scalar, space, domain, codomain, eltype, scalartype, numin, numout, numind, leftorth, rightorth, leftnull, rightnull, tsvd, adjoint, AdjointTensorMap, normalize!, normalize, norm, axpy!, axpby!, add!, add!!, dot, mul!, rmul!, NoTruncation, fuse, zerovector!, zerovector, scale, scale!, scale!!
@reexport import LinearAlgebra: BLAS, rank, qr, diag, I, diagm

# utils
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
export environment

# Method
export sign_first_element, SimpleIteration
export pushleft, pushright, pushmid
export leftFixedPoint, rightFixedPoint, midFixedPoint

# iMPS & iMPO
export AbstractInfiniteMPS, DenseInfiniteMPS, issparse, coef, Center
export InfiniteMPS, iMPS

export InfiniteMPO, iMPO

# SparseMPO

# Impurity tensor

# Algorithm

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

include("Method/SimpleIteration.jl")
include("Method/iterate.jl")
include("Method/pushleft.jl")
include("Method/pushright.jl")
include("Method/pushmid.jl")
include("Method/leftFixedPoint.jl")
include("Method/rightFixedPoint.jl")
include("Method/midFixedPoint.jl")

include("MPS/AbstractInfiniteMPS.jl")
include("MPS/InfiniteMPS.jl")

include("MPO/InfiniteMPO.jl")

end # module iMPSForClassicalModels
