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
export isAdjoint, isLeftIsometric, isRightIsometric, leftorth, rightorth

# Environment
export AbstractEnvironment
export TransferMatrix, MPSTransferMatrix, MPOTransferMatrix
export ChannelEnvironment, MPSChannelEnvironment, MPOChannelEnvironment
export BondEnvironment, CenterEnvironment
export environment

# Method
export SimpleIteration
export pushleft, pushright
export leftFixedPoint, rightFixedPoint, sign_first_element

# MPS

# MPO

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
include("Method/leftFixedPoint.jl")
include("Method/rightFixedPoint.jl")

end # module iMPSForClassicalModels
