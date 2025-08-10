module iMPSForClassicalModels

using Reexport
@reexport using KrylovKit, TensorKit, TensorKit.TensorOperations, OptimKit
@reexport import Base: +, -, *, /, ==, iterate, promote_rule, convert, length, size, show, getindex, setindex!, lastindex, keys, similar, merge, merge!, iterate, complex
@reexport import TensorKit: ×, one, zero, dim, inner, scalar, space, domain, codomain, eltype, scalartype, numin, numout, numind, leftorth, rightorth, leftnull, rightnull, tsvd, adjoint, AdjointTensorMap, normalize!, norm, axpy!, axpby!, add!, add!!, dot, mul!, rmul!, NoTruncation, fuse, zerovector!, zerovector, scale, scale!, scale!!
@reexport import LinearAlgebra: BLAS, rank, qr, diag, I, diagm

# utils
include("utils/Defaults.jl")
include("utils/TensorMap.jl")
export trivial, istrivial
include("utils/trivial.jl")
export SimpleIterationInfo, LanczosInfo, ArnoldiInfo
export BondInfo, FixedPointInfo
include("utils/Info.jl")
export SimpleIteration
include("utils/SimpleIteration.jl")
include("utils/iterate.jl")

# Tensor wrapper
export AbstractTensorWrapper
include("TensorWrapper/TensorWrapper.jl")
export AbstractLocalTensor, AbstractBondTensor, AbstractMPSTensor, AbstractMPOTensor
export LocalTensor, BondTensor, MPSTensor, MPOTensor
export AdjointLocalTensor, AdjointBondTensor, AdjointMPSTensor, AdjointMPOTensor
include("TensorWrapper/LocalTensor.jl")
export AbstractEnvironmentTensor, LeftEnvironmentTensor, RightEnvironmentTensor
include("TensorWrapper/EnvironmentTensor.jl")

export isAdjoint, isLeftIsometric, isRightIsometric, leftorth, rightorth
include("TensorWrapper/utils.jl")

# Transfer matrix
export AbstractTransferMatrix, AbstractMPSTransferMatrix, AbstractMPOTransferMatrix
include("TransferMatrix/AbstractTransferMatrix.jl")
export TransferMatrix, MPSTransferMatrix, MPOTransferMatrix
include("TransferMatrix/TransferMatrix.jl")

# Environment
export AbstractEnvironment
include("Environment/AbstractEnvironment.jl")
export IsometricEnvironment, IsometricMPSEnvironment, IsometricMPOEnvironment
include("Environment/IsometricEnvironment.jl")
export BondEnvironment, CenterEnvironment
include("Environment/BondEnvironment.jl")
include("Environment/CenterEnvironment.jl")

export environment
include("Environment/environment.jl")

# Method
export pushleft, pushright
include("Method/push.jl")

# MPS

# MPO

# Impurity tensor

# Algorithm

end # module iMPSForClassicalModel
