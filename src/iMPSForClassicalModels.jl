module iMPSForClassicalModels

using Reexport
@reexport using TensorKit, TensorKit.TensorOperations, OptimKit, KrylovKit
@reexport import Base: +, -, *, /, ==, promote_rule, convert, length, size, show, getindex, setindex!, lastindex, keys, similar, merge, merge!, iterate, complex
@reexport import TensorKit: ×, one, zero, dim, inner, scalar, space, domain, codomain, eltype, scalartype, numin, numout, numind, leftorth, rightorth, leftnull, rightnull, tsvd, adjoint, AdjointTensorMap, normalize!, norm, axpy!, axpby!, add!, add!!, dot, mul!, rmul!, NoTruncation, fuse, zerovector!, zerovector, scale, scale!, scale!!
@reexport import LinearAlgebra: BLAS, rank, qr, diag, I, diagm

# utils
include("utils/Defaults.jl")
include("utils/TensorMap.jl")
export trivial, istrivial
include("utils/trivial.jl")
export BondInfo
include("utils/Info.jl")

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

# Environment
export AbstractEnvironment
include("Environment/AbstractEnvironment.jl")
export AbstractTransferMatrix, AbstractMPSTransferMatrix, AbstractMPOTransferMatrix
export TransferMatrix, MPSTransferMatrix, MPOTransferMatrix
include("Environment/TransferMatrix.jl")
export IsometricEnvironment, IsometricMPSEnvironment, IsometricMPOEnvironment
include("Environment/IsometricEnvironment.jl")
export BondEnvironment, CenterEnvironment
include("Environment/BondEnvironment.jl")
include("Environment/CenterEnvironment.jl")

export environment
include("Environment/environment.jl")

# Method

# MPS

# MPO

# Impurity tensor

# Algorithm

end # module iMPSForClassicalModel
