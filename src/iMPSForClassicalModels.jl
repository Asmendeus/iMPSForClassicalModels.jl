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
export AbstractLocalTensor, AbstractBondTensor, AbstractMPSTensor, AbstractMPOTensor, AbstractMPSOrMPOTensor
export LocalTensor, BondTensor, MPSTensor, MPOTensor, MPSOrMPOTensor
export AdjointLocalTensor, AdjointBondTensor, AdjointMPSTensor, AdjointMPOTensor, AdjointMPSOrMPOTensor
include("TensorWrapper/LocalTensor.jl")
export LeftIsometricTensor, LeftIsometricBondTensor, LeftIsometricMPSTensor, LeftIsometricMPOTensor, LeftIsometricMPSOrMPOTensor
export AdjointLeftIsometricTensor, AdjointLeftIsometricBondTensor, AdjointLeftIsometricMPSTensor, AdjointLeftIsometricMPOTensor, AdjointLeftIsometricMPSOrMPOTensor
include("TensorWrapper/LeftIsometricTensor.jl")
export RightIsometricTensor, RightIsometricBondTensor, RightIsometricMPSTensor, RightIsometricMPOTensor, RightIsometricMPSOrMPOTensor
export AdjointRightIsometricTensor, AdjointRightIsometricBondTensor, AdjointRightIsometricMPSTensor, AdjointRightIsometricMPOTensor, AdjointRightIsometricMPSOrMPOTensor
include("TensorWrapper/RightIsometricTensor.jl")
export AbstractEnvironmentTensor, LeftEnvironmentTensor, RightEnvironmentTensor
include("TensorWrapper/EnvironmentTensor.jl")

export isAdjoint, isLeftIsometric, isRightIsometric, leftorth, rightorth
include("TensorWrapper/utils.jl")

# Transfer matrix
export AbstractTransferMatrix, AbstractMPSTransferMatrix, AbstractMPOTransferMatrix, AbstractMPSOrMPOTransferMatrix
include("TransferMatrix/AbstractTransferMatrix.jl")
export TransferMatrix, MPSTransferMatrix, MPOTransferMatrix, MPSOrMPOTransferMatrix
include("TransferMatrix/TransferMatrix.jl")

# Environment
export AbstractEnvironment
include("Environment/AbstractEnvironment.jl")
export BondEnvironment
include("Environment/BondEnvironment.jl")
export CenterEnvironment
include("Environment/CenterEnvironment.jl")
export MultirowCenterEnvironment
include("Environment/MultirowCenterEnvironment.jl")
export LeftEnvironment
include("Environment/LeftEnvironment.jl")

export Environment
include("Environment/utils.jl")

# MPS

# MPO

# Impurity tensor

# Method

# Algorithm

end # module iMPSForClassicalModel
