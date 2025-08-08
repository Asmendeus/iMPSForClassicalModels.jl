module iMPSForClassicalModels

using Reexport
@reexport using TensorKit, TensorKit.TensorOperations, OptimKit, KrylovKit
@reexport import Base: +, -, *, /, ==, promote_rule, convert, length, size, show, getindex, setindex!, lastindex, keys, similar, merge, merge!, iterate, complex
@reexport import TensorKit: ×, one, zero, dim, inner, scalar, space, domain, codomain, eltype, scalartype, numin, numout, numind, leftorth, rightorth, leftnull, rightnull, tsvd, adjoint, normalize!, norm, axpy!, axpby!, add!, add!!, dot, mul!, rmul!, NoTruncation, fuse, zerovector!, zerovector, scale, scale!, scale!!
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
export LeftIsometricTensor, LeftIsometricMPSTensor, LeftIsometricMPOTensor, LeftIsometricMPSOrMPOTensor
export AdjointLeftIsometricTensor, AdjointLeftIsometricMPSTensor, AdjointLeftIsometricMPOTensor, AdjointLeftIsometricMPSOrMPOTensor
include("TensorWrapper/LeftIsometricTensor.jl")
export RightIsometricTensor, RightIsometricMPSTensor, RightIsometricMPOTensor, RightIsometricMPSOrMPOTensor
export AdjointRightIsometricTensor, AdjointRightIsometricMPSTensor, AdjointRightIsometricMPOTensor, AdjointRightIsometricMPSOrMPOTensor
include("TensorWrapper/RightIsometricTensor.jl")
export AbstractEnvironmentTensor, LeftEnvironmentTensor, RightEnvironmentTensor
include("TensorWrapper/EnvironmentTensor.jl")

export isLeftIsometric, isRightIsometric, leftorth, rightorth
include("TensorWrapper/utils.jl")

# Transfer matrix
export AbstractTransferMatrix
include("TransferMatrix/AbstractTransferMatrix.jl")
export TransferMatrix
include("TransferMatrix/TransferMatrix.jl")

# Environment
export AbstractEnvironment
include("Environment/AbstractEnvironment.jl")
export BondEnvironment
include("Environment/BondEnvironment.jl")


# Impurity tensor

# MPS

# MPO

# SparseMPO

# Algorithm

end # module iMPSForClassicalModel
