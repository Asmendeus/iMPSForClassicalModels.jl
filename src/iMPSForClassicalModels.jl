module iMPSForClassicalModels

using Reexport
@reexport using TensorKit, TensorKit.TensorOperations, OptimKit, KrylovKit
@reexport import Base: +, -, *, /, ==, promote_rule, convert, length, size, show, getindex, setindex!, lastindex, keys, similar, merge, merge!, iterate, complex
@reexport import TensorKit: ×, one, zero, dim, inner, scalar, domain, codomain, eltype, scalartype, leftorth, rightorth, leftnull, rightnull, tsvd, adjoint, normalize!, norm, axpy!, axpby!, add!, add!!, dot, mul!, rmul!, NoTruncation, fuse, zerovector!, zerovector, scale, scale!, scale!!
@reexport import LinearAlgebra: BLAS, rank, qr, diag, I, diagm

# utils
include("utils/Defaults.jl")

# Tensor wrapper
export AbstractTensorWrapper
include("TensorWrapper/TensorWrapper.jl")
export AbstractMPSTensor, MPSTensor, LeftIsometricTensor, RightIsometricTensor, AdjointMPSTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor
include("TensorWrapper/MPSTensor.jl")
export AbstractMPSBondTensor, MPSBondTensor
include("TensorWrapper/MPSBondTensor.jl")
export AbstractMPOTensor, MPOTensor
include("TensorWrapper/MPOTensor.jl")
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
export LeftEnvironment, RightEnvironment, MidEnvironment
include("Environment/LeftEnvironment.jl")
include("Environment/RightEnvironment.jl")
include("Environment/MidEnvironment.jl")
export MultirowLeftEnvironment, MultirowRightEnvironment, MultirowMidEnvironment
include("Environment/MultirowLeftEnvironment.jl")
include("Environment/MultirowRightEnvironment.jl")
include("Environment/MultirowMidEnvironment.jl")

export Environment
include("Environment/utils.jl")

# MPS

# MPO

# Impurity tensor

# Algorithm

end # module iMPSForClassicalModel
