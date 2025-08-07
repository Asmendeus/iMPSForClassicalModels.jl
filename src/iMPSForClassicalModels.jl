module iMPSForClassicalModels

using Reexport
@reexport using TensorKit, TensorKit.TensorOperations, OptimKit, KrylovKit
@reexport import Base: +, -, *, /, ==, promote_rule, convert, length, size, show, getindex, setindex!, lastindex, keys, similar, merge, merge!, iterate, complex
@reexport import TensorKit: ×, one, zero, dim, inner, scalar, space, domain, codomain, eltype, scalartype, leftorth, rightorth, leftnull, rightnull, tsvd, adjoint, normalize!, norm, axpy!, axpby!, add!, add!!, dot, mul!, rmul!, NoTruncation, fuse, zerovector!, zerovector, scale, scale!, scale!!
@reexport import LinearAlgebra: BLAS, rank, qr, diag, I, diagm

# utils
include("utils/Defaults.jl")

# Tensor wrapper
export AbstractTensorWrapper
include("TensorWrapper/TensorWrapper.jl")
export AbstractMPSTensor, MPSTensor, AdjointMPSTensor
include("TensorWrapper/MPSTensor.jl")
export LeftIsometricTensor, AdjointLeftIsometricTensor
include("TensorWrapper/LeftIsometricTensor.jl")
export RightIsometricTensor, AdjointRightIsometricTensor
include("TensorWrapper/RightIsometricTensor.jl")
export AbstractEnvironmentTensor, LeftEnvironmentTensor, RightEnvironmentTensor
include("TensorWrapper/EnvironmentTensor.jl")

# Transfer matrix

# Environment

# Impurity tensor

# MPS

# MPO

# Algorithm

end # module iMPSForClassicalModel
