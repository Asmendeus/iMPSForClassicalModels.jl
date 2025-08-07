"""
    abstract type AbstractEnvironmentTensor <: AbstractTensorWrapper

The right or left variational environment tensors.
"""
abstract type AbstractEnvironmentTensor <: AbstractTensorWrapper end

"""
    struct LeftEnvironmentTensor{N} <: AbstractEnvironmentTensor
        A::AbstractTensorMap
    end

Wrapper type for N-leg left environment tensor.

Convention (' marks codomain):

     __
    |  | —— 2
    |  |
    |  | —— 3
    |FL|    ⋮
    |  | —— N
    |  |
    |  | —— 1'
     ‾‾

# Constructors
    LeftEnvironmentTensor(::AbstractTensorMap)
    LeftEnvironmentTensor{N}(::AbstractTensorMap)
"""
struct LeftEnvironmentTensor{N} <: AbstractEnvironmentTensor
    A::AbstractTensorMap

    function LeftEnvironmentTensor(A::AbstractTensorMap)
        (length(codomain(A)) == 1 && length(domain(A)) > 0) || throw(ArgumentError("The codomain $(codomain(A)) does not conform to `LeftEnvironmentTensor` codomain"))
        N = length(domain(A)) + 1
        return new{N}(A)
    end
    function LeftEnvironmentTensor{N}(A::AbstractTensorMap) where N
        (length(codomain(A)) == 1 && length(domain(A)) == N-1 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `LeftEnvironmentTensor{$N}` space"))
        return new{N}(A)
    end
end

"""
    struct RightEnvironmentTensor{N} <: AbstractEnvironmentTensor
        A::AbstractTensorMap
    end

Wrapper type for N-leg right environment tensor.

Convention (' marks codomain):
               __
        1' —— |  |
              |  |
        2' —— |  |
        ⋮     |FR|
    (N-1)' —— |  |
              |  |
         N —— |  |
               ‾‾
# Constructors
    RightEnvironmentTensor(::AbstractTensorMap)
    RightEnvironmentTensor{N}(::AbstractTensorMap)
"""
struct RightEnvironmentTensor{N} <: AbstractEnvironmentTensor
    A::AbstractTensorMap

    function RightEnvironmentTensor(A::AbstractTensorMap)
        (length(domain(A)) == 1 && length(codomain(A) > 0)) || throw(ArgumentError("The codomain $(codomain(A)) does not conform to `RightEnvironmentTensor` codomain"))
        N = length(codomain(A)) + 1
        return new{N}(A)
    end
    function RightEnvironmentTensor{N}(A::AbstractTensorMap) where N
        (length(domain(A)) == 1 && length(codomain(A)) == N-1 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `RightEnvironmentTensor{$N}` space"))
        return new{N}(A)
    end
end
