"""
    abstract type AbstractEnvironmentTensor <: AbstractTensorWrapper

The right or left variational environment tensors.
"""
abstract type AbstractEnvironmentTensor{N} <: AbstractTensorWrapper end

"""
    struct LeftEnvironmentTensor{N} <: AbstractEnvironmentTensor{N}
        A::AbstractTensorMap
    end

Wrapper type for rank-N variational left environment tensor.

Convention (' marks codomain):

     __
    |  | —— N'
    |  |
    |  | —— (N-1)
    |FL|    ⋮
    |  | —— 2
    |  |
    |  | —— 1
     ‾‾

# Constructors
    LeftEnvironmentTensor{N}(::AbstractTensorMap)
    LeftEnvironmentTensor(::AbstractTensorMap)
"""
struct LeftEnvironmentTensor{N} <: AbstractEnvironmentTensor{N}
    A::AbstractTensorMap

    function LeftEnvironmentTensor{N}(A::AbstractTensorMap) where N
        (length(codomain(A)) == 1 && length(domain(A)) == N-1) || throw(ArgumentError("The space $(space(A)) does not conform to `LeftEnvironmentTensor` space"))
        return new{N}(A)
    end
    function LeftEnvironmentTensor(A::AbstractTensorMap)
        length(codomain(A)) == 1 || throw(ArgumentError("The codomain $(codomain(A)) does not conform to `LeftEnvironmentTensor` codomain"))
        N = length(domain(A)) + 1
        return new{N}(A)
    end
end

"""
    struct RightEnvironmentTensor{N} <: AbstractEnvironmentTensor{N}
        A::AbstractTensorMap
    end

Wrapper type for rank-N variational right environment tensor.

Convention (' marks codomain):
               __
        N  —— |  |
              |  |
    (N-1)' —— |  |
        ⋮     |FR|
        2' —— |  |
              |  |
        1' —— |  |
               ‾‾
# Constructors
    RightEnvironmentTensor{N}(::AbstractTensorMap)
    RightEnvironmentTensor(::AbstractTensorMap)
"""
struct RightEnvironmentTensor{N} <: AbstractEnvironmentTensor{N}
    A::AbstractTensorMap

    function RightEnvironmentTensor{N}(A::AbstractTensorMap) where N
        (length(domain(A)) == 1 && length(codomain(A)) == N-1) || throw(ArgumentError("The space $(space(A)) does not conform to `RightEnvironmentTensor` space"))
        return new{N}(A)
    end
    function RightEnvironmentTensor(A::AbstractTensorMap)
        length(domain(A)) == 1 || throw(ArgumentError("The codomain $(codomain(A)) does not conform to `RightEnvironmentTensor` codomain"))
        N = length(codomain(A)) + 1
        return new{N}(A)
    end
end
