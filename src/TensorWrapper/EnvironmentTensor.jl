"""
    abstract type AbstractEnvironmentTensor{N} <: AbstractTensorWrapper

N-leg environment tensors.
"""
abstract type AbstractEnvironmentTensor{N} <: AbstractTensorWrapper end
numind(::AbstractEnvironmentTensor{N}) where N = N

"""
    struct LeftEnvironmentTensor{N} <: AbstractEnvironmentTensor{N}
        A::AbstractTensorMap
    end

Wrapper type for N-leg left environment tensor.

Convention (' marks codomain):

     __
    |  | -- N
    |  |
    |  | —— N-1
    |FL|    ⋮
    |  | —— 2
    |  |
    |  | -- 1'
     ‾‾

# Constructors
    LeftEnvironmentTensor(::AbstractTensorMap)
    LeftEnvironmentTensor{N}(::AbstractTensorMap)
"""
struct LeftEnvironmentTensor{N} <: AbstractEnvironmentTensor{N}
    A::AbstractTensorMap

    function LeftEnvironmentTensor(A::AbstractTensorMap)
        (numout(A) == 1 && numin(A) > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `LeftEnvironmentTensor` space"))
        N = numind(A)
        return new{N}(A)
    end
    function LeftEnvironmentTensor{N}(A::AbstractTensorMap) where N
        (numout(A) == 1 && numin(A) == N-1 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `LeftEnvironmentTensor{$N}` space"))
        return new{N}(A)
    end
end

"""
    struct RightEnvironmentTensor{N} <: AbstractEnvironmentTensor{N}
        A::AbstractTensorMap
    end

Wrapper type for N-leg right environment tensor.

Convention (' marks codomain):
               __
        1' -- |  |
              |  |
        2' —— |  |
        ⋮     |FR|
    (N-1)' —— |  |
              |  |
         N -- |  |
               ‾‾
# Constructors
    RightEnvironmentTensor(::AbstractTensorMap)
    RightEnvironmentTensor{N}(::AbstractTensorMap)
"""
struct RightEnvironmentTensor{N} <: AbstractEnvironmentTensor{N}
    A::AbstractTensorMap

    function RightEnvironmentTensor(A::AbstractTensorMap)
        (numin(A) == 1 && numout(A) > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `RightEnvironmentTensor` space"))
        N = numind(A)
        return new{N}(A)
    end
    function RightEnvironmentTensor{N}(A::AbstractTensorMap) where N
        (numin(A) == 1 && numout(A) == N-1 > 0) || throw(ArgumentError("The space $(space(A)) does not conform to `RightEnvironmentTensor{$N}` space"))
        return new{N}(A)
    end
end
