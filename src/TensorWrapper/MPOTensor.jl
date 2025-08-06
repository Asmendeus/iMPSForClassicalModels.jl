"""
    abstract type AbstractMPOTensor <: AbstractTensorWrapper

Elements of rank-4 MPO.
"""
abstract type AbstractMPOTensor <: AbstractTensorWrapper end

"""
    struct MPOTensor <: AbstractMPOTensor
        A::AbstractTensorMap
    end

Wrapper type for rank-4 MPO local tensors.

Convention (' marks codomain):

          2'
          |
    1' —— O —— 4
          |
          3

# Constructor(s)
    MPOTensor(::AbstractTensorMap)
"""
struct MPOTensor <: AbstractMPOTensor
    A::AbstractTensorMap

    function MPOTensor(A::AbstractTensorMap)
        (length(codomain(A)) == 2 && length(domain(A)) == 2) || throw(ArgumentError("The space $(space(A)) does not conform to `MPOTensor` space"))
        return new(A)
    end
end
