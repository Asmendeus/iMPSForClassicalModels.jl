"""
     abstract type AbstractTensorWrapper

Wrapper type for classifying different tensors.

Note each concrete subtype must have a field `A::AbstractTensorMap` to save the tensor.
"""
abstract type AbstractTensorWrapper end
