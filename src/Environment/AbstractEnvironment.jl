"""
     abstract type AbstractEnvironment{N}

Wrapper environment type for solving fixed point equations to
    1. generate environment tensor;
    2. generate center tensor and center bond tensor.

Note `N = 2` corresponds to transfer matrix.
"""
abstract type AbstractEnvironment{N} end
