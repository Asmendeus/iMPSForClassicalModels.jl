"""
     abstract type AbstractEnvironment{N}

Wrapper environment type for
    1. generating environment tensor;
    2. updating center tensor and center bond tensor.

Note `N = 2` corresponds to transfer matrix.
"""
abstract type AbstractEnvironment{N} end
