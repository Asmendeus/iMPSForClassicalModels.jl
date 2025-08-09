"""
     abstract type AbstractEnvironment{N}

Wrapper environment type for
    1. generating environment tensor;
    2. updating center tensor and center bond tensor;
    3. realizing methods on transfer matrix (transfer matrix is a special `N = 2` environment).
"""
abstract type AbstractEnvironment{N} end
