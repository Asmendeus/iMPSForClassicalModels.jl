abstract type AbstractEnvironmentTensor end

"""
     __
    |  | ——
    |  |
    |FL| ——
    |  |
    |  | ——
     ‾‾
"""
struct LeftEnvironmentTensor <: AbstractEnvironmentTensor
end

"""
        __
    —— |  |
       |  |
    —— |FR|
       |  |
    —— |  |
        ‾‾
"""
struct RightEnvironmentTensor <: AbstractEnvironmentTensor
end