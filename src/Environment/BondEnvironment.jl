"""
    mutable struct BondEnvironment{N} <: AbstractEnvironment{0, N}
        FL::LeftEnvironmentTensor{N}
        FR::RightEnvironmentTensor{N}
    end

Wrapper type for center bond tensor's generating environment.

Graphic presentation:

     __            __
    |  | ——    —— |  |
    |  |          |  |
    |  | ———————— |  |
    |FL|     ⋮    |FR|
    |  | ———————— |  |
    |  |          |  |
    |  | ——    —— |  |
     ‾‾            ‾‾

# Constructors
    BondEnvironment{N}(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N})
    BondEnvironment(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N})
"""
mutable struct BondEnvironment{N} <: AbstractEnvironment{0, N}
    FL::LeftEnvironmentTensor{N}
    FR::RightEnvironmentTensor{N}

    BondEnvironment{N}(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N}) where N = new{N}(FL, FR)
    BondEnvironment(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N}) where N = new{N}(FL, FR)
end
