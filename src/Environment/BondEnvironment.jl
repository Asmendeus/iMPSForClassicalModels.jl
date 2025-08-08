"""
    mutable struct BondEnvironment{N} <: AbstractEnvironment{N}
        FL::LeftEnvironmentTensor{N}
        FR::RightEnvironmentTensor{N}
    end

Wrapper type for center bond tensor's generating environment.

Graphic presentation:

     __            __
    |  | --    -- |  |
    |  |          |  |
    |  | ———————— |  |
    |FL|     ⋮    |FR|
    |  | ———————— |  |
    |  |          |  |
    |  | --    -- |  |
     ‾‾            ‾‾

# Constructors
    BondEnvironment{N}(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N})
    BondEnvironment(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N})
"""
mutable struct BondEnvironment{N} <: AbstractEnvironment{N}
    FL::LeftEnvironmentTensor{N}
    FR::RightEnvironmentTensor{N}

    function BondEnvironment{N}(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N}) where N
        aspace_FL = otimes([space(FL, i) for i in 2:N-1]...)
        aspace_FR = otimes([space(FR, i) for i in N-1:-1:2]...)
        aspace_FL == aspace_FR || throw(SpaceMismatch("Mismatched environment tensors' auxiliary spaces: $(aspace_FL) ≠ $(aspace_FR))"))
        new{N}(FL, FR)
    end
    BondEnvironment(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N}) where N = BondEnvironment{N}(FL, FR)
end
