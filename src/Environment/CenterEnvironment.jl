"""
    mutable struct CenterEnvironment <: AbstractEnvironment{3}
        FL::LeftEnvironmentTensor{3}
        O::MPOTensor
        FR::RightEnvironmentTensor{3}
    end

Wrapper type for center bond tensor's generating environment.

Graphic presentation:

     __           __
    |  | --   -- |  |
    |  |    |    |  |
    |FL| —— O —— |FR|
    |  |    |    |  |
    |  | --   -- |  |
     ‾‾           ‾‾

# Constructors
    CenterEnvironment(FL::LeftEnvironmentTensor{3}, O::MPOTensor, FR::RightEnvironmentTensor{3})
"""
mutable struct CenterEnvironment <: AbstractEnvironment{3}
    FL::LeftEnvironmentTensor{3}
    O::MPOTensor
    FR::RightEnvironmentTensor{3}

    function CenterEnvironment(FL::LeftEnvironmentTensor{3}, O::MPOTensor, FR::RightEnvironmentTensor{3})
        (space(FL, 2) == space(O, 1)) || throw(SpaceMismatch("Mismatched auxiliary spaces of left environment tensor and MPO tensor : $(space(FL, 2)) ≠ $(space(O, 1)))"))
        (space(FR, 2) == space(O, 4)) || throw(SpaceMismatch("Mismatched auxiliary spaces of right environment tensor and MPO tensor : $(space(FR, 2)) ≠ $(space(O, 4)))"))
        return new(FL, O, FR)
    end
end
