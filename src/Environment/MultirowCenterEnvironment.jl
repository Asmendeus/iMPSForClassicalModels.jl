"""
    mutable struct MultirowCenterEnvironment{N} <: AbstractEnvironment{N}
        FL::LeftEnvironmentTensor{N}
        const O::AbstractVector{MPOTensor}
        FR::RightEnvironmentTensor{N}
    end

Wrapper type for center bond tensor's generating environment.

Graphic presentation:

     __            __
    |  | --    -- |  |
    |  |    |     |  |
    |  | —— O1 —— |  |
    |  |    |     |  |
    |FL|    ⋮     |FR|
    |  |    |     |  |
    |  | —— OW —— |  |  W = N - 2
    |  |    |     |  |
    |  | --    -- |  |
     ‾‾            ‾‾

# Constructors
    MultirowCenterEnvironment{N}(FL::LeftEnvironmentTensor{N}, O::AbstractVector{MPOTensor}, FR::RightEnvironmentTensor{N})
    MultirowCenterEnvironment(FL::LeftEnvironmentTensor{N}, O::AbstractVector{MPOTensor}, FR::RightEnvironmentTensor{N})
"""
mutable struct MultirowCenterEnvironment{N} <: AbstractEnvironment{N}
    FL::LeftEnvironmentTensor{N}
    const O::AbstractVector{MPOTensor}
    FR::RightEnvironmentTensor{N}

    function MultirowCenterEnvironment{N}(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N}) where N
        length(O) == N - 2 || throw(ArgumentError("The number of MPO tensors does not match the dimension of environment tensors"))

        aspace_FL = otimes([space(FL, i) for i in 2:N-1]...)
        aspace_OL = otimes([space(O[l], 1) for l in N-2:-1:1]...)
        aspace_FL == aspace_OL || throw(SpaceMismatch("Mismatched auxiliary spaces of left environment tensor and MPO tensors: $(aspace_FL) ≠ $(aspace_OL))"))

        aspace_FR = otimes([space(FR, i) for i in N-1:-1:2]...)
        aspace_OR = otimes([space(O[l], 4) for l in N-2:-1:1]...)
        aspace_FR == aspace_OR || throw(SpaceMismatch("Mismatched auxiliary spaces of right environment tensor and MPO tensors: $(aspace_FR) ≠ $(aspace_OR))"))

        pspace_Odn = [space(O[i], 2) for i in 1:N-3]
        pspace_Oup = [space(O[i], 3) for i in 2:N-2]
        pspace_Odn == pspace_Oup || throw(SpaceMismatch("Mismatched intra-MPO physical spaces of MPO tensors: $(pspace_Oup) ≠ $(pspace_Odn))"))

        return new{N}(FL, O, FR)
    end
    MultirowCenterEnvironment(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N}) where N = MultirowCenterEnvironment{N}(FL, O, FR)
end
