"""
    mutable struct CenterEnvironment{N} <: AbstractEnvironment{N}
        FL::LeftEnvironmentTensor{N}
        const O::AbstractVector{MPOTensor}
        FR::RightEnvironmentTensor{N}
    end

Wrapper type for center tensor's generating environment.

Graphic presentation:

     __              __
    |  | --   ----- |  |
    |  |    |       |  |
    |  | —— O[1] —— |  |
    |FL|    ⋮       |FR|
    |  | —— O[W] —— |  |      W = N - 2
    |  |    |       |  |
    |  | --   ----- |  |
     ‾‾              ‾‾

# Constructors
    CenterEnvironment{N}(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N})
    CenterEnvironment(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N})
    CenterEnvironment(FL::LeftEnvironmentTensor{N}, O::MPOTensor, FR::RightEnvironmentTensor{N})
"""
mutable struct CenterEnvironment{N} <: AbstractEnvironment{N}
    FL::LeftEnvironmentTensor{N}
    const O::AbstractVector{MPOTensor}
    FR::RightEnvironmentTensor{N}

    function CenterEnvironment{N}(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N}) where N

        aspace_FL = [domain(FL, i) for i in 1:N-2]
        aspace_OL = [codomain(O[l], 1) for l in 1:N-2]
        (aspace_FL == aspace_OL) || throw(SpaceMismatch("Mismatched auxiliary spaces of left environment tensor and MPO tensor: $(aspace_FL) ≠ $(aspace_OL))"))

        aspace_FR = [codomain(FR, i) for i in N-1:-1:2]
        aspace_OR = [domain(O[l], 2) for l in 1:N-2]
        (aspace_FR == aspace_OR) || throw(SpaceMismatch("Mismatched auxiliary spaces of MPO tensor and right environment tensor: $(aspace_OR) ≠ $(aspace_FR))"))

        return new(FL, O, FR)
    end
    CenterEnvironment(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N}) where N = CenterEnvironment{N}(FL, O, FR)
    CenterEnvironment(FL::LeftEnvironmentTensor{3}, O::MPOTensor, FR::RightEnvironmentTensor{3}) = CenterEnvironment{3}(FL, [O,], FR)
end
