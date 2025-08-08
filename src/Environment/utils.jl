"""
    Environment(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N}) -> BondEnvironment{N}(FL, FR)
    Environment(FL::LeftEnvironmentTensor{3}, O::MPOTensor, FR::RightEnvironmentTensor{3}) -> CenterEnvironment(FL, O, FR)
    Environment(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N}) -> MultirowCenterEnvironment(FL, O, FR)
    Environment(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}}) -> LeftEnvironment(AL, O, BL)
    Environment(AL::LeftIsometricTensor{R}, O::MPOTensor, BL::AdjointLeftIsometricTensor{R}) -> LeftEnvironment(AL, O, BL)
"""
Environment(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N}) where N = BondEnvironment(FL, FR)
Environment(FL::LeftEnvironmentTensor{3}, O::MPOTensor, FR::RightEnvironmentTensor{3}) = CenterEnvironment(FL, O, FR)
Environment(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N}) where N = MultirowCenterEnvironment(FL, O, FR)
Environment(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}}) where R = LeftEnvironment(AL, O, BL)
Environment(AL::LeftIsometricTensor{R}, O::MPOTensor, BL::AdjointLeftIsometricTensor{R}) where R = LeftEnvironment(AL, O, BL)
