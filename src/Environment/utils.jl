"""
    Environment(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N}) -> BondEnvironment{N}(FL, FR)
    Environment(FL::LeftEnvironmentTensor{3}, O::MPOTensor, FR::RightEnvironmentTensor{3}) -> CenterEnvironment(FL, O, FR)
    Environment(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}}) -> LeftEnvironment(AL, O, BL)
    Environment(AL::LeftIsometricTensor{R}, O::MPOTensor, BL::AdjointLeftIsometricTensor{R}) -> LeftEnvironment(AL, O, BL)
    Environment(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}}) -> RightEnvironment(AR, O, BR)
    Environment(AR::RightIsometricTensor{R}, O::MPOTensor, BR::AdjointRightIsometricTensor{R}) -> RightEnvironment(AR, O, BR)
    Environment(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N}) -> MultirowCenterEnvironment(FL, O, FR)
    Environment(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}}) -> MultirowLeftEnvironment(AL, O, BL)
    Environment(AL::LeftIsometricTensor{R}, O::AbstractVecOrMat{<:MPOTensor}, BL::AdjointLeftIsometricTensor{R}) -> MultirowLeftEnvironment(AL, O, BL)
    Environment(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}}) -> MultirowRightEnvironment(AR, O, BR)
    Environment(AR::RightIsometricTensor{R}, O::AbstractVecOrMat{<:MPOTensor}, BR::AdjointRightIsometricTensor{R}) -> MultirowRightEnvironment(AR, O, BR)

"""
Environment(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N}) where N = BondEnvironment(FL, FR)

Environment(FL::LeftEnvironmentTensor{3}, O::MPOTensor, FR::RightEnvironmentTensor{3}) = CenterEnvironment(FL, O, FR)

Environment(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}}) where R = LeftEnvironment(AL, O, BL)
Environment(AL::LeftIsometricTensor{R}, O::MPOTensor, BL::AdjointLeftIsometricTensor{R}) where R = LeftEnvironment(AL, O, BL)

Environment(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractVector{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}}) where R = RightEnvironment(AR, O, BR)
Environment(AR::RightIsometricTensor{R}, O::MPOTensor, BR::AdjointRightIsometricTensor{R}) where R = RightEnvironment(AR, O, BR)

Environment(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N}) where N = MultirowCenterEnvironment(FL, O, FR)

Environment(AL::AbstractVector{<:LeftIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor{R}}) where R = MultirowLeftEnvironment(AL, O, BL)
Environment(AL::LeftIsometricTensor{R}, O::AbstractVecOrMat{<:MPOTensor}, BL::AdjointLeftIsometricTensor{R}) where R = MultirowLeftEnvironment(AL, O, BL)

Environment(AR::AbstractVector{<:RightIsometricTensor{R}}, O::AbstractMatrix{<:MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor{R}}) where R = MultirowRightEnvironment(AR, O, BR)
Environment(AR::RightIsometricTensor{R}, O::AbstractVecOrMat{<:MPOTensor}, BR::AdjointRightIsometricTensor{R}) where R = MultirowRightEnvironment(AR, O, BR)
