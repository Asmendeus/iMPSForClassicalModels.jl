"""
environment(A::LocalTensor{R}, B::LocalTensor{R}) -> TransferMatrix(A, B)
    environment(A::LocalTensor{R}, O::MPOTensor, B::LocalTensor{R}) -> IsometricEnvironment(A, O, B)
    environment(A::LocalTensor{R}, O::AbstractVector{<:MPOTensor}, B::LocalTensor{R}) -> IsometricEnvironment(A, O, B)
    environment(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N}) -> BondEnvironment(FL, FR)
    environment(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N}) -> CenterEnvironment(FL, O, FR)
    environment(FL::LeftEnvironmentTensor{N}, O::MPOTensor, FR::RightEnvironmentTensor{N}) -> CenterEnvironment(FL, O, FR)
"""
environment(A::LocalTensor{R}, B::AdjointLocalTensor{R}) where R = TransferMatrix(A, B)
environment(A::LocalTensor{R}, O::AbstractVector{<:MPOTensor}, B::AdjointLocalTensor{R}) where R = IsometricEnvironment(A, O, B)
environment(A::LocalTensor{R}, O::MPOTensor, B::AdjointLocalTensor{R}) where R = IsometricEnvironment(A, O, B)
environment(FL::LeftEnvironmentTensor{N}, FR::RightEnvironmentTensor{N}) where N = BondEnvironment(FL, FR)
environment(FL::LeftEnvironmentTensor{N}, O::AbstractVector{<:MPOTensor}, FR::RightEnvironmentTensor{N}) where N = CenterEnvironment(FL, O, FR)
environment(FL::LeftEnvironmentTensor{3}, O::MPOTensor, FR::RightEnvironmentTensor{3}) = CenterEnvironment(FL, O, FR)
