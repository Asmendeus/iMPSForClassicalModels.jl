"""
    contract(FL::LeftEnvironmentTensor{N}, env::ChannelEnvironment{N, L, R}, FR::RightEnvironmentTensor{N}) -> Number
    contract(AC::LocalTensor{R}, env::CenterEnvironment{N}, AC′::AdjointLocalTensor{R}) -> Number
    contract(C::BondTensor, env::BondEnvironment{N}, C′::AdjointBondTensor) -> Number

Contract losed tensor network into a single number
1.
         __                                   __
        |  | -- A[1] ----- ... -- A[L] ----- |  |
        |  |    |                 |          |  |
        |  | —— O[1, 1] —— ... —— O[1, L] —— |  |
        |  |    |                 |          |  |
        |FL|    ⋮                  ⋮          |FR|
        |  |    |                 |          |  |
        |  | —— O[W, 1] —— ... —— O[W, L] —— |  |
        |  |    |                 |          |  |
        |  | -- B[1] ----- ... -- B[L] ----- |  |
         ‾‾                                   ‾‾
2.
         __              __
        |  | -- A ----- |  |
        |  |    |       |  |
        |  | —— O[1] —— |  |
        |  |    |       |  |
        |FL|    ⋮        |FR|
        |  |    |       |  |
        |  | —— O[W] —— |  |
        |  |    |       |  |
        |  | -- B ----- |  |
         ‾‾              ‾‾
3.
         __            __
        |  | -- CA -- |  |
        |  |          |  |
        |  | ———————— |  |
        |  |          |  |
        |FL|    ⋮      |FR|
        |  |          |  |
        |  | ———————— |  |
        |  |          |  |
        |  | -- CB -- |  |
         ‾‾            ‾‾
"""
function contract(FL::LeftEnvironmentTensor{N}, env::ChannelEnvironment{N, L, R}, FR::RightEnvironmentTensor{N}) where {N, L, R}
end
function contract(AC::LocalTensor{R}, env::CenterEnvironment{N}, AC′::AdjointLocalTensor{R}) where {N, R}
end
function contract(C::BondTensor, env::BondEnvironment{N}, C′::AdjointBondTensor) where N
end
