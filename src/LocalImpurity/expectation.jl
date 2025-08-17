#! Dev memo: Approximate treatment for large W
"""
    expectation(obj::DenseInfiniteMPS{L}, H::SparseInfiniteMPO{W, L}, obj′::AdjointInfiniteMPS{L}, M::LocalImpurity{W, L};
            alg::EigenAlgorithm=Defaults.alg_eig) -> Number

Calculate 1-site physical quantity expectation.
"""
function expectation(obj::DenseInfiniteMPS{L}, H::SparseInfiniteMPO{W, L}, obj′::AdjointInfiniteMPS{L}, M::LocalImpurity{W, L};
            alg::EigenAlgorithm=Defaults.alg_eig, kwargs...) where {W, L}

    AL, AR, _, _ = getAllCanonicalFormTensors(obj; kwargs...)
    _, _, BL, _ = leftFixedPoint(adjoint.(obj′.parent.A), alg)
    _, _, BR, _ = rightFixedPoint(adjoint.(obj′.parent.A), alg)

    _, FL, _ = leftFixedPoint(environment(AL, H.A, BL))
    _, FR, _ = rightFixedPoint(environment(AR, H.A, BR))

    return contract(FL[1], environment(obj.A, M.A, adjoint.(obj′.parent.A)), FR[end])
end

"""
    expectation(obj::DenseInfiniteMPS{L}, H::SparseInfiniteMPO{W, L}, obj′::AdjointInfiniteMPS{L},
            M1::LocalImpurity{W, L}, M2::LocalImpurity{W, L}, r::Int64;
            alg::EigenAlgorithm=Defaults.alg_eig) -> Vector{<:Number}
Calculate 2-site correlation function expectation, return length `r` Vector.
"""
function expectation(obj::DenseInfiniteMPS{L}, H::SparseInfiniteMPO{W, L}, obj′::AdjointInfiniteMPS{L},
            M1::LocalImpurity{W, L}, M2::LocalImpurity{W, L}, r::Int64;
            alg::EigenAlgorithm=Defaults.alg_eig) where {W, L}

    _, _, AL, _ = leftFixedPoint(obj.A, alg)
    _, _, AR, _ = rightFixedPoint(obj.A, alg)
    _, _, BL, _ = leftFixedPoint(adjoint.(obj′.parent.A), alg)
    _, _, BR, _ = rightFixedPoint(adjoint.(obj′.parent.A), alg)

    _, FL, _ = leftFixedPoint(environment(AL, H.A, BL))
    _, FR, _ = rightFixedPoint(environment(AR, H.A, BR))

    tmp = pushleft(FL[1], AL, M1.A, BL)

    MM = [contract(tmp, environment(obj.A, M.A, adjoint.(obj′.parent.A)), FR[end]),]

    for l in 1:L
        tmp = pushleft(tmp, AL, H.A, BL)
        push!(MM, contract(tmp, environment(obj.A, M.A, adjoint.(obj′.parent.A)), FR[end]))
    end

    return MM
end
