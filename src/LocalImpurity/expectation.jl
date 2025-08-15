"""
    expectation(ψ::DenseInfiniteMPS{L}, H::SparseMPO{W, L}, ψ′::AdjointInfiniteMPS{L}, M::LocalImpurity{W, L};
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Defaults.alg_eig()) -> Number

Calculate 1-site physical quantity expectation.
"""
function expectation(ψ::DenseInfiniteMPS{L}, H::SparseMPO{W, L}, ψ′::AdjointInfiniteMPS{L}, M::LocalImpurity{W, L};
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Defaults.alg_eig()) where {W, L}

    _, _, AL, _ = leftFixedPoint(ψ.A, alg)
    _, _, AR, _ = rightFixedPoint(ψ.A, alg)
    _, _, BL, _ = leftFixedPoint(adjoint.(ψ′.parent.A), alg)
    _, _, BR, _ = rightFixedPoint(adjoint.(ψ′.parent.A), alg)

    _, FL, _ = leftFixedPoint(environment(AL, H.A, BL))
    _, FR, _ = rightFixedPoint(environment(AR, H.A, BR))

    return contract(FL[1], environment(ψ.A, M.A, adjoint.(ψ′.parent.A)), FR[end])
end

"""
    expectation(ψ::DenseInfiniteMPS{L}, H::SparseMPO{W, L}, ψ′::AdjointInfiniteMPS{L},
            M1::LocalImpurity{W, L}, M2::LocalImpurity{W, L}, r::Int64;
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Defaults.alg_eig()) -> Vector{<:Number}
Calculate 2-site correlation function expectation, return length `r` Vector.
"""
function expectation(ψ::DenseInfiniteMPS{L}, H::SparseMPO{W, L}, ψ′::AdjointInfiniteMPS{L},
            M1::LocalImpurity{W, L}, M2::LocalImpurity{W, L}, r::Int64;
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Defaults.alg_eig()) where {W, L}

    _, _, AL, _ = leftFixedPoint(ψ.A, alg)
    _, _, AR, _ = rightFixedPoint(ψ.A, alg)
    _, _, BL, _ = leftFixedPoint(adjoint.(ψ′.parent.A), alg)
    _, _, BR, _ = rightFixedPoint(adjoint.(ψ′.parent.A), alg)

    _, FL, _ = leftFixedPoint(environment(AL, H.A, BL))
    _, FR, _ = rightFixedPoint(environment(AR, H.A, BR))

    tmp = pushleft(FL[1], AL, M1.A, BL)

    MM = [contract(tmp, environment(ψ.A, M.A, adjoint.(ψ′.parent.A)), FR[end]),]

    for l in 1:L
        tmp = pushleft(tmp, AL, H.A, BL)
        push!(MM, contract(tmp, environment(ψ.A, M.A, adjoint.(ψ′.parent.A)), FR[end]))
    end

    return MM
end
