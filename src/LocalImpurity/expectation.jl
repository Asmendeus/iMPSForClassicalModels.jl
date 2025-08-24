#! Dev memo: Approximate treatment for large W
#! Dev memo: Allow keyword arguments
"""
    expectation(obj::DenseInfiniteMPS{L}, H::SparseUMPO{W, L}, obj′::AdjointInfiniteMPS{L}, M::LocalImpurity{W, L}) -> Number

Calculate 1-site physical quantity expectation.
"""
function expectation(obj::DenseInfiniteMPS{L}, H::SparseUMPO{W, L}, obj′::AdjointInfiniteMPS{L}, M::LocalImpurity{W, L}) where {W, L}

    ψ = canonicalize(obj)
    AL = getAL(ψ)
    AR = getAR(ψ)

    ψ′ = canonicalize(obj′)
    BL = getAL(ψ′)
    BR = getAR(ψ′)

    λ, FL, _ = leftFixedPoint(environment(AL, H.A, BL))
    _, FR, _ = rightFixedPoint(environment(AR, H.A, BR))

    M_val = contract(FL[1], environment(obj.A, M.A, adjoint.(obj′.parent.A)), FR[end])
    Z_val = prod(λ) * contract(environment(FL[1], FR[end]))

    return M_val / Z_val
end

"""
    expectation(obj::DenseInfiniteMPS{L}, H::SparseUMPO{W, L}, obj′::AdjointInfiniteMPS{L},
            M1::LocalImpurity{W, L}, M2::LocalImpurity{W, L}, r::Int64) -> Vector{<:Number}
Calculate 2-site correlation function expectation, return length `r` Vector.
"""
function expectation(obj::DenseInfiniteMPS{L}, H::SparseUMPO{W, L}, obj′::AdjointInfiniteMPS{L},
            M1::LocalImpurity{W, L}, M2::LocalImpurity{W, L}, r::Int64) where {W, L}

    ψ = canonicalize(obj)
    AL = getAL(ψ)
    AR = getAR(ψ)

    ψ′ = canonicalize(obj′)
    BL = getAL(ψ′)
    BR = getAR(ψ′)

    _, FL, _ = leftFixedPoint(environment(AL, H.A, BL))
    _, FR, _ = rightFixedPoint(environment(AR, H.A, BR))

    tmp = pushleft(FL[1], AL, M1.A, BL)

    MM_val = [contract(tmp, environment(obj.A, M2.A, adjoint.(obj′.parent.A)), FR[end]),]

    for _ in 1:r
        tmp = pushleft(tmp, AL, H.A, BL)
        push!(MM_val, contract(tmp, environment(obj.A, M2.A, adjoint.(obj′.parent.A)), FR[end]))
    end

    Z_val = prod(λ) * contract(environment(FL[1], FR[end]))
    ZZ_val = Z_val * (prod(λ) .^ (1:r))

    return MM_val ./ ZZ_val
end
