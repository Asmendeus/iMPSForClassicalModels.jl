"""
    expectation(obj::DenseInfiniteMPS{L}, H::SparseGeneralMPO{W, L}, obj′::AdjointInfiniteMPS{L}, M::LocalImpurity{W, L}) -> Number

Calculate 1-site physical quantity expectation.
"""
function expectation(obj::DenseInfiniteMPS{L}, H::SparseGeneralMPO{W, L}, obj′::AdjointInfiniteMPS{L}, M::LocalImpurity{W, L}) where {W, L}

    λ, FL, _ = leftFixedPoint(environment(getAL(obj), H.A, getAL(obj′)))
    _, FR, _ = rightFixedPoint(environment(getAR(obj), H.A, getAR(obj′)))

    M_val = contract(FL[1], environment([getAC(obj)[1], getAR(obj)[2:end]...], M.A, [getAC(obj′)[1], getAR(obj′)[2:end]...]), FR[end])
    Z_val = prod(λ) * contract(getC(obj)[end], environment(FL[1], FR[end]), getC(obj′)[end])

    return M_val / Z_val
end

"""
    expectation(obj::DenseInfiniteMPS{L}, H::SparseGeneralMPO{W, L}, obj′::AdjointInfiniteMPS{L},
            M1::LocalImpurity{W, L}, M2::LocalImpurity{W, L}, r::Int64) -> Vector{<:Number}
Calculate 2-site correlation function expectation, return length `r` Vector.
"""
function expectation(obj::DenseInfiniteMPS{L}, H::SparseGeneralMPO{W, L}, obj′::AdjointInfiniteMPS{L},
            M1::LocalImpurity{W, L}, M2::LocalImpurity{W, L}, r::Int64) where {W, L}

    λ, FL, _ = leftFixedPoint(environment(getAL(obj), H.A, getAL(obj′)))
    _, FR, _ = rightFixedPoint(environment(getAR(obj), H.A, getAR(obj′)))

    tmp = pushleft(FL[1], AL, M1.A, BL)

    MM_val = [contract(tmp, environment([getAC(obj)[1], getAR(obj)[2:end]...], M.A, [getAC(obj′)[1], getAR(obj′)[2:end]...]), FR[end]),]

    for _ in 1:r
        tmp = pushleft(tmp, getAL(obj), H.A, getAL(obj′))
        push!(MM_val, contract(tmp, environment([getAC(obj)[1], getAR(obj)[2:end]...], M.A, [getAC(obj′)[1], getAR(obj′)[2:end]...]), FR[end]))
    end

    Z_val = prod(λ) * contract(getC(obj)[end], environment(FL[1], FR[end]), getC(obj′)[end])
    ZZ_val = Z_val * (prod(λ) .^ (1:r))

    return MM_val ./ ZZ_val
end
