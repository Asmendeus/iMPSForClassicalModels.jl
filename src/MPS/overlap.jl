"""
    overlap(obj::DenseInfiniteMPS{L, T}, obj′::DenseInfiniteMPS{L, T};
            XL₀::AbstractVector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(TransferMatrix(getAL(obj), adjoint.(getAL(obj′)))),
            alg::EigenAlgorithm=Defaults.alg_eig) where {L, T}

Overlap of two iMPSs, i.e., `∏_l λ[l]` solved by
     __                                              __
    |  | -- AL[1] -- ... -- AL[L] --                |  | --
    |FL|    |               |         =  ∏_l λ[l] * |FL|
    |  | -- BL[1] -- ... -- BL[L] --                |  | --
     ‾‾                                              ‾‾
`abs(overlap(obj, obj′)) = 1` means the two are equivalent in the sense of gauge degrees of freedom.
Note: Based on numerical experience, overlap less than 1e-6 indicates that the two iMPS are approximate
"""
function overlap(obj::DenseInfiniteMPS{L, T}, obj′::DenseInfiniteMPS{L, T};
            XL₀::AbstractVector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(TransferMatrix(getAL(obj), adjoint.(getAL(obj′)))),
            alg::EigenAlgorithm=Defaults.alg_eig) where {L, T}
    t = TransferMatrix(getAL(obj), adjoint.(getAL(obj′)))
    λ, _, _ = leftFixedPoint(t, XL₀, alg)
    return prod(λ)
end
